#!/usr/bin/env python3
"""
Post-process a VCF for pooled samples:
  1. Add AF = AD[1] / (AD[0] + AD[1]) to FORMAT per genotype
  2. Add ##SAMPLE pooled flag to header
  3. Mask GT to ./. where allele balance is inconsistent with
     a diploid call (not near 0%, 50%, or 100%)

Usage:
    pool_annotate.py <input.vcf.gz> <output.vcf.gz> <sample_id> \
                     [het_lo] [het_hi] [hom_alt_lo] [min_dp]

    het_lo     lower bound for het allele balance  [default: 0.35]
    het_hi     upper bound for het allele balance  [default: 0.65]
    hom_alt_lo lower bound for hom-alt allele balance [default: 0.90]
    
    Genotypes with allele balance outside:
        0.0 - (1 - hom_alt_lo)  → hom-ref zone
        het_lo - het_hi          → het zone
        hom_alt_lo - 1.0         → hom-alt zone
    are masked to ./.
"""
import sys, gzip, re

in_vcf      = sys.argv[1]
out_vcf     = sys.argv[2]
sample_id   = sys.argv[3]
het_lo      = float(sys.argv[4]) if len(sys.argv) > 4 else 0.35
het_hi      = float(sys.argv[5]) if len(sys.argv) > 5 else 0.65
hom_alt_lo  = float(sys.argv[6]) if len(sys.argv) > 6 else 0.90
min_dp      = int(sys.argv[7])   if len(sys.argv) > 7 else 0
hom_ref_hi  = 1.0 - hom_alt_lo  # upper bound for hom-ref zone

def calc_af(ad_field):
    """Calculate allele frequency from AD field. Returns None if missing."""
    try:
        parts = ad_field.split(',')
        if len(parts) < 2:
            return None
        ref_count = int(parts[0])
        alt_count = int(parts[1])
        total = ref_count + alt_count
        return alt_count / total if total > 0 else None
    except (ValueError, ZeroDivisionError):
        return None

def mask_gt(gt, af):
    """
    Return masked GT ./. if AF is inconsistent with the called diploid genotype.
    Phased genotypes (|) are also masked.
    """
    if af is None:
        return gt
    # Normalise phased to unphased for comparison
    gt_norm = gt.replace('|', '/')
    if gt_norm in ('./.', '.'):
        return gt
    # Check if AF falls in an expected diploid zone
    in_hom_ref  = af <= hom_ref_hi
    in_het      = het_lo <= af <= het_hi
    in_hom_alt  = af >= hom_alt_lo
    if in_hom_ref or in_het or in_hom_alt:
        return gt   # consistent with diploid — keep
    # AF is in an ambiguous zone — mask
    return './.'

n_total  = 0
n_masked = 0
n_af_added = 0

opener  = gzip.open if in_vcf.endswith('.gz') else open
wopener = open  # always write uncompressed — caller handles BGZF compression

header_sample_added = False

with opener(in_vcf, 'rt') as fin, wopener(out_vcf, 'wt') as fout:
    for line in fin:
        # ── header lines ──────────────────────────────────────────────────────
        if line.startswith('#'):
            if line.startswith('##FORMAT=<ID=GT') and not header_sample_added:
                # Insert AF FORMAT definition and SAMPLE tag before first FORMAT line
                fout.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency calculated from AD (alt / total)">\n')
                fout.write(f'##SAMPLE=<ID={sample_id},Pooled=true,HetLo={het_lo},HetHi={het_hi},HomAltLo={hom_alt_lo}>\n')
                header_sample_added = True
            fout.write(line)
            continue

        # ── variant lines ─────────────────────────────────────────────────────
        n_total += 1
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 10:
            fout.write(line)
            continue

        fmt    = fields[8]
        fmt_keys = fmt.split(':')

        # Add AF to FORMAT if not already present
        if 'AF' not in fmt_keys:
            fields[8] = fmt + ':AF'
            fmt_keys.append('AF')

        ad_idx = fmt_keys.index('AD') if 'AD' in fmt_keys else None
        dp_idx = fmt_keys.index('DP') if 'DP' in fmt_keys else None
        gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else None
        af_idx = fmt_keys.index('AF')

        # Process each sample (columns 9+)
        new_samples = []
        for sample_field in fields[9:]:
            sample_vals = sample_field.split(':')

            # Pad sample values if shorter than format keys
            while len(sample_vals) < len(fmt_keys):
                sample_vals.append('.')

            # Calculate AF from AD — only if depth meets the minimum threshold
            af = None
            dp_ok = True
            if dp_idx is not None and dp_idx < len(sample_vals):
                try:
                    dp_val = sample_vals[dp_idx]
                    if dp_val in ('.', ''):
                        dp_ok = False
                    elif int(dp_val) < min_dp:
                        dp_ok = False
                except ValueError:
                    dp_ok = False

            if dp_ok and ad_idx is not None and ad_idx < len(sample_vals):
                af = calc_af(sample_vals[ad_idx])

            # Set AF value
            if af is not None:
                sample_vals[af_idx] = f'{af:.4f}'
                n_af_added += 1
            else:
                sample_vals[af_idx] = '.'

            # Mask ambiguous GT
            if gt_idx is not None and gt_idx < len(sample_vals):
                original_gt = sample_vals[gt_idx]
                sample_vals[gt_idx] = mask_gt(original_gt, af)
                if sample_vals[gt_idx] != original_gt:
                    n_masked += 1

            new_samples.append(':'.join(sample_vals))

        fields[9:] = new_samples
        fout.write('\t'.join(fields) + '\n')

print(f"Processed {n_total} variants")
print(f"AF added: {n_af_added}")
print(f"GT masked: {n_masked} ({100*n_masked/n_total:.1f}% of variants)" if n_total else "GT masked: 0")
