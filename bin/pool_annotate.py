#!/usr/bin/env python3
"""
Post-process a cohort VCF of pooled samples.

HaplotypeCaller and GenotypeGVCFs call every sample as a diploid
individual. A pool is a mix of many individuals whose allele frequency can
take any value, so those diploid genotypes are not meaningful. This script:
  1. sets every GT to missing (./.), so no downstream tool treats the
     diploid calls as real genotypes
  2. adds FORMAT/AF: the frequency of each ALT allele in the pool from AD,
     AD[alt] / sum(AD), one value per ALT allele; '.' when FORMAT/DP is
     below min_dp or AD is missing
  3. records in the header that the samples are pools
All other FORMAT fields (AD, DP, ...) are left unchanged.

Usage:
    pool_annotate.py <input.vcf[.gz]> <output.vcf | -> [min_dp]

Output is uncompressed VCF ('-' for stdout); the caller BGZF-compresses it.
Progress messages go to stderr.
"""
import sys, gzip

in_vcf  = sys.argv[1]
out_vcf = sys.argv[2]
min_dp  = int(sys.argv[3]) if len(sys.argv) > 3 else 0

def allele_freqs(ad_field, n_alt):
    """Per-ALT allele frequencies from AD, or None if AD is unusable."""
    try:
        counts = [int(x) for x in ad_field.split(',')]
    except ValueError:
        return None
    total = sum(counts)
    if len(counts) != n_alt + 1 or total == 0:
        return None
    return [c / total for c in counts[1:]]

def depth_ok(dp_field):
    if min_dp <= 0:
        return True
    try:
        return int(dp_field) >= min_dp
    except ValueError:
        return False

def missing_gt(gt):
    """Missing genotype with the same ploidy: 0/1 -> ./., 1 -> ."""
    return '/'.join('.' for _ in gt.replace('|', '/').split('/'))

opener = gzip.open if in_vcf.endswith('.gz') else open
n_sites = n_masked = n_af = 0
has_af_header = False

with opener(in_vcf, 'rt') as fin, (sys.stdout if out_vcf == '-' else open(out_vcf, 'w')) as fout:
    for line in fin:
        # ── header ────────────────────────────────────────────────────────────
        if line.startswith('##'):
            if line.startswith('##FORMAT=<ID=AF,'):
                has_af_header = True
            fout.write(line)
            continue
        if line.startswith('#CHROM'):
            if not has_af_header:
                fout.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description='
                           '"Frequency of each ALT allele in the pool, AD[alt] / sum(AD); '
                           f'missing when DP < {min_dp}">\n')
            fout.write('##poolAnnotate=<Description="Pooled samples: GT set to missing '
                       'because diploid genotype calls are not meaningful for pools; '
                       f'use FORMAT/AD or FORMAT/AF",MinDP={min_dp}>\n')
            for sample in line.rstrip('\n').split('\t')[9:]:
                fout.write(f'##SAMPLE=<ID={sample},Pooled=true>\n')
            fout.write(line)
            continue

        # ── records ───────────────────────────────────────────────────────────
        n_sites += 1
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 10:
            fout.write(line)
            continue

        n_alt    = 0 if fields[4] == '.' else len(fields[4].split(','))
        fmt_keys = fields[8].split(':')
        if 'AF' not in fmt_keys:
            fmt_keys.append('AF')
            fields[8] = ':'.join(fmt_keys)
        gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else None
        ad_idx = fmt_keys.index('AD') if 'AD' in fmt_keys else None
        dp_idx = fmt_keys.index('DP') if 'DP' in fmt_keys else None
        af_idx = fmt_keys.index('AF')

        for s in range(9, len(fields)):
            vals = fields[s].split(':')
            while len(vals) < len(fmt_keys):
                vals.append('.')

            if gt_idx is not None:
                masked = missing_gt(vals[gt_idx])
                if masked != vals[gt_idx]:
                    vals[gt_idx] = masked
                    n_masked += 1

            freqs = None
            if n_alt and ad_idx is not None and (dp_idx is None or depth_ok(vals[dp_idx])):
                freqs = allele_freqs(vals[ad_idx], n_alt)
            if freqs is not None:
                vals[af_idx] = ','.join(f'{f:.4f}' for f in freqs)
                n_af += 1
            else:
                vals[af_idx] = '.'

            fields[s] = ':'.join(vals)

        fout.write('\t'.join(fields) + '\n')

print(f'Processed {n_sites} sites: {n_masked} genotypes set to missing, '
      f'AF added for {n_af} sample-sites', file=sys.stderr)
