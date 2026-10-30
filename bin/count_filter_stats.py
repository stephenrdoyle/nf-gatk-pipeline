#!/usr/bin/env python3
"""
Count variants before and after filtering from a filtered VCF,
and merge with a filter threshold TSV to produce a combined summary.

Usage:
    count_filter_stats.py <raw_vcf> <filtered_vcf> <variant_type> [<threshold_tsv>]

    threshold_tsv is optional — for hard/custom filters where thresholds
    are fixed, pass a pre-built TSV. For optimised filters, pass the TSV
    produced by optimised_filter.py. If omitted, threshold columns are NA.

Output:
    {variant_type}_counts.tsv   — written to current directory
"""
import sys, subprocess, csv, os

raw_vcf      = sys.argv[1]
filtered_vcf = sys.argv[2]
variant_type = sys.argv[3]
threshold_tsv = sys.argv[4] if len(sys.argv) > 4 else None

def count_variants(vcf):
    """Count PASS and total variants using bcftools stats — memory efficient."""
    result = subprocess.run(
        ['bcftools', 'stats', vcf],
        capture_output=True, text=True
    )
    total   = 0
    passing = 0
    for line in result.stdout.splitlines():
        # SN section: summary numbers
        if line.startswith('SN') and 'number of records' in line:
            total = int(line.strip().split()[-1])
        if line.startswith('SN') and 'number of no-ALTs' in line:
            # subtract ref-only records if any
            pass
    # Count PASS using shell pipe to wc -l — avoids loading records into Python memory
    pass_result = subprocess.run(
        'bcftools view -H -f PASS,. ' + vcf + ' | wc -l',
        shell=True, capture_output=True, text=True
    )
    passing = int(pass_result.stdout.strip()) if pass_result.stdout.strip() else 0
    flagged = total - passing
    return total, passing, flagged

total_before, pass_before, _ = count_variants(raw_vcf)
total_after,  pass_after, flagged_after = count_variants(filtered_vcf)

print(f"{variant_type}: {total_before} variants before, {pass_after} PASS after filtering ({flagged_after} flagged)")

# Load threshold data if provided
thresholds = {}
if threshold_tsv and os.path.exists(threshold_tsv):
    with open(threshold_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ann = row['annotation']
            thresholds[ann] = row

# Write combined counts TSV
out_file = f"{variant_type}_counts.tsv"
with open(out_file, 'w') as f:
    header = [
        'variant_type', 'annotation', 'operator', 'threshold',
        'n_total_before', 'n_pass_before',
        'n_total_after',  'n_pass_after', 'n_flagged',
        'pct_flagged_by_filter', 'pct_flagged_overall',
        'mean', 'median'
    ]
    f.write('\t'.join(header) + '\n')

    if thresholds:
        for ann, row in thresholds.items():
            # Hard/custom TSVs carry 'NA' here (only optimised mode counts
            # per-annotation flags), so pass the value through untouched
            n_flag_ann = row.get('n_flagged') or 'NA'
            pct_ann    = row.get('pct_flagged', 'NA')
            pct_overall = round(100 * flagged_after / total_before, 2) if total_before else 0
            f.write('\t'.join([
                variant_type,
                ann,
                row.get('operator', 'NA'),
                row.get('threshold', 'NA'),
                str(total_before),
                str(pass_before),
                str(total_after),
                str(pass_after),
                str(n_flag_ann),
                str(pct_ann),
                str(pct_overall),
                row.get('mean', 'NA'),
                row.get('median', 'NA'),
            ]) + '\n')
    else:
        # No per-annotation data (hard/custom without threshold TSV)
        pct_overall = round(100 * flagged_after / total_before, 2) if total_before else 0
        f.write('\t'.join([
            variant_type, 'all', 'NA', 'NA',
            str(total_before), str(pass_before),
            str(total_after),  str(pass_after),
            str(flagged_after),
            'NA', str(pct_overall),
            'NA', 'NA'
        ]) + '\n')

print(f"Counts written: {out_file}")
