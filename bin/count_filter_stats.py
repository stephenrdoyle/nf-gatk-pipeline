#!/usr/bin/env python3
"""
Summarise site filtering from the filtered SNP and indel VCFs.

Counts come from each VCF's FILTER column, and each filter's expression from
its ##FILTER header line (VariantFiltration writes the expression as the
description), so the summary is the same whichever filter mode produced the
VCF: hard, optimised, custom or DRAGEN.

Output (stdout, TSV), one row per variant type and filter:
  id            unique row key (MultiQC tables are keyed on the first column)
  variant_type  SNPs | Indels (incl. mixed)
  filter        "(all sites)", PASS, or a filter name
  expression    the filter's expression from the VCF header
  n_sites       sites with this FILTER value (PASS: sites passing everything)
  pct_sites     n_sites as % of all sites of this type
  n_only        sites failed by this filter and no other
Filters that flagged no sites are not listed.

Usage:
    count_filter_stats.py snps=<snps.vcf.gz> indels=<indels.vcf.gz> > filter_summary.tsv
"""
import re
import subprocess
import sys
from collections import Counter

TYPE_NAMES = {'snps': 'SNPs', 'indels': 'Indels (incl. mixed)'}

def filter_expressions(vcf):
    """FILTER ID -> description (the filter expression) from the header."""
    header = subprocess.run(['bcftools', 'view', '--header-only', vcf],
                            check=True, capture_output=True, text=True).stdout
    found = {}
    for line in header.splitlines():
        m = re.match(r'##FILTER=<ID=([^,>]+)(?:,Description="(.*)")?>', line)
        if m:
            found[m.group(1)] = m.group(2) or ''
    return found

def filter_combinations(vcf):
    """Count of each FILTER value ('PASS', 'A', 'A;B', ...) across all sites."""
    counts = Counter()
    proc = subprocess.Popen(['bcftools', 'query', '--format', '%FILTER\n', vcf],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        counts[line.rstrip('\n')] += 1
    if proc.wait() != 0:
        # a silent bcftools failure would otherwise report zero sites
        sys.exit(f'ERROR: bcftools query failed on {vcf} (exit {proc.returncode})')
    return counts

def rows_for(vtype, vcf):
    name = TYPE_NAMES.get(vtype, vtype)
    expressions = filter_expressions(vcf)
    combos = filter_combinations(vcf)
    total = sum(combos.values())
    pct = lambda n: f'{100 * n / total:.2f}' if total else '0.00'

    n_pass = combos.get('PASS', 0) + combos.get('.', 0)
    per_filter, only = Counter(), Counter()
    for value, n in combos.items():
        if value in ('PASS', '.'):
            continue
        names = value.split(';')
        for f in names:
            per_filter[f] += n
        if len(names) == 1:
            only[names[0]] += n

    yield [f'{name}: all sites', name, '(all sites)', '', total, pct(total), '']
    yield [f'{name}: PASS', name, 'PASS', 'passes every filter', n_pass, pct(n_pass), '']
    for f, n in sorted(per_filter.items(), key=lambda kv: -kv[1]):
        yield [f'{name}: {f}', name, f, expressions.get(f, ''), n, pct(n), only[f]]

out = sys.stdout
out.write('\t'.join(['id', 'variant_type', 'filter', 'expression',
                     'n_sites', 'pct_sites', 'n_only']) + '\n')
for arg in sys.argv[1:]:
    vtype, vcf = arg.split('=', 1)
    for row in rows_for(vtype, vcf):
        out.write('\t'.join(str(x) for x in row) + '\n')
