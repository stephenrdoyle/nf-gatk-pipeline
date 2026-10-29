#!/usr/bin/env python3
"""
Compute percentile-based filter thresholds from GATK annotation TSV,
apply filters via VariantFiltration, and generate an HTML report.

Usage: optimised_filter.py <annotations.tsv> <variant_type> <lo_pct> <annotations> <lo_bad> <hi_bad>
"""
import csv, math, json, sys, os, subprocess, statistics

tsv_file     = sys.argv[1]
variant_type = sys.argv[2]
lo_pct       = int(sys.argv[3])
hi_pct       = 100 - lo_pct
annotations  = sys.argv[4].split(',')
lo_bad       = sys.argv[5].split(',')
hi_bad       = sys.argv[6].split(',')

def percentile(data, pct):
    if not data: return None
    data = sorted(data)
    k    = (len(data) - 1) * pct / 100
    f, c = math.floor(k), math.ceil(k)
    return data[f] if f == c else data[f] * (c - k) + data[c] * (k - f)

# Read annotation values
values = {a: [] for a in annotations}
with open(tsv_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        for ann in annotations:
            try:
                v = row.get(ann, '')
                if v and v not in ('.', 'NA', 'na'):
                    values[ann].append(float(v))
            except ValueError:
                pass

# Derive thresholds
thresholds = {}
for ann in annotations:
    if ann in lo_bad:
        t = percentile(values[ann], lo_pct)
        if t is not None:
            thresholds[ann] = ('lt', t)
    elif ann in hi_bad:
        t = percentile(values[ann], hi_pct)
        if t is not None:
            thresholds[ann] = ('gt', t)

# Write thresholds JSON
with open('thresholds.json', 'w') as f:
    json.dump({k: {'direction': v[0], 'value': v[1]} for k, v in thresholds.items()}, f, indent=2)

# Write a shell script fragment that can be sourced to build the GATK command.
# Using an array avoids all quoting and word-splitting issues with negative numbers.
with open('filter_args.sh', 'w') as f:
    f.write('FILTER_ARGS=()\n')
    for ann, (direction, value) in thresholds.items():
        op   = '<' if direction == 'lt' else '>'
        # Prefix with variant type so SNP and indel filters are distinguishable in the VCF header
        vtype_prefix = variant_type.upper().rstrip('S')  # snps -> SNP, indels -> INDEL
        op_str = 'min' if direction == 'lt' else 'max'   # min/max clearer than LT/GT
        name = f"{vtype_prefix}_{ann}_{op_str}{value:.2f}".replace('.', '_').replace('-', 'neg')
        f.write(f'FILTER_ARGS+=(--filter-expression "{ann} {op} {value:.4f}")\n')
        f.write(f'FILTER_ARGS+=(--filter-name "{name}")\n')

print(f"Thresholds computed for: {list(thresholds.keys())}")

# ── per-type TSV summary ──────────────────────────────────────────────────────
# Written as {variant_type}_filter_summary.tsv for downstream merging.
# Columns: variant_type, annotation, operator, threshold, n_total,
#          n_flagged, pct_flagged, mean, median
import statistics as _stats

tsv_out = f"{variant_type}_filter_summary.tsv"
with open(tsv_out, 'w') as f:
    f.write("variant_type\tannotation\toperator\tthreshold\tn_total\tn_flagged\tpct_flagged\tmean\tmedian\n")
    for ann, (direction, threshold) in thresholds.items():
        vals   = values.get(ann, [])
        op     = '<' if direction == 'lt' else '>'
        n_tot  = len(vals)
        n_flag = sum(1 for v in vals if (v < threshold if direction == 'lt' else v > threshold))
        pct_f  = round(100 * n_flag / n_tot, 2) if n_tot else 0
        mean   = round(_stats.mean(vals), 4)   if vals else 'NA'
        median = round(_stats.median(vals), 4) if vals else 'NA'
        f.write(f"{variant_type}\t{ann}\t{op}\t{threshold:.4f}\t{n_tot}\t{n_flag}\t{pct_f}\t{mean}\t{median}\n")
print(f"Filter summary written: {tsv_out}")

# ── HTML report ───────────────────────────────────────────────────────────────
html = []
html.append(f"""<!DOCTYPE html><html><head>
<meta charset="utf-8">
<title>Filter report: {variant_type}</title>
<style>
  body  {{ font-family: sans-serif; max-width: 1000px; margin: 40px auto; padding: 0 20px; color: #333; }}
  h1    {{ color: #2c3e50; }}
  h2    {{ color: #34495e; border-bottom: 1px solid #eee; padding-bottom: 6px; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 30px; }}
  th,td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
  th    {{ background: #f5f5f5; }}
  .hist {{ font-family: monospace; font-size: 12px; white-space: pre; line-height: 1.4; }}
  .fail {{ color: #e74c3c; }}
</style></head><body>""")

html.append(f"<h1>Variant Filter Report: {variant_type.upper()}</h1>")
html.append(f"<p>Thresholds derived at the {lo_pct}th / {hi_pct}th percentiles of each annotation distribution.</p>")

html.append("<h2>Derived Thresholds</h2>")
html.append("<table><tr><th>Annotation</th><th>Operator</th><th>Threshold</th>"
            "<th>N values</th><th>Mean</th><th>Median</th><th>% flagged</th></tr>")
for ann, (direction, threshold) in thresholds.items():
    vals   = values.get(ann, [])
    op     = '<' if direction == 'lt' else '>'
    n_fail = sum(1 for v in vals if (v < threshold if direction == 'lt' else v > threshold))
    pct_f  = 100 * n_fail / len(vals) if vals else 0
    mean   = f"{statistics.mean(vals):.3f}"   if vals else 'NA'
    median = f"{statistics.median(vals):.3f}" if vals else 'NA'
    html.append(f"<tr><td>{ann}</td><td>{op}</td><td>{threshold:.4f}</td>"
                f"<td>{len(vals)}</td><td>{mean}</td><td>{median}</td><td>{pct_f:.1f}%</td></tr>")
html.append("</table>")

html.append("<h2>Annotation Distributions</h2>")
for ann, (direction, threshold) in thresholds.items():
    vals = values.get(ann, [])
    if not vals:
        continue
    lo, hi   = min(vals), max(vals)
    n_bins   = 40
    bin_size = (hi - lo) / n_bins if hi != lo else 1
    bins     = [0] * n_bins
    for v in vals:
        bins[min(int((v - lo) / bin_size), n_bins - 1)] += 1
    max_count = max(bins) or 1
    bar_width = 40
    op = '<' if direction == 'lt' else '>'
    html.append(f"<h3>{ann}  (filter: {ann} {op} {threshold:.4f})</h3>")
    html.append("<div class='hist'>")
    html.append(f"{'Value':>10}  {'Count':>6}  {'Bar':<40}")
    for i, count in enumerate(bins):
        bin_lo = lo + i * bin_size
        bin_hi = bin_lo + bin_size
        bar    = '█' * int(bar_width * count / max_count)
        in_filter = (bin_hi <= threshold) if direction == 'lt' else (bin_lo >= threshold)
        marker = ' &lt;-- FILTERED' if in_filter else ''
        cls    = " class='fail'" if in_filter else ''
        html.append(f"<span{cls}>{bin_lo:10.3f}  {count:6d}  {bar:<40}{marker}</span>")
    html.append("</div>")

html.append("</body></html>")

report_name = f"{variant_type}_filter_report.html"
with open(report_name, 'w') as f:
    f.write('\n'.join(html))
print(f"Report written: {report_name}")
