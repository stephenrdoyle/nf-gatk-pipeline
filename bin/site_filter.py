#!/usr/bin/env python3
"""
Site filters for the cohort VCF, decided per site from a table of its
annotations, plus a diagnostic report on how well the thresholds fit.

SITE_FILTER (modules/site_filter.nf) runs it in two steps:

  apply   reads the per-site table written by bcftools query (one line per
          record: CHROM POS REF ALT F_MISSING then the INFO tags in
          --columns), decides each record's FILTER and writes a sites-only VCF
          of those decisions to stdout; bcftools annotate copies them into the
          cohort VCF. Also writes:
            filter_summary.tsv        sites per filter (MultiQC, QC report)
            expected_filters.tsv      count of each FILTER value, for 'verify'
            filter_diagnostics.tsv    per annotation: distribution, threshold,
                                      fraction flagged, warning
            filter_diagnostics.html   the same, with histograms
  verify  compares the FILTER values bcftools actually wrote with
          expected_filters.tsv, so a mis-matched annotation can't pass silently

Variant types follow the pipeline's split rule: a record is a SNP if REF
and every ALT (ignoring the spanning-deletion '*') are single bases;
anything else (indels, mixed sites) is filtered with the indel rules.

A missing annotation passes, as in GATK VariantFiltration: an expression on
a value the site doesn't have does not flag it.

Standard library only.
"""

from __future__ import annotations

import argparse
import html
import math
import re
import subprocess
import sys
from collections import Counter

# GATK hard-filter thresholds: (annotation, operator, value, filter name)
HARD = {
    "snp": [("QD", "<", 2.0, "SNP_QD2"), ("MQ", "<", 40.0, "SNP_MQ40"), ("FS", ">", 60.0, "SNP_FS60"),
            ("SOR", ">", 3.0, "SNP_SOR3"), ("MQRankSum", "<", -12.5, "SNP_MQRankSum-12.5"),
            ("ReadPosRankSum", "<", -8.0, "SNP_ReadPosRankSum-8")],
    "indel": [("QD", "<", 2.0, "INDEL_QD2"), ("FS", ">", 200.0, "INDEL_FS200"), ("SOR", ">", 10.0, "INDEL_SOR10"),
              ("ReadPosRankSum", "<", -20.0, "INDEL_ReadPosRankSum-20")],
}
EXCESS_HET_DEFAULT = 54.69          # GATK germline short-variant guidance
DEPTH_MULTIPLES = (1.5, 2.0, 3.0, 5.0)
MISSING_LEVELS = (0.1, 0.2, 0.5)
WARN_PCT = 10.0                     # one filter flagging more than this: warn
TYPE_NAMES = {"snp": "SNPs", "indel": "Indels (incl. mixed)"}
OPS = {"<": lambda a, b: a < b, "<=": lambda a, b: a <= b, ">": lambda a, b: a > b,
       ">=": lambda a, b: a >= b, "==": lambda a, b: a == b, "!=": lambda a, b: a != b}
EXPR_RE = re.compile(r"^\s*(?:INFO/)?([A-Za-z_][A-Za-z0-9_.]*)\s*(<=|>=|==|!=|<|>)\s*"
                     r"(-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*$")
FIXED_COLS = ["CHROM", "POS", "REF", "ALT", "F_MISSING"]


# ── rules ────────────────────────────────────────────────────────────────────
def custom_name(vtype: str, ann: str, op: str, value: str) -> str:
    """'QD < 2.0' -> SNP_QD_lt_2.0: distinct for distinct expressions."""
    word = {"<": "lt", "<=": "le", ">": "gt", ">=": "ge", "==": "eq", "!=": "ne"}[op]
    return f"{'SNP' if vtype == 'snp' else 'INDEL'}_{ann}_{word}_{value.replace('-', 'neg')}"


def parse_expressions(text: str, vtype: str) -> list:
    """'QD < 2.0,FS > 60.0' -> [(ann, op, value, name)]; exits on anything else."""
    rules = []
    for expr in (e for e in (text or "").split(",") if e.strip()):
        m = EXPR_RE.match(expr)
        if not m:
            sys.exit(f"ERROR: can't parse filter expression '{expr.strip()}': "
                     "expected '<INFO tag> <op> <number>', e.g. 'QD < 2.0' (op: < <= > >= == !=)")
        ann, op, value = m.groups()
        rules.append((ann, op, float(value), custom_name(vtype, ann, op, value)))
    return rules


def classify(ref: str, alt: str) -> str:
    alts = [a for a in alt.split(",") if a not in ("*", ".")]
    if alts and len(ref) == 1 and all(len(a) == 1 for a in alts):
        return "snp"
    return "indel"


def num(v: str):
    if v in (".", ""):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def parse_depth(text: str | None):
    """None | ('x', 2.0) for '2x' | ('abs', 15000.0)"""
    if not text or str(text).lower() in ("false", "null", "none"):
        return None
    t = str(text).strip().lower()
    if t.endswith("x"):
        return ("x", float(t[:-1]))
    return ("abs", float(t))


# ── reading ──────────────────────────────────────────────────────────────────
def read_rows(path: str):
    if path.endswith(".gz"):
        proc = subprocess.Popen(["gzip", "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        stream = proc.stdout
    else:
        proc, stream = None, open(path)
    try:
        for line in stream:
            if line.strip():
                yield line.rstrip("\n").split("\t")
    finally:
        stream.close()
    if proc and proc.wait() != 0:
        sys.exit(f"ERROR: could not read {path}")


def median_of(counts: Counter):
    total = sum(counts.values())
    if not total:
        return None
    seen = 0
    for v in sorted(counts):
        seen += counts[v]
        if seen >= total / 2:
            return v


def quantile(counts: Counter, q: float):
    total = sum(counts.values())
    if not total:
        return None
    seen = 0
    for v in sorted(counts):
        seen += counts[v]
        if seen >= q * total:
            return v


# ── apply ────────────────────────────────────────────────────────────────────
def apply(a) -> None:
    cols = FIXED_COLS + a.columns.split(",")
    idx = {c: i for i, c in enumerate(cols)}
    rules = {"snp": list(HARD["snp"]), "indel": list(HARD["indel"])} if a.mode == "hard" else \
            {"snp": parse_expressions(a.snp_filters, "snp"), "indel": parse_expressions(a.indel_filters, "indel")}
    for vtype in rules:
        for ann, *_ in rules[vtype]:
            if ann not in idx:
                sys.exit(f"ERROR: annotation {ann} is not in --columns ({a.columns})")
    depth = parse_depth(a.max_depth)
    max_missing = float(a.max_missing) if a.max_missing not in (None, "", "null", "false") else None
    excess_het = None
    if a.excess_het not in (None, "", "false", "null"):
        excess_het = EXCESS_HET_DEFAULT if a.excess_het == "true" else float(a.excess_het)

    # pass 1: distributions (depth median for --max_depth Nx, and the diagnostics)
    dist = {t: {c: Counter() for c in cols[4:]} for t in ("snp", "indel")}
    n_type = Counter()
    for f in read_rows(a.sites):
        t = classify(f[2], f[3])
        n_type[t] += 1
        for c in cols[4:]:
            v = num(f[idx[c]])
            if v is not None:
                dist[t][c][int(v) if c == "DP" else round(v, 3 if c == "F_MISSING" else 2)] += 1
    if a.expect_records is not None and sum(n_type.values()) != a.expect_records:
        sys.exit(f"ERROR: site table has {sum(n_type.values())} records, the VCF has {a.expect_records}")
    dp_median = {t: median_of(dist[t].get("DP", Counter())) for t in dist}
    depth_max = {}
    if depth:
        if "DP" not in idx:
            sys.exit("ERROR: --max-depth needs DP in --columns")
        for t in dist:
            depth_max[t] = depth[1] * dp_median[t] if depth[0] == "x" and dp_median[t] else \
                           depth[1] if depth[0] == "abs" else None

    # every rule, with its FILTER name and VCF header description
    defs = {}   # name -> description
    for t in rules:
        for ann, op, val, name in rules[t]:
            defs[name] = f"{ann} {op} {val:g}"
    if depth:
        desc = (f"INFO/DP > {depth[1]:g} x median site depth (SNPs {fmt_n(depth_max.get('snp'))}, "
                f"indels {fmt_n(depth_max.get('indel'))})" if depth[0] == "x" else f"INFO/DP > {depth[1]:g}")
        defs["HighDepth"] = desc
    if max_missing is not None:
        defs["HighMissing"] = f"fraction of samples with no genotype call > {max_missing:g}"
    if excess_het is not None:
        defs["ExcessHet"] = f"ExcessHet > {excess_het:g}"

    # pass 2: decide, write the decisions VCF, tally
    combos = {t: Counter() for t in rules}
    out = sys.stdout
    out.write("##fileformat=VCFv4.2\n")
    for name, desc in defs.items():
        out.write(f'##FILTER=<ID={name},Description="{desc}">\n')
    with open(a.contigs) as fh:
        for line in fh:
            if line.startswith("##contig"):
                out.write(line)
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for f in read_rows(a.sites):
        t = classify(f[2], f[3])
        failed = []
        for ann, op, val, name in rules[t]:
            v = num(f[idx[ann]])
            if v is not None and OPS[op](v, val):
                failed.append(name)
        if depth and depth_max.get(t) is not None:
            v = num(f[idx["DP"]])
            if v is not None and v > depth_max[t]:
                failed.append("HighDepth")
        if max_missing is not None:
            v = num(f[idx["F_MISSING"]])
            if v is not None and v > max_missing:
                failed.append("HighMissing")
        if excess_het is not None:
            v = num(f[idx["ExcessHet"]]) if "ExcessHet" in idx else None
            if v is not None and v > excess_het:
                failed.append("ExcessHet")
        value = ";".join(failed) if failed else "PASS"
        combos[t][value] += 1
        out.write(f"{f[0]}\t{f[1]}\t.\t{f[2]}\t{f[3]}\t.\t{value}\t.\n")
    out.flush()

    write_summary(a.summary, combos, defs)
    with open(a.expected, "w") as fh:
        total = Counter()
        for c in combos.values():
            total.update(c)
        for value, n in sorted(total.items()):
            fh.write(f"{value}\t{n}\n")
    diag = diagnostics(rules, dist, n_type, combos, dp_median, depth, depth_max, max_missing, excess_het)
    write_diag_tsv(a.diag_tsv, diag)
    with open(a.diag_html, "w") as fh:
        fh.write(diag_html(diag, rules, dist, dp_median, depth, depth_max, max_missing, excess_het, a.mode))
    print(f"site filter: {sum(n_type.values()):,} records; PASS "
          f"{sum(c.get('PASS', 0) for c in combos.values()):,}", file=sys.stderr)


def fmt_n(v):
    return "–" if v is None else f"{v:,.0f}"


def write_summary(path, combos, defs):
    """filter_summary.tsv, same columns as before: id variant_type filter
    expression n_sites pct_sites n_only."""
    with open(path, "w") as fh:
        fh.write("\t".join(["id", "variant_type", "filter", "expression", "n_sites", "pct_sites", "n_only"]) + "\n")
        for t in ("snp", "indel"):
            name = TYPE_NAMES[t]
            c = combos[t]
            total = sum(c.values())
            pct = lambda n: f"{100 * n / total:.2f}" if total else "0.00"
            per, only = Counter(), Counter()
            for value, n in c.items():
                if value == "PASS":
                    continue
                names = value.split(";")
                for x in names:
                    per[x] += n
                if len(names) == 1:
                    only[names[0]] += n
            rows = [[f"{name}: all sites", name, "(all sites)", "", total, pct(total), ""],
                    [f"{name}: PASS", name, "PASS", "passes every filter", c.get("PASS", 0), pct(c.get("PASS", 0)), ""]]
            rows += [[f"{name}: {x}", name, x, defs.get(x, ""), n, pct(n), only[x]]
                     for x, n in sorted(per.items(), key=lambda kv: -kv[1])]
            for r in rows:
                fh.write("\t".join(str(v) for v in r) + "\n")


# ── diagnostics ──────────────────────────────────────────────────────────────
def flagged_by(counts: Counter, op: str, val: float) -> int:
    return sum(n for v, n in counts.items() if OPS[op](v, val))


def diagnostics(rules, dist, n_type, combos, dp_median, depth, depth_max, max_missing, excess_het):
    """One row per (type, check): what it flags and whether that looks off."""
    rows = []
    for t in ("snp", "indel"):
        total = n_type[t]
        if not total:
            continue
        only = Counter()
        for value, n in combos[t].items():
            if value != "PASS" and ";" not in value:
                only[value] += n
        checks = [(ann, op, val, name, True) for ann, op, val, name in rules[t]]
        if "DP" in dist[t] and dp_median[t]:
            for k in DEPTH_MULTIPLES:
                applied = bool(depth and depth[0] == "x" and abs(depth[1] - k) < 1e-9)
                checks.append(("DP", ">", k * dp_median[t], "HighDepth" if applied else None, applied, f"{k:g}x median"))
            if depth and depth[0] == "abs":
                checks.append(("DP", ">", depth[1], "HighDepth", True))
        for lvl in MISSING_LEVELS:
            applied = max_missing is not None and abs(max_missing - lvl) < 1e-9
            checks.append(("F_MISSING", ">", lvl, "HighMissing" if applied else None, applied))
        if max_missing is not None and max_missing not in MISSING_LEVELS:
            checks.append(("F_MISSING", ">", max_missing, "HighMissing", True))
        if "ExcessHet" in dist[t]:
            eh = excess_het if excess_het is not None else EXCESS_HET_DEFAULT
            checks.append(("ExcessHet", ">", eh, "ExcessHet" if excess_het is not None else None, excess_het is not None))
        for c in checks:
            ann, op, val, name, applied = c[:5]
            label = c[5] if len(c) > 5 else ""
            counts = dist[t].get(ann, Counter())
            n_val = sum(counts.values())
            n_flag = flagged_by(counts, op, val)
            pct = 100 * n_flag / total
            warn = ""
            if applied and pct > WARN_PCT:
                warn = (f"flags {pct:.1f}% of {TYPE_NAMES[t].lower()}: check that this threshold suits this cohort "
                        "(see the histogram)")
            rows.append({
                "variant_type": TYPE_NAMES[t], "annotation": ann, "rule": f"{ann} {op} {val:g}" + (f" ({label})" if label else ""),
                "applied": "yes" if applied else "no", "filter": name or "", "n_sites": total,
                "n_with_value": n_val, "n_flagged": n_flag, "pct_flagged": round(pct, 2),
                "n_only": only.get(name, 0) if applied and name else "",
                "p01": quantile(counts, .01), "p05": quantile(counts, .05), "p50": quantile(counts, .5),
                "p95": quantile(counts, .95), "p99": quantile(counts, .99), "warning": warn,
            })
    return rows


def write_diag_tsv(path, rows):
    cols = ["variant_type", "annotation", "rule", "applied", "filter", "n_sites", "n_with_value", "n_flagged",
            "pct_flagged", "n_only", "p01", "p05", "p50", "p95", "p99", "warning"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")


def svg_hist(counts: Counter, marks: list, xlab: str, w=440, h=200, bins=50) -> str:
    """Histogram over the 0.5-99.5% range (stretched to show the marks);
    values outside are piled into the end bins, which are labelled."""
    if not counts:
        return "<p class='note'>No values.</p>"
    lo, hi = quantile(counts, .005), quantile(counts, .995)
    for m, _, _ in marks:
        lo, hi = min(lo, m), max(hi, m)
    if hi == lo:
        hi = lo + 1
    width = (hi - lo) / bins
    hist = [0] * bins
    for v, n in counts.items():
        hist[max(0, min(int((v - lo) / width), bins - 1))] += n
    top = max(math.log10(c + 1) for c in hist) or 1
    ml, mr, mt, mb = 40, 8, 8, 34
    sx = lambda x: ml + (x - lo) / (hi - lo) * (w - ml - mr)
    bw = (w - ml - mr) / bins
    out = [f"<svg viewBox='0 0 {w} {h}' width='100%' style='max-width:{w}px' role='img'>",
           f"<rect x='{ml}' y='{mt}' width='{w-ml-mr}' height='{h-mt-mb}' class='frame'/>"]
    for i, c in enumerate(hist):
        bh = math.log10(c + 1) / top * (h - mt - mb)
        out.append(f"<rect x='{ml+i*bw:.1f}' y='{h-mb-bh:.1f}' width='{max(bw-0.5, .5):.1f}' height='{bh:.1f}' fill='#4a7ea8'>"
                   f"<title>{lo+i*width:.4g} to {lo+(i+1)*width:.4g}: {c:,}</title></rect>")
    for m, label, applied in marks:
        cls = "mark" if applied else "mark2"
        out.append(f"<line x1='{sx(m):.1f}' x2='{sx(m):.1f}' y1='{mt}' y2='{h-mb}' class='{cls}'/>"
                   f"<text x='{sx(m)+2:.1f}' y='{mt+10}' class='tick'>{html.escape(label)}</text>")
    for i in range(5):
        t = lo + i * (hi - lo) / 4
        out.append(f"<text x='{sx(t):.1f}' y='{h-mb+13}' class='tick' text-anchor='middle'>{f'{t:,.0f}' if abs(t) >= 1000 else f'{t:.3g}'}</text>")
    out.append(f"<text x='{(w+ml)/2:.0f}' y='{h-4}' class='axis' text-anchor='middle'>{html.escape(xlab)} (log count)</text></svg>")
    return "".join(out)


def diag_html(rows, rules, dist, dp_median, depth, depth_max, max_missing, excess_het, mode) -> str:
    e = html.escape
    warns = [r for r in rows if r["warning"]]
    p = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
         f"<title>Site filter diagnostics</title><style>{CSS}</style></head><body>",
         "<h1>Site filter diagnostics</h1>",
         f"<p class='lede'>Filter mode <b>{e(mode)}</b>. Each histogram shows one annotation's distribution in this cohort, "
         "with applied thresholds as solid red lines and options not applied as dashed grey lines. Thresholds are fixed, "
         "so this page is where to check that they suit this cohort: a threshold inside the bulk of a distribution, "
         f"or one filter flagging more than {WARN_PCT:g}% of sites, deserves a look.</p>"]
    p.append("<h2>Checks</h2><ul class='checks'>")
    if warns:
        for r in warns:
            p.append(f"<li class='bad'><span class='m'>!</span> <b>{e(r['variant_type'])}, {e(r['rule'])}</b>: {e(r['warning'])}</li>")
    else:
        p.append(f"<li class='ok'><span class='m'>✓</span> No applied filter flags more than {WARN_PCT:g}% of sites on its own count.</li>")
    p.append("</ul>")
    for t in ("snp", "indel"):
        trows = [r for r in rows if r["variant_type"] == TYPE_NAMES[t]]
        if not trows:
            continue
        p.append(f"<h2>{e(TYPE_NAMES[t])}</h2>")
        p.append("<div class='scroll'><table><thead><tr><th>rule</th><th>applied</th><th>% flagged</th><th>flagged</th>"
                 "<th>by this filter alone</th><th>sites with a value</th><th>1%</th><th>5%</th><th>median</th><th>95%</th><th>99%</th></tr></thead><tbody>")
        for r in trows:
            cls = " class='warnrow'" if r["warning"] else ""
            p.append(f"<tr{cls}><td>{e(r['rule'])}</td><td>{e(r['applied'])}</td><td>{r['pct_flagged']:.2f}</td>"
                     f"<td>{r['n_flagged']:,}</td><td>{'' if r['n_only'] == '' else format(r['n_only'], ',')}</td>"
                     f"<td>{r['n_with_value']:,}</td>"
                     + "".join(f"<td>{'' if r[q] is None else format(r[q], 'g')}</td>" for q in ("p01", "p05", "p50", "p95", "p99"))
                     + "</tr>")
        p.append("</tbody></table></div><div class='grid'>")
        anns = []
        for ann, *_ in rules[t]:
            if ann not in anns:
                anns.append(ann)
        anns += [x for x in ("DP", "F_MISSING", "ExcessHet") if x in dist[t] and x not in anns]
        for ann in anns:
            marks = [(val, f"{op}{val:g}", True) for a2, op, val, _ in rules[t] if a2 == ann]
            if ann == "DP" and dp_median[t]:
                for k in DEPTH_MULTIPLES:
                    applied = bool(depth and depth[0] == "x" and abs(depth[1] - k) < 1e-9)
                    marks.append((k * dp_median[t], f"{k:g}x", applied))
                if depth and depth[0] == "abs":
                    marks.append((depth[1], f">{depth[1]:g}", True))
            if ann == "F_MISSING":
                marks += [(lvl, f"{lvl:g}", max_missing is not None and abs(max_missing - lvl) < 1e-9) for lvl in MISSING_LEVELS]
            if ann == "ExcessHet":
                marks.append((excess_het or EXCESS_HET_DEFAULT, f"{excess_het or EXCESS_HET_DEFAULT:g}", excess_het is not None))
            label = {"DP": f"INFO/DP (median {fmt_n(dp_median[t])})", "F_MISSING": "fraction of samples with no call"}.get(ann, ann)
            p.append(f"<div><h3>{e(ann)}</h3>{svg_hist(dist[t][ann], marks, label)}</div>")
        p.append("</div>")
    p.append("<h2>Notes</h2><ul class='methods'>"
             "<li>Sites are flagged, never removed: <code>bcftools view -f PASS</code> keeps the passing ones.</li>"
             "<li>A site missing an annotation passes that filter, as in GATK VariantFiltration.</li>"
             "<li>Depth: INFO/DP is summed over all samples, so the median scales with cohort size and coverage; a ceiling "
             "relative to it catches collapsed repeats and paralogs. Rows for every multiple are shown whether or not one "
             "is applied (<code>--max_depth</code>).</li>"
             "<li>Missingness is measured after the genotype filter. How much to allow depends on the analysis, so it is "
             "off unless <code>--max_missing</code> is set; check per-sample missingness first, since a few poor samples "
             "make every site look incomplete.</li>"
             "<li>ExcessHet flags sites where too many samples are heterozygous (paralogs). Don't apply it to F1 crosses, "
             "where every site that differs between the parents is heterozygous in all offspring.</li></ul></body></html>")
    return "\n".join(p)


CSS = """
:root{--fg:#1d2329;--muted:#5f6b76;--line:#d9dee3;--bg:#fff;--ok:#1a7f37;--bad:#c62828;--warnbg:#fff4e5}
@media (prefers-color-scheme:dark){:root{--fg:#e6e9ec;--muted:#9aa5af;--line:#3a424a;--bg:#14181c;--warnbg:#3a2a12}}
body{font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--fg);background:var(--bg);max-width:1000px;margin:24px auto;padding:0 16px}
h1{font-size:22px}h2{font-size:17px;margin-top:28px;border-top:1px solid var(--line);padding-top:8px}h3{font-size:13px;margin:8px 0 2px}
.lede,.note{color:var(--muted)}table{border-collapse:collapse;font-size:13px}th,td{padding:3px 8px;border-bottom:1px solid var(--line);text-align:left}
.scroll{overflow-x:auto}.warnrow{background:var(--warnbg)}ul.checks{list-style:none;padding:0}.m{display:inline-block;width:1.2em;font-weight:700}
li.ok .m{color:var(--ok)}li.bad .m{color:var(--bad)}.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:12px}
svg text{fill:var(--fg)}.tick{font-size:10px;fill:var(--muted)}.axis{font-size:11px}.frame{fill:none;stroke:var(--line)}
.mark{stroke:var(--bad);stroke-width:1.5}.mark2{stroke:var(--muted);stroke-dasharray:3 3}code{font-size:12px}
"""


# ── verify ───────────────────────────────────────────────────────────────────
def verify(a) -> None:
    norm = lambda v: "PASS" if v in ("PASS", ".") else ";".join(sorted(v.split(";")))
    expected, observed = Counter(), Counter()
    with open(a.expected) as fh:
        for line in fh:
            value, n = line.rstrip("\n").split("\t")
            expected[norm(value)] += int(n)
    with open(a.observed) as fh:
        for line in fh:
            if line.strip():
                observed[norm(line.strip())] += 1
    if expected != observed:
        diff = sorted(set(expected) | set(observed), key=lambda v: -abs(expected[v] - observed[v]))[:5]
        sys.exit("ERROR: FILTER values in the written VCF don't match the decisions "
                 "(value: expected / written): " + "; ".join(f"{v}: {expected[v]} / {observed[v]}" for v in diff))
    print(f"verified: {sum(observed.values()):,} records, FILTER values as decided", file=sys.stderr)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("apply")
    s.add_argument("--sites", required=True)
    s.add_argument("--columns", required=True, help="INFO tags after the fixed columns, comma-separated")
    s.add_argument("--contigs", required=True, help="file with the input VCF's ##contig lines")
    s.add_argument("--mode", choices=["hard", "custom"], required=True)
    s.add_argument("--snp-filters", dest="snp_filters", default="")
    s.add_argument("--indel-filters", dest="indel_filters", default="")
    s.add_argument("--max-depth", dest="max_depth", default=None)
    s.add_argument("--max-missing", dest="max_missing", default=None)
    s.add_argument("--excess-het", dest="excess_het", default=None)
    s.add_argument("--expect-records", dest="expect_records", type=int, default=None)
    s.add_argument("--summary", default="filter_summary.tsv")
    s.add_argument("--expected", default="expected_filters.tsv")
    s.add_argument("--diag-tsv", dest="diag_tsv", default="filter_diagnostics.tsv")
    s.add_argument("--diag-html", dest="diag_html", default="filter_diagnostics.html")
    v = sub.add_parser("verify")
    v.add_argument("expected")
    v.add_argument("observed")
    a = ap.parse_args(argv)
    apply(a) if a.cmd == "apply" else verify(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
