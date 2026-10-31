#!/usr/bin/env python3
"""
Cohort QC report: one self-contained HTML file per run.

Self-contained on purpose: no CDN links, web fonts or JavaScript libraries,
so it opens anywhere. Charts are inline SVG; hovering a point shows the
sample name (an SVG <title>, no script needed).

The report is ordered so that what would invalidate the cohort VCF comes
first: missing samples or genome pieces, duplicated or swapped samples, sex
mismatches and per-sample outliers. The PCA, relatedness and variant
summaries follow.

Nothing here re-analyses the VCF. Every number is read from files the
pipeline already wrote (bcftools stats, missingness, filter summary,
SNPRelate output, software versions), so the report cannot disagree with
them.

Standard library only.

Usage:
    qc_report.py --inputs DIR --samples samples.tsv --meta meta.json --out report.html
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import statistics
import sys

# ── thresholds (documented in the report's Methods section) ─────────────────
ROBUST_Z = 4.0           # |x - median| / (1.4826 * MAD) above this is an outlier
MIN_SAMPLES_FOR_OUTLIERS = 5
KING_DUPLICATE = 0.354   # KING kinship: duplicate / identical twin
KING_FIRST = 0.177       # first degree (parent-offspring, full sibs)
KING_SECOND = 0.0884
KING_THIRD = 0.0442
SEX_MALE_MAX = 0.25      # X/autosome heterozygosity ratio: <= male (one X)
SEX_FEMALE_MIN = 0.5     #                                   >= female
MITO_HET_MAX = 0.2       # heterozygous fraction of non-ref mito genotypes
DENSITY_MAX_CONTIGS = 40 # density plot: contigs with the most sites, in genome order

# Okabe-Ito: distinguishable under common colour vision deficiencies
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
           "#56B4E9", "#F0E442", "#000000"]
GREY = "#9aa0a6"


# ── parsing ──────────────────────────────────────────────────────────────────
def parse_bcftools_stats(path: str) -> dict:
    """SN counts, overall Ts/Tv and per-sample PSC rows, columns by name."""
    out = {"SN": {}, "TSTV": None, "PSC": {}}
    if not path or not os.path.exists(path):
        return out
    psc_cols = None
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if line.startswith("SN\t"):
                out["SN"][parts[2].rstrip(":")] = _num(parts[3])
            elif line.startswith("TSTV\t") and out["TSTV"] is None:
                out["TSTV"] = _num(parts[4])
            elif line.startswith("# PSC"):
                psc_cols = [h.split("]", 1)[-1] for h in parts]
            elif line.startswith("PSC\t") and psc_cols:
                row = dict(zip(psc_cols, parts))
                out["PSC"][row["sample"]] = {k: _num(v) for k, v in row.items()
                                             if k not in ("# PSC", "id", "sample")}
    return out


def read_tsv(path: str) -> list[dict]:
    if not path or not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_text(path: str) -> str:
    if not os.path.exists(path):
        return ""
    with open(path) as f:
        return f.read()


def _num(v):
    try:
        f = float(v)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return v


# ── per-sample metrics ───────────────────────────────────────────────────────
def het_rate(psc: dict) -> float | None:
    """Heterozygous fraction of called genotypes at variant sites."""
    called = psc.get("nRefHom", 0) + psc.get("nNonRefHom", 0) + psc.get("nHets", 0)
    return psc.get("nHets", 0) / called if called else None


def robust_outliers(values: dict, side: str = "both") -> dict:
    """{sample: robust z} for samples beyond ROBUST_Z (median/MAD based)."""
    vals = {k: v for k, v in values.items() if isinstance(v, (int, float))}
    if len(vals) < MIN_SAMPLES_FOR_OUTLIERS:
        return {}
    med = statistics.median(vals.values())
    mad = statistics.median(abs(v - med) for v in vals.values()) * 1.4826
    if mad == 0:
        return {}
    out = {}
    for k, v in vals.items():
        z = (v - med) / mad
        if (side in ("both", "high") and z > ROBUST_Z) or (side in ("both", "low") and z < -ROBUST_Z):
            out[k] = z
    return out


def infer_sex(ratio: float | None) -> str:
    if ratio is None:
        return "?"
    if ratio <= SEX_MALE_MAX:
        return "M"
    if ratio >= SEX_FEMALE_MIN:
        return "F"
    return "ambiguous"


def norm_sex(label: str) -> str:
    s = (label or "").strip().lower()
    if s in ("m", "male", "1"):
        return "M"
    if s in ("f", "female", "2"):
        return "F"
    return ""


def king_degree(k: float) -> str:
    if k > KING_DUPLICATE:
        return "duplicate"
    if k > KING_FIRST:
        return "first degree"
    if k > KING_SECOND:
        return "second degree"
    if k > KING_THIRD:
        return "third degree"
    return "unrelated"


# ── formatting ───────────────────────────────────────────────────────────────
def esc(x) -> str:
    return html.escape(str(x), quote=True)


def fmt(v, digits=3) -> str:
    if v is None or v == "":
        return "–"
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, float):
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        return f"{v:.{digits}g}"
    return esc(v)


def status(ok: bool | None, text: str) -> str:
    cls, mark = {True: ("ok", "✓"), False: ("bad", "✗"), None: ("warn", "!")}[ok]
    return f'<li class="{cls}"><span class="mark">{mark}</span> {text}</li>'


def colour_map(keys) -> dict:
    keys = sorted({k for k in keys if k})
    return {k: PALETTE[i % len(PALETTE)] for i, k in enumerate(keys)}


# ── SVG charts ───────────────────────────────────────────────────────────────
def svg_scatter(points: list, xlab: str, ylab: str, colours: dict,
                label: set = frozenset(), w: int = 460, h: int = 380) -> str:
    """points: (x, y, sample, colour_key). Hover shows the sample name."""
    pts = [p for p in points if isinstance(p[0], (int, float)) and isinstance(p[1], (int, float))]
    if not pts:
        return '<p class="note">No data to plot.</p>'
    ml, mr, mt, mb = 56, 14, 12, 44
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    x0, x1 = _pad(min(xs), max(xs))
    y0, y1 = _pad(min(ys), max(ys))
    sx = lambda x: ml + (x - x0) / (x1 - x0) * (w - ml - mr)
    sy = lambda y: h - mb - (y - y0) / (y1 - y0) * (h - mt - mb)
    out = [f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" role="img">',
           f'<rect x="{ml}" y="{mt}" width="{w-ml-mr}" height="{h-mt-mb}" class="frame"/>']
    for t in _ticks(x0, x1):
        out.append(f'<text x="{sx(t):.1f}" y="{h-mb+16}" class="tick" text-anchor="middle">{_tick(t, x0, x1)}</text>')
    for t in _ticks(y0, y1):
        out.append(f'<text x="{ml-6}" y="{sy(t)+4:.1f}" class="tick" text-anchor="end">{_tick(t, y0, y1)}</text>')
    out.append(f'<text x="{(ml + w - mr) / 2}" y="{h-8}" class="axis" text-anchor="middle">{esc(xlab)}</text>')
    out.append(f'<text x="14" y="{(mt + h - mb) / 2}" class="axis" text-anchor="middle" '
               f'transform="rotate(-90 14 {(mt + h - mb) / 2})">{esc(ylab)}</text>')
    for x, y, name, key in pts:
        c = colours.get(key, GREY)
        out.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{c}" fill-opacity="0.8">'
                   f'<title>{esc(name)}{" (" + esc(key) + ")" if key else ""}</title></circle>')
        if name in label:
            out.append(f'<text x="{sx(x)+6:.1f}" y="{sy(y)-6:.1f}" class="lab">{esc(name)}</text>')
    out.append("</svg>")
    return "".join(out)


def svg_hist(values: list, xlab: str, bins: int = 30, marks: tuple = (),
             w: int = 460, h: int = 220, log_y: bool = False) -> str:
    vals = [v for v in values if isinstance(v, (int, float))]
    if not vals:
        return '<p class="note">No data to plot.</p>'
    ml, mr, mt, mb = 48, 10, 10, 40
    lo, hi = min(vals + [m for m, _ in marks]), max(vals + [m for m, _ in marks])
    if hi == lo:
        hi = lo + 1
    width = (hi - lo) / bins
    counts = [0] * bins
    for v in vals:
        counts[min(int((v - lo) / width), bins - 1)] += 1
    tr = (lambda c: math.log10(c + 1)) if log_y else (lambda c: c)
    top = max(tr(c) for c in counts) or 1
    sx = lambda x: ml + (x - lo) / (hi - lo) * (w - ml - mr)
    out = [f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" role="img">',
           f'<rect x="{ml}" y="{mt}" width="{w-ml-mr}" height="{h-mt-mb}" class="frame"/>']
    bw = (w - ml - mr) / bins
    for i, c in enumerate(counts):
        bh = tr(c) / top * (h - mt - mb)
        out.append(f'<rect x="{ml + i*bw:.1f}" y="{h-mb-bh:.1f}" width="{max(bw-1, 1):.1f}" '
                   f'height="{bh:.1f}" fill="#4a7ea8"><title>{c:,}</title></rect>')
    for m, name in marks:
        out.append(f'<line x1="{sx(m):.1f}" x2="{sx(m):.1f}" y1="{mt}" y2="{h-mb}" class="mark"/>'
                   f'<text x="{sx(m)+3:.1f}" y="{mt+11}" class="tick">{esc(name)}</text>')
    for t in _ticks(lo, hi):
        out.append(f'<text x="{sx(t):.1f}" y="{h-mb+15}" class="tick" text-anchor="middle">{_tick(t, lo, hi)}</text>')
    out.append(f'<text x="{(ml + w - mr) / 2}" y="{h-6}" class="axis" text-anchor="middle">'
               f'{esc(xlab)}{" (log count)" if log_y else ""}</text></svg>')
    return "".join(out)


def svg_density(rows: list, w: int = 920, row_h: int = 34) -> str:
    """Variant sites per Mb along each contig, one strip per contig. Only the
    DENSITY_MAX_CONTIGS contigs with most sites are drawn (in genome order),
    so a fragmented assembly doesn't give thousands of strips."""
    by = {}
    for r in rows:
        by.setdefault(r["contig"], []).append((int(r["mb"]), int(r["n"])))
    if not by:
        return '<p class="note">No density data.</p>'
    dropped = 0
    if len(by) > DENSITY_MAX_CONTIGS:
        keep = set(sorted(by, key=lambda c: -sum(n for _, n in by[c]))[:DENSITY_MAX_CONTIGS])
        dropped = len(by) - len(keep)
        by = {c: v for c, v in by.items() if c in keep}
    longest = max(max(b for b, _ in v) for v in by.values()) + 1
    top = max(n for v in by.values() for _, n in v) or 1
    lw = 190
    bw = (w - lw - 10) / longest
    h = row_h * len(by) + 22
    out = [f'<svg viewBox="0 0 {w} {h}" width="100%" role="img">']
    for i, (contig, bins) in enumerate(by.items()):
        y = i * row_h + 4
        out.append(f'<text x="{lw-6}" y="{y + row_h/2 + 4:.1f}" class="tick" text-anchor="end">{esc(contig[:28])}</text>')
        span = max(b for b, _ in bins) + 1
        out.append(f'<rect x="{lw}" y="{y}" width="{span*bw:.1f}" height="{row_h-6}" class="frame"/>')
        for b, n in bins:
            bh = n / top * (row_h - 6)
            out.append(f'<rect x="{lw + b*bw:.1f}" y="{y + row_h - 6 - bh:.1f}" width="{max(bw, 0.6):.2f}" '
                       f'height="{bh:.1f}" fill="#4a7ea8"><title>{esc(contig)} {b}-{b+1} Mb: {n:,} sites</title></rect>')
    out.append(f'<text x="{lw}" y="{h-4}" class="tick">0 Mb</text>'
               f'<text x="{w-10}" y="{h-4}" class="tick" text-anchor="end">{longest} Mb</text></svg>')
    if dropped:
        out.append(f'<p class="note">{DENSITY_MAX_CONTIGS} contigs with the most sites shown; '
                   f'{dropped:,} smaller contigs not drawn.</p>')
    return "".join(out)


def _pad(lo, hi):
    if hi == lo:
        return lo - 1, hi + 1
    d = (hi - lo) * 0.06
    return lo - d, hi + d


def _ticks(lo, hi, n=5):
    step = (hi - lo) / (n - 1)
    return [lo + i * step for i in range(n)]


def _tick(t, lo=None, hi=None):
    """Tick label with enough decimals to tell neighbouring ticks apart."""
    if lo is not None and hi is not None and hi > lo:
        step = (hi - lo) / 4
        decimals = max(0, 1 - math.floor(math.log10(step)))
        return f"{t:,.{min(decimals, 6)}f}"
    return f"{t:.2g}"


def table(headers: list, rows: list, cls: str = "") -> str:
    out = [f'<table class="{cls}"><thead><tr>' + "".join(f"<th>{esc(h)}</th>" for h in headers) + "</tr></thead><tbody>"]
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


# ── report ───────────────────────────────────────────────────────────────────
def build(inputs: str, samples_tsv: str, meta: dict) -> str:
    p = lambda name: os.path.join(inputs, name)
    nuclear = parse_bcftools_stats(p("stats_nuclear.txt"))
    sexstats = parse_bcftools_stats(p("stats_sex.txt")) if meta.get("sex_chrom") else None
    mitostats = parse_bcftools_stats(p("stats_mito.txt")) if meta.get("mito_chrom") else None
    before = parse_bcftools_stats(p("before_filter_bcftools_stats.txt"))
    after = parse_bcftools_stats(p("after_filter_bcftools_stats.txt"))
    miss_after = {r["sample"]: r for r in read_tsv(p("after_filter_missingness.tsv"))}
    filt = read_tsv(p("filter_summary.tsv"))
    kin = read_tsv(p("kinship.tsv"))
    pca_n = read_tsv(p("pca_nuclear.tsv"))
    eig_n = read_tsv(p("pca_nuclear_eigen.tsv"))
    pca_m = read_tsv(p("pca_mito.tsv"))
    eig_m = read_tsv(p("pca_mito_eigen.tsv"))
    density = read_tsv(p("density.tsv"))
    snp_counts = {r["set"]: int(r["n_snps"]) for r in read_tsv(p("snp_counts.tsv"))}
    notes = [l for l in read_text(p("notes.txt")).splitlines() if l.strip()]
    meta_rows = {r["sample"]: r for r in read_tsv(samples_tsv)}
    versions = read_text(p("software_versions.txt"))

    samples = sorted(set(nuclear["PSC"]) | set(meta_rows))
    group = {s: meta_rows.get(s, {}).get("group", "") for s in samples}
    label_sex = {s: norm_sex(meta_rows.get(s, {}).get("sex", "")) for s in samples}

    # per-sample metrics
    m = {s: {} for s in samples}
    for s in samples:
        psc = nuclear["PSC"].get(s, {})
        m[s]["depth"] = psc.get("average depth")
        m[s]["het"] = het_rate(psc)
        m[s]["singletons"] = psc.get("nSingletons")
        ma = miss_after.get(s)
        m[s]["missing"] = float(ma["pct_missing"]) if ma and ma.get("pct_missing") not in (None, "") else None
        if sexstats is not None:
            xs = het_rate(sexstats["PSC"].get(s, {}))
            m[s]["x_ratio"] = xs / m[s]["het"] if xs is not None and m[s]["het"] else None
            m[s]["sex_inferred"] = infer_sex(m[s]["x_ratio"])
        if mitostats is not None:
            mp = mitostats["PSC"].get(s, {})
            nonref = mp.get("nHets", 0) + mp.get("nNonRefHom", 0)
            m[s]["mito_het"] = mp.get("nHets", 0) / nonref if nonref else None
            md = mp.get("average depth")
            m[s]["mito_depth_ratio"] = md / m[s]["depth"] if md and m[s]["depth"] else None

    flags = {s: [] for s in samples}
    for key, side, what in (("missing", "high", "high missingness"), ("depth", "both", "depth outlier"),
                            ("het", "both", "heterozygosity outlier"), ("singletons", "high", "many singletons")):
        for s, z in robust_outliers({s: m[s].get(key) for s in samples}, side).items():
            flags[s].append(f"{what} (z = {z:+.1f})")
    if sexstats is not None:
        for s in samples:
            inf, lab = m[s].get("sex_inferred"), label_sex[s]
            if lab and inf in ("M", "F") and inf != lab:
                flags[s].append(f"sex mismatch (labelled {lab}, chrX says {inf})")
            elif inf == "ambiguous":
                flags[s].append("ambiguous chrX heterozygosity")
    if mitostats is not None:
        for s in samples:
            h = m[s].get("mito_het")
            if h is not None and h > MITO_HET_MAX:
                flags[s].append(f"mito heteroplasmy / contamination ({h:.0%} of non-ref mito calls heterozygous)")
        for s, z in robust_outliers({s: m[s].get("mito_depth_ratio") for s in samples}).items():
            flags[s].append(f"mito depth ratio outlier (z = {z:+.1f})")

    kin_rows = [(r["ID1"], r["ID2"], float(r["kinship"]), float(r["IBS0"])) for r in kin]
    dups = [k for k in kin_rows if k[2] > KING_DUPLICATE]
    for a, b, k, _ in dups:
        flags.setdefault(a, []).append(f"duplicate of {b} (kinship {k:.2f})")
        flags.setdefault(b, []).append(f"duplicate of {a} (kinship {k:.2f})")
    flagged = {s: f for s, f in flags.items() if f}

    # ── page ──
    run = meta.get("run_name", "")
    parts = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
             f"<meta name='viewport' content='width=device-width, initial-scale=1'>"
             f"<title>Cohort QC: {esc(run)}</title><style>{CSS}</style></head><body>",
             f"<h1>Cohort QC report <span class='sub'>{esc(run)}</span></h1>",
             f"<p class='lede'>{len(samples):,} samples · "
             f"{fmt(after['SN'].get('number of records'))} PASS sites after filtering · "
             f"sex contig: {esc(meta.get('sex_chrom') or 'not set')} · "
             f"mito contig: {esc(meta.get('mito_chrom') or 'not set')}</p>"]

    # 1. before the results
    n_exp, n_vcf = meta.get("n_samples_expected"), len(nuclear["PSC"])
    pieces_exp, pieces_got = meta.get("n_pieces_expected"), meta.get("n_pieces_gathered")
    checks = [
        status(n_vcf == n_exp, f"<b>Samples:</b> {n_vcf:,} in the cohort VCF, {fmt(n_exp)} in the samplesheet"),
        status(pieces_got == pieces_exp, f"<b>Genome:</b> {fmt(pieces_got)} of {fmt(pieces_exp)} joint-genotyping pieces gathered"),
        status(not dups, f"<b>Duplicates:</b> {len(dups)} sample pair(s) with KING kinship &gt; {KING_DUPLICATE}"),
    ]
    if sexstats is not None:
        mism = [s for s in samples if any(f.startswith("sex mismatch") for f in flags[s])]
        labelled = sum(1 for s in samples if label_sex[s])
        checks.append(status(not mism if labelled else None,
                             f"<b>Sex:</b> {len(mism)} mismatch(es) between labels and chrX"
                             + ("" if labelled else " (no <code>sex</code> column in the samplesheet: inferred only)")))
    if mitostats is not None:
        het_flag = [s for s in samples if any(f.startswith("mito heteroplasmy") for f in flags[s])]
        checks.append(status(not het_flag, f"<b>Mitochondria:</b> {len(het_flag)} sample(s) with heterozygous mito calls above {MITO_HET_MAX:.0%}"))
    checks.append(status(not flagged if samples else None, f"<b>Per-sample flags:</b> {len(flagged)} sample(s) flagged"))
    parts.append("<h2>Before the results</h2><ul class='checks'>" + "".join(checks) + "</ul>")
    if flagged:
        parts.append(table(["sample", "group", "why"],
                           [[esc(s), esc(group.get(s, "")), "<br>".join(esc(f) for f in fl)]
                            for s, fl in sorted(flagged.items())]))
    for n in notes:
        parts.append(f"<p class='note'>{esc(n)}</p>")

    # 2. PCA
    parts.append("<h2>Population structure (nuclear PCA)</h2>")
    parts.append(f"<p class='note'>PASS biallelic SNPs, minor allele frequency ≥ {esc(meta.get('qc_maf'))}, "
                 f"one per kb, LD-pruned at r² {esc(meta.get('qc_ld_threshold'))}"
                 + (", sex and mito contigs excluded" if meta.get("sex_chrom") or meta.get("mito_chrom") else "")
                 + f": {fmt(snp_counts.get('nuclear_ld_pruned'))} of {fmt(snp_counts.get('nuclear_input'))} SNPs used. "
                 "Colour = samplesheet <code>group</code> (or sex label); points far from the rest are labelled.</p>")
    parts.append(pca_panels(pca_n, eig_n, group, label_sex))

    # 3. relatedness
    parts.append("<h2>Relatedness (KING-robust kinship)</h2>")
    deg = {}
    for _, _, k, _ in kin_rows:
        deg[king_degree(k)] = deg.get(king_degree(k), 0) + 1
    parts.append("<div class='row'><div>" + svg_hist([k for _, _, k, _ in kin_rows], "kinship, all sample pairs", bins=40,
                 marks=((KING_DUPLICATE, "dup"), (KING_FIRST, "1st"), (KING_SECOND, "2nd")), log_y=True)
                 + "</div><div>" + table(["relationship", "pairs"],
                 [[esc(d), fmt(deg.get(d, 0))] for d in ("duplicate", "first degree", "second degree", "third degree", "unrelated")])
                 + "<p class='note'>Full sibs and parent–offspring are first degree (≈ 0.25). In a cross, check the "
                 "pattern matches the pedigree.</p></div></div>")
    if dups:
        parts.append(table(["sample 1", "sample 2", "kinship", "IBS0"],
                           [[esc(a), esc(b), f"{k:.3f}", f"{i:.4f}"] for a, b, k, i in sorted(dups, key=lambda x: -x[2])]))

    # 4. sex chromosome
    if sexstats is not None:
        parts.append(f"<h2>Sex chromosome ({esc(meta['sex_chrom'])})</h2>")
        pts = [(m[s]["het"], m[s].get("x_ratio"), s, label_sex[s] or "unlabelled") for s in samples]
        cols = {"M": PALETTE[0], "F": PALETTE[1], "unlabelled": GREY}
        parts.append("<div class='row'><div>" + svg_scatter(pts, "autosomal heterozygosity", "chrX / autosome heterozygosity",
                     cols, {s for s in samples if any("sex" in f or "chrX" in f for f in flags[s])}) + "</div><div>"
                     + f"<p class='note'>One X (males, XO) gives almost no heterozygous chrX calls: ratio ≤ {SEX_MALE_MAX} "
                     f"is called M, ≥ {SEX_FEMALE_MIN} F, in between ambiguous. Colour = samplesheet <code>sex</code>.</p>"
                     + table(["inferred", "samples"], [[esc(k), fmt(sum(1 for s in samples if m[s].get('sex_inferred') == k))]
                                                       for k in ("M", "F", "ambiguous", "?")]) + "</div></div>")

    # 5. mitochondria
    if mitostats is not None:
        parts.append(f"<h2>Mitochondria ({esc(meta['mito_chrom'])})</h2>")
        parts.append(f"<p class='note'>Maternal lineages: all {fmt(snp_counts.get('mito_input'))} biallelic PASS mito SNPs, no LD pruning. "
                     "Heteroplasmy = heterozygous fraction of non-reference mito calls; mito depth ratio = mito ÷ nuclear depth.</p>")
        parts.append(pca_panels(pca_m, eig_m, group, label_sex, pcs=((1, 2),)) if pca_m else "<p class='note'>No mito PCA (see notes above).</p>")
        parts.append("<div class='row'><div>" + svg_hist([m[s].get("mito_het") for s in samples], "heterozygous fraction of non-ref mito calls",
                     marks=((MITO_HET_MAX, "flag"),)) + "</div><div>"
                     + svg_hist([m[s].get("mito_depth_ratio") for s in samples], "mito / nuclear depth") + "</div></div>")

    # 6. per-sample table
    parts.append("<h2>Per sample</h2>")
    hdr = ["sample", "group", "depth", "% missing (after)", "het rate", "singletons"]
    if sexstats is not None:
        hdr += ["sex label", "chrX ratio", "sex (chrX)"]
    if mitostats is not None:
        hdr += ["mito het", "mito depth ratio"]
    hdr.append("flags")
    rows = []
    for s in samples:
        r = [esc(s), esc(group.get(s, "")), fmt(m[s]["depth"]), fmt(m[s]["missing"]), fmt(m[s]["het"]), fmt(m[s]["singletons"])]
        if sexstats is not None:
            r += [esc(label_sex[s]), fmt(m[s].get("x_ratio")), esc(m[s].get("sex_inferred", ""))]
        if mitostats is not None:
            r += [fmt(m[s].get("mito_het")), fmt(m[s].get("mito_depth_ratio"))]
        r.append(f"<span class='bad'>{len(flags[s])}</span>" if flags[s] else "")
        rows.append(r)
    parts.append("<div class='scroll'>" + table(hdr, rows, "persample") + "</div>")

    # 7. variants
    parts.append("<h2>Variants</h2>")
    keys = ("number of records", "number of SNPs", "number of indels", "number of multiallelic sites", "number of MNPs")
    vrows = [[esc(k.replace("number of ", "")), fmt(before["SN"].get(k)), fmt(after["SN"].get(k))] for k in keys]
    vrows.append(["Ts/Tv", fmt(before["TSTV"]), fmt(after["TSTV"])])
    parts.append("<div class='row'><div>" + table(["", "before filtering (all sites)", "after (PASS)"], vrows) + "</div><div>"
                 + (table(["type", "filter", "expression", "sites", "% of type", "only this filter"],
                          [[esc(r.get("variant_type")), esc(r.get("filter")), esc(r.get("expression")), fmt(_num(r.get("n_sites"))),
                            fmt(_num(r.get("pct_sites"))), fmt(_num(r.get("n_only")))] for r in filt]) if filt else "")
                 + "</div></div>")
    parts.append("<h3>PASS sites per Mb along the genome</h3><p class='note'>Gaps are regions with no PASS sites "
                 "(e.g. assembly gaps, or a missing window); spikes often mark collapsed repeats.</p>")
    parts.append(svg_density(density))

    # 8. methods and provenance
    parts.append("<h2>Methods and definitions</h2><ul class='methods'>"
                 f"<li><b>Outliers</b>: robust z = (x − median) ÷ (1.4826 × MAD); flagged beyond ±{ROBUST_Z} (needs ≥ {MIN_SAMPLES_FOR_OUTLIERS} samples). "
                 "Missingness and singletons are flagged on the high side only.</li>"
                 "<li><b>Het rate</b>: heterozygous fraction of called genotypes at PASS variant sites on nuclear contigs (not per bp).</li>"
                 f"<li><b>Kinship</b>: SNPRelate KING-robust on the LD-pruned SNPs; duplicate &gt; {KING_DUPLICATE}, first degree &gt; {KING_FIRST}, "
                 f"second &gt; {KING_SECOND}, third &gt; {KING_THIRD}.</li>"
                 "<li><b>Per-sample depth, singletons</b>: bcftools stats (PSC) over PASS sites on nuclear contigs.</li></ul>")
    if versions:
        parts.append(f"<details><summary>Software versions and command</summary><pre>{esc(versions)}</pre></details>")
    parts.append("</body></html>")
    return "\n".join(parts)


def pca_panels(pca: list, eig: list, group: dict, label_sex: dict, pcs=((1, 2), (3, 4))) -> str:
    if not pca:
        return '<p class="note">No PCA output.</p>'
    ve = {r["pc"]: float(r["variance_explained"]) for r in eig}
    use_group = any(group.values())
    key = (lambda s: group.get(s, "")) if use_group else (lambda s: label_sex.get(s, ""))
    cols = colour_map(key(r["sample"]) for r in pca)
    out = ["<div class='row'>"]
    for a, b in pcs:
        ca, cb = f"PC{a}", f"PC{b}"
        if ca not in pca[0] or cb not in pca[0]:
            continue
        pts = [(float(r[ca]), float(r[cb]), r["sample"], key(r["sample"])) for r in pca]
        far = set(robust_outliers({p[2]: p[0] for p in pts})) | set(robust_outliers({p[2]: p[1] for p in pts}))
        out.append("<div>" + svg_scatter(pts, f"{ca} ({ve.get(ca, 0):.1%})", f"{cb} ({ve.get(cb, 0):.1%})", cols, far) + "</div>")
    out.append("</div>")
    if cols:
        out.append("<p class='legend'>" + " ".join(f"<span style='color:{c}'>●</span> {esc(k)}" for k, c in cols.items()) + "</p>")
    return "".join(out)


CSS = """
:root{--fg:#1d2329;--muted:#5f6b76;--line:#d9dee3;--bg:#fff;--ok:#1a7f37;--bad:#c62828;--warn:#b26a00}
@media (prefers-color-scheme:dark){:root{--fg:#e6e9ec;--muted:#9aa5af;--line:#3a424a;--bg:#14181c}}
body{font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--fg);background:var(--bg);max-width:1000px;margin:24px auto;padding:0 16px}
h1{font-size:22px;margin:0 0 4px}h1 .sub{color:var(--muted);font-weight:400}
h2{font-size:17px;margin:28px 0 8px;padding-top:8px;border-top:1px solid var(--line)}h3{font-size:14px}
.lede,.note,.legend{color:var(--muted)}.note{font-size:13px}
table{border-collapse:collapse;font-size:13px;margin:8px 0}th,td{padding:3px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{font-weight:600}.scroll{overflow-x:auto}
ul.checks{list-style:none;padding:0}ul.checks li{margin:4px 0}.mark{display:inline-block;width:1.2em;font-weight:700}
li.ok .mark{color:var(--ok)}li.bad .mark,span.bad{color:var(--bad)}li.warn .mark{color:var(--warn)}
.row{display:flex;flex-wrap:wrap;gap:16px}.row>div{flex:1 1 300px;min-width:0}
svg text{fill:var(--fg)}.tick{font-size:10px;fill:var(--muted)}.axis{font-size:11px}.lab{font-size:10px}
.frame{fill:none;stroke:var(--line)}.mark{stroke:var(--bad);stroke-dasharray:3 3}
pre{white-space:pre-wrap;font-size:12px}code{font-size:12px}
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs", required=True, help="directory with the QC input files")
    ap.add_argument("--samples", required=True, help="TSV: sample, sex, group")
    ap.add_argument("--meta", required=True, help="JSON: run metadata and parameters")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    with open(a.meta) as f:
        meta = json.load(f)
    html_out = build(a.inputs, a.samples, meta)
    with open(a.out, "w") as f:
        f.write(html_out)
    print(f"QC report written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
