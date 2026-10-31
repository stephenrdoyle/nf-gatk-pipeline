#!/usr/bin/env python3
"""
Compare site-filtering schemes on an existing cohort VCF, without rerunning
anything: what each scheme keeps, and what each individual filter removes.

Input is the per-site table written by run_comparison.sh (bcftools query on
cohort.filtered.vcf.gz after +fill-tags recomputes AC/AN/F_MISSING from the
genotype-filtered calls), columns:

    CHROM POS REF ALT QUAL FILTER QD FS SOR MQ MQRankSum ReadPosRankSum
    DP ExcessHet F_MISSING AC AN

Two passes over the table: the first finds the median site depth (per
variant type) for the depth ceilings; the second evaluates every filter at
every site and tallies the combinations. Every scheme is then computed from
those tallies, so adding a scheme costs nothing.

How to read the output:
  - Ts/Tv of PASS SNPs: higher is cleaner. Random errors give about 0.5;
    real variation in nematodes is typically around 2.
  - "Removed by this filter alone": sites failing this filter and no other
    in the scheme. Their Ts/Tv says what the filter is removing: close to the
    PASS Ts/Tv means mostly real sites; near 0.5-1 means mostly errors.

Missing annotations pass, as in GATK VariantFiltration (an expression on a
missing value does not flag the site).

Standard library only.

Usage: compare_filters.py sites.tsv.gz > filter_comparison.md
"""

import subprocess
import sys
from collections import Counter

# GATK hard-filter thresholds (the pipeline's 'hard' mode); missing = pass
HARD = {
    True:  [("QD", "<", 2.0), ("MQ", "<", 40.0), ("FS", ">", 60.0), ("SOR", ">", 3.0),
            ("MQRankSum", "<", -12.5), ("ReadPosRankSum", "<", -8.0)],
    False: [("QD", "<", 2.0), ("FS", ">", 200.0), ("SOR", ">", 10.0), ("ReadPosRankSum", "<", -20.0)],
}
QUAL_MIN = 30.0
EXCESS_HET_MAX = 54.69          # GATK germline short-variant guidance
DP_MULTIPLES = (1.5, 2.0, 3.0)  # depth ceiling = multiple x median site DP (per type)
FMISS_MAX = (0.1, 0.2, 0.5)     # fraction of samples with no call after the genotype filter

COLS = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "QD", "FS", "SOR", "MQ",
        "MQRankSum", "ReadPosRankSum", "DP", "ExcessHet", "F_MISSING", "AC", "AN"]
IDX = {c: i for i, c in enumerate(COLS)}
TRANSITIONS = {frozenset("AG"), frozenset("CT")}


def rows(path):
    """Fields of each line; reads .gz through gzip -dc (much faster than the gzip module)."""
    if path.endswith(".gz"):
        proc = subprocess.Popen(["gzip", "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        stream = proc.stdout
    else:
        proc, stream = None, open(path)
    for line in stream:
        yield line.rstrip("\n").split("\t")
    if proc and proc.wait() != 0:
        sys.exit(f"ERROR: gzip -dc {path} failed")


def num(v):
    if v in (".", ""):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def classify(ref, alt):
    """(is_snp, tstv) with the pipeline's split rule: SNP if REF and every
    ALT (ignoring '*') are single bases; anything else goes with the indels.
    tstv is 'ts'/'tv' for biallelic SNPs, else None."""
    alts = [a for a in alt.split(",") if a not in ("*", ".")]
    if not alts:
        return None, None
    is_snp = len(ref) == 1 and all(len(a) == 1 for a in alts)
    if is_snp and len(alts) == 1:
        return True, "ts" if frozenset((ref.upper(), alts[0].upper())) in TRANSITIONS else "tv"
    return is_snp, None


def median_from_counts(counts: Counter):
    total = sum(counts.values())
    if not total:
        return None
    half, seen = total / 2, 0
    for v in sorted(counts):
        seen += counts[v]
        if seen >= half:
            return v


def main(path):
    # ── pass 1: site depth distribution per type ──
    dp_counts = {True: Counter(), False: Counter()}
    for f in rows(path):
        is_snp, _ = classify(f[IDX["REF"]], f[IDX["ALT"]])
        dp = num(f[IDX["DP"]])
        if is_snp is not None and dp is not None:
            dp_counts[is_snp][int(dp)] += 1
    med = {t: median_from_counts(c) for t, c in dp_counts.items()}

    # ── pass 2: every atomic filter at every site, tallied by combination ──
    names, bit = [], {}

    def b(name):
        if name not in bit:
            bit[name] = 1 << len(names)
            names.append(name)
        return bit[name]

    tally = Counter()   # (is_snp, mask, tstv, has_alt) -> sites
    skipped = 0
    for f in rows(path):
        is_snp, tstv = classify(f[IDX["REF"]], f[IDX["ALT"]])
        if is_snp is None:
            skipped += 1
            continue
        mask = 0
        # current pipeline result: each optimised filter by annotation
        for flt in f[IDX["FILTER"]].split(";"):
            if flt not in ("PASS", "."):
                parts = flt.split("_")
                mask |= b("opt:" + (parts[1] if len(parts) > 2 and parts[0] in ("SNP", "INDEL") else flt))
        for ann, op, thr in HARD[is_snp]:
            v = num(f[IDX[ann]])
            if v is not None and (v < thr if op == "<" else v > thr):
                mask |= b("hard:" + ann)
        q = num(f[IDX["QUAL"]])
        if q is not None and q < QUAL_MIN:
            mask |= b("QUAL<30")
        dp = num(f[IDX["DP"]])
        if dp is not None and med[is_snp]:
            for k in DP_MULTIPLES:
                if dp > k * med[is_snp]:
                    mask |= b(f"DP>{k:g}x")
        fm = num(f[IDX["F_MISSING"]])
        if fm is not None:
            for k in FMISS_MAX:
                if fm > k:
                    mask |= b(f"missing>{k:g}")
        eh = num(f[IDX["ExcessHet"]])
        if eh is not None and eh > EXCESS_HET_MAX:
            mask |= b("ExcessHet")
        ac = f[IDX["AC"]]
        has_alt = ac not in (".", "") and sum(int(x) for x in ac.split(",") if x not in (".", "")) > 0
        tally[(is_snp, mask, tstv, has_alt)] += 1

    report(tally, names, bit, med, skipped)


def report(tally, names, bit, med, skipped):
    opt = [n for n in names if n.startswith("opt:")]
    hard = [n for n in names if n.startswith("hard:")]
    q = [n for n in names if n == "QUAL<30"]
    base = hard + q
    have = set(names)
    sel = lambda *ns: [n for n in ns if n in have]

    schemes = [
        ("current (optimised, 5%)", opt),
        ("GATK hard", hard),
        ("hard + QUAL<30", base),
    ]
    for k in DP_MULTIPLES:
        schemes.append((f"hard + QUAL + DP>{k:g}x median", base + sel(f"DP>{k:g}x")))
    for k in FMISS_MAX:
        schemes.append((f"hard + QUAL + DP>2x + missing>{k:g}", base + sel("DP>2x", f"missing>{k:g}")))
    proposed = base + sel("DP>2x", "missing>0.2")
    schemes.append(("proposed + ExcessHet", proposed + sel("ExcessHet")))

    mask_of = lambda ns: sum(bit[n] for n in ns)
    opt_mask = mask_of(opt)

    def stats(is_snp, fmask):
        n = n_pass = n_pass_alt = ts = tv = lost = gained = 0
        for (s, m, tstv, has_alt), c in tally.items():
            if s != is_snp:
                continue
            n += c
            ok = not (m & fmask)
            cur_ok = not (m & opt_mask)
            if ok:
                n_pass += c
                n_pass_alt += c if has_alt else 0
                if tstv == "ts":
                    ts += c
                elif tstv == "tv":
                    tv += c
            if cur_ok and not ok:
                lost += c
            if ok and not cur_ok:
                gained += c
        return n, n_pass, n_pass_alt, (ts / tv if tv else None), lost, gained

    def alone(is_snp, ns):
        """Per filter in the scheme: flagged, flagged by it alone, Ts/Tv of those alone."""
        out = []
        fmask = mask_of(ns)
        for name in ns:
            fb = bit[name]
            flagged = only = ts = tv = 0
            for (s, m, tstv, _), c in tally.items():
                if s != is_snp or not (m & fb):
                    continue
                flagged += c
                if (m & fmask) == fb:
                    only += c
                    ts += c if tstv == "ts" else 0
                    tv += c if tstv == "tv" else 0
            if flagged:
                out.append((name, flagged, only, ts / tv if tv else None))
        return out

    f = lambda x: "–" if x is None else f"{x:,}" if isinstance(x, int) else f"{x:.2f}"
    p = print
    p("# Site filter comparison\n")
    p(f"Median site depth (INFO/DP): SNPs {f(med[True])}, indels {f(med[False])}. "
      f"Records without a usable ALT skipped: {skipped:,}.\n")
    p("Ts/Tv is for biallelic PASS SNPs. 'with a called ALT' = at least one non-reference "
      "genotype left after the genotype filter. lost / gained are relative to the current PASS set.\n")
    for is_snp, label in ((True, "SNPs"), (False, "Indels (incl. mixed)")):
        p(f"## {label}\n")
        p("| scheme | PASS | % | PASS with a called ALT | Ts/Tv | lost vs current | gained vs current |")
        p("|---|---|---|---|---|---|---|")
        for name, ns in schemes:
            n, np_, npa, tt, lost, gained = stats(is_snp, mask_of(ns))
            p(f"| {name} | {f(np_)} | {100 * np_ / n:.1f} | {f(npa)} | {f(tt) if is_snp else '–'} "
              f"| {f(lost)} | {f(gained)} |")
        n_all = stats(is_snp, 0)
        p(f"\nAll sites: {f(n_all[0])}; with a called ALT after the genotype filter: {f(n_all[2])}.\n")
        for title, ns in (("current (optimised, 5%)", opt), ("proposed + ExcessHet", proposed + sel("ExcessHet"))):
            p(f"### What each filter removes: {title}\n")
            p("| filter | flagged | removed by this filter alone | Ts/Tv of those |")
            p("|---|---|---|---|")
            for name, flagged, only, tt in alone(is_snp, ns):
                p(f"| {name} | {f(flagged)} | {f(only)} | {f(tt) if is_snp else '–'} |")
            p("")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
