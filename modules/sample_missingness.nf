process SAMPLE_MISSINGNESS {
    tag "${label}"
    label 'process_low'

    publishDir "${params.outdir}/missingness", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    val  label

    output:
    path "${label}_missingness.tsv", emit: missingness

    script:
    """
    # ── 1. per-sample missingness from bcftools stats (single fast pass) ──────
    bcftools stats -s - ${vcf} > bcftools_stats.txt

    # ── 2. LowDP counts per sample per type — single streaming pass ───────────
    # Decompress once, pipe through two awk processes simultaneously using
    # process substitution to avoid reading the file twice.
    # awk extracts: SAMPLE, FORMAT tags, FILTER, TYPE
    # We only need sites where at least one FORMAT field contains LowDP
    bcftools view -H ${vcf} | \\
        awk '
        BEGIN { OFS="\\t" }
        {
            # Field 4 = REF, 5 = ALT — determine type
            ref = \$4; alt = \$5
            type = (length(ref) == length(alt) && length(ref) == 1) ? "SNP" : "INDEL"

            # FORMAT field is column 9
            fmt = \$9
            n = split(fmt, fmt_keys, ":")

            # Find FT index in FORMAT
            ft_idx = 0
            for (i = 1; i <= n; i++) {
                if (fmt_keys[i] == "FT") { ft_idx = i; break }
            }
            if (ft_idx == 0) next  # no FT field, skip

            # Iterate sample columns (10+)
            for (s = 10; s <= NF; s++) {
                split(\$s, vals, ":")
                ft = (ft_idx <= length(vals)) ? vals[ft_idx] : "."
                if (ft ~ /LowDP/) {
                    print s-9, type  # sample index (1-based), type
                }
            }
        }
        ' > lowdp_hits.txt

    # ── 3. get sample names in order ──────────────────────────────────────────
    bcftools query -l ${vcf} > sample_names.txt

    # ── 4. combine into final TSV ──────────────────────────────────────────────
    python3 << 'PYEOF'
from collections import defaultdict

label       = "${label}"
out_file    = f"{label}_missingness.tsv"
total_sites = 0
samples     = {}

# Read sample names (ordered)
with open('sample_names.txt') as f:
    sample_names = [l.strip() for l in f if l.strip()]

# Parse bcftools stats
with open('bcftools_stats.txt') as f:
    for line in f:
        if line.startswith('SN') and 'number of records' in line:
            total_sites = int(line.strip().split()[-1])
        if line.startswith('PSC') and not line.startswith('# PSC'):
            parts = line.strip().split('\\t')
            if len(parts) >= 9:
                sample      = parts[2]
                missing     = int(parts[10]) if len(parts) > 10 else 0
                hom_rr      = int(parts[3])
                hom_aa      = int(parts[4])
                het         = int(parts[5])
                called      = hom_rr + hom_aa + het
                pct_missing = round(100 * missing / total_sites, 4) if total_sites else 0
                pct_called  = round(100 * called  / total_sites, 4) if total_sites else 0
                samples[sample] = {
                    'total_sites':   total_sites,
                    'n_called':      called,
                    'n_missing':     missing,
                    'pct_called':    pct_called,
                    'pct_missing':   pct_missing,
                    'n_lowdp_snp':   0,
                    'n_lowdp_indel': 0,
                }

# Count LowDP hits from awk output (sample_index is 1-based)
with open('lowdp_hits.txt') as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) < 2:
            continue
        try:
            sidx  = int(parts[0]) - 1
            vtype = parts[1]
        except ValueError:
            continue
        if sidx < 0 or sidx >= len(sample_names):
            continue
        sample = sample_names[sidx]
        if sample not in samples:
            continue
        if vtype == 'SNP':
            samples[sample]['n_lowdp_snp']   += 1
        else:
            samples[sample]['n_lowdp_indel'] += 1

with open(out_file, 'w') as f:
    header = ('sample\\ttotal_sites\\tn_called\\tn_missing\\tpct_called\\t'
              'pct_missing\\tn_lowdp_snp\\tn_lowdp_indel\\tlabel\\n')
    f.write(header)
    for sample, d in samples.items():
        f.write('\\t'.join(str(x) for x in [
            sample,
            d['total_sites'], d['n_called'],    d['n_missing'],
            d['pct_called'],  d['pct_missing'],
            d['n_lowdp_snp'], d['n_lowdp_indel'],
            label
        ]) + '\\n')

print(f"Missingness report written: {out_file} ({len(samples)} samples)")
PYEOF
    """
}
