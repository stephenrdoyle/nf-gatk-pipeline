process SAMPLE_MISSINGNESS {
    /*
     * Per-sample genotype missingness for one stage of the cohort VCF.
     *
     * Run twice (see main.nf):
     *   before_filter  cohort.genotyped.vcf.gz, all sites: genotypes as
     *                  GATK called them
     *   after_filter   cohort.filtered.vcf.gz, PASS sites only: the data you
     *                  keep, after per-genotype DP/GQ no-calls
     * Low-DP counts are genotypes with FORMAT/DP below --genotype_dp_threshold
     * (the genotype filter's DP rule), counted the same way at both stages.
     * pass_only restricts both the counts and the low-DP tally to PASS sites,
     * so failed sites (which are discarded anyway) don't inflate missingness.
     *
     * Called/missing counts come from BCFTOOLS_STATS for the same stage
     * (same VCF, same pass_only) rather than a second bcftools stats pass.
     */
    tag "${label}"
    label 'process_low'

    publishDir "${params.outdir}/missingness", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path stats       // BCFTOOLS_STATS output for this stage (same pass_only)
    val  label       // "before_filter" or "after_filter"
    val  pass_only   // true: only sites whose FILTER is PASS (or '.')

    output:
    path "${label}_missingness.tsv", emit: missingness

    script:
    def filter_arg = pass_only ? "--apply-filters 'PASS,.'" : ''
    """
    # Nextflow runs bash -ue only: without pipefail a bcftools failure in the
    # pipe below would leave empty low-DP counts and a task that "succeeded"
    set -o pipefail

    # ── 1. per-sample called / missing genotypes: from ${stats} (input) ──────
    # ── 2. Low-depth genotypes (FORMAT/DP < --genotype_dp_threshold) per
    #       sample, split SNP / indel; nothing is counted if no threshold ──
    # A site is a SNP when REF and every ALT allele are single bases
    # (multiallelic SNPs included); anything else counts as an indel.
    bcftools view --threads ${task.cpus} ${filter_arg} --no-header ${vcf} | \\
        awk -v thr='${params.genotype_dp_threshold ?: 0}' '
        BEGIN { OFS="\\t"; skip = (thr + 0 <= 0) }   # read all input (no early exit: SIGPIPE under pipefail)
        {
            if (skip) next
            n_alt = split(\$5, alts, ",")
            type = (length(\$4) == 1) ? "SNP" : "INDEL"
            for (a = 1; a <= n_alt; a++)
                if (length(alts[a]) != 1 || alts[a] == "*") type = "INDEL"

            n = split(\$9, fmt_keys, ":")
            dp_idx = 0
            for (i = 1; i <= n; i++) if (fmt_keys[i] == "DP") { dp_idx = i; break }
            if (dp_idx == 0) next   # no FORMAT/DP on this record

            for (s = 10; s <= NF; s++) {
                n_vals = split(\$s, vals, ":")
                if (dp_idx <= n_vals && vals[dp_idx] != "." && vals[dp_idx] + 0 < thr + 0) print s - 9, type
            }
        }
        ' > lowdp_hits.txt

    # ── 3. sample names in VCF order ──────────────────────────────────────────
    bcftools query --list-samples ${vcf} > sample_names.txt

    # ── 4. combine into the per-sample TSV ────────────────────────────────────
    python3 << 'PYEOF'
label       = "${label}"
out_file    = f"{label}_missingness.tsv"
total_sites = None
samples     = {}

with open('sample_names.txt') as f:
    sample_names = [l.strip() for l in f if l.strip()]

# PSC columns are located by name from the '# PSC' header line rather than
# by position, because the layout varies between bcftools versions
# (in 1.19 nMissing is column [14], not [11] = nSingletons)
psc_cols = None
with open('${stats}') as f:
    for line in f:
        parts = line.rstrip('\\n').split('\\t')
        if line.startswith('SN') and 'number of records' in line:
            total_sites = int(parts[-1])
        elif line.startswith('# PSC'):
            psc_cols = {h.split(']', 1)[-1]: i for i, h in enumerate(parts)}   # '[4]nRefHom' -> 'nRefHom'
        elif line.startswith('PSC'):
            if psc_cols is None or 'nMissing' not in psc_cols:
                raise SystemExit('ERROR: bcftools stats PSC section has no nMissing column')
            called  = sum(int(parts[psc_cols[k]]) for k in ('nRefHom', 'nNonRefHom', 'nHets'))
            missing = int(parts[psc_cols['nMissing']])
            samples[parts[psc_cols['sample']]] = {
                'n_called': called, 'n_missing': missing,
                'n_lowdp_snp': 0, 'n_lowdp_indel': 0,
            }

# Fail rather than write an empty or wrong table
if total_sites is None:
    raise SystemExit('ERROR: no "number of records" line in bcftools stats output')
if set(samples) != set(sample_names):
    raise SystemExit(f'ERROR: bcftools stats has {len(samples)} samples, the VCF has {len(sample_names)}')

with open('lowdp_hits.txt') as f:
    for line in f:
        idx, vtype = line.rstrip('\\n').split('\\t')
        key = 'n_lowdp_snp' if vtype == 'SNP' else 'n_lowdp_indel'
        samples[sample_names[int(idx) - 1]][key] += 1

pct = lambda n: round(100 * n / total_sites, 4) if total_sites else 0
with open(out_file, 'w') as f:
    f.write('\\t'.join(['sample', 'total_sites', 'n_called', 'n_missing', 'pct_called',
                       'pct_missing', 'n_lowdp_snp', 'n_lowdp_indel', 'label']) + '\\n')
    for sample in sample_names:
        d = samples[sample]
        f.write('\\t'.join(str(x) for x in [
            sample, total_sites, d['n_called'], d['n_missing'],
            pct(d['n_called']), pct(d['n_missing']),
            d['n_lowdp_snp'], d['n_lowdp_indel'], label,
        ]) + '\\n')

print(f'Missingness written: {out_file} ({len(sample_names)} samples, {total_sites} sites)')
PYEOF
    """
}
