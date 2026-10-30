/*
 * --heterozygosity auto: estimate the cohort's heterozygosity from the BAMs,
 * for GenotypeGVCFs' --heterozygosity / --indel-heterozygosity priors
 * (GATK's defaults, 0.001 / 1.25e-4, assume human-like diversity).
 *
 * ESTIMATE_HETEROZYGOSITY  per sample, ANGSD in a container (not a farm
 *                          module): site frequency spectrum from genotype
 *                          likelihoods over a fixed sample of the genome, so
 *                          no variant calls and no caller prior are involved.
 *                          Heterozygosity = het fraction of the 1-sample SFS.
 * HETEROZYGOSITY_SUMMARY   median across samples (samples with under half the
 *                          median number of covered sites are excluded), rounded
 *                          to 3 significant figures so the value - and so
 *                          GenotypeGVCFs' cache - is stable between runs.
 *
 * For a randomly mating population an individual's expected heterozygosity
 * equals the nucleotide diversity GATK's prior stands for. Per-sample
 * estimates also keep memory flat for any cohort size. Requires paired-end
 * reads (ANGSD's default proper-pair filter).
 */

process ESTIMATE_HETEROZYGOSITY {
    tag "${sample_id}"
    label 'process_low'
    // container (params.angsd_container) and a beforeScript that loads
    // Singularity are set in nextflow.config

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    val  sample_bp     // genome sampled, bp (--heterozygosity_sample_bp)

    output:
    path "${sample_id}.sfs", emit: sfs

    script:
    """
    # Deterministic, evenly spaced genome sample of ~sample_bp: K windows at
    # equal steps along the concatenated genome, each at most 1 Mb, at most
    # the step (so they never overlap) and cut at contig ends. On fragmented
    # assemblies cut windows fall short, so K grows until the sample reaches
    # 95% of sample_bp. Every contig whole if the genome is smaller than
    # sample_bp. The same .fai always gives the same regions.
    awk -F'\\t' -v target=${sample_bp} -v W=1000000 '
        function sample(K,   k, pos, c, offset, s, e, w, step, total) {
            step = G / K; w = int(W < step ? W : step); c = 1; offset = 0; total = 0; n = 0
            for (k = 0; k < K; k++) {
                pos = int(k * step)
                while (pos >= offset + len[c]) { offset += len[c]; c++ }
                s = pos - offset + 1; e = s + w - 1; if (e > len[c]) e = len[c]
                out[++n] = name[c] ":" s "-" e; total += e - s + 1
            }
            return total
        }
        { name[NR] = \$1; len[NR] = \$2; G += \$2 }
        END {
            if (G <= target) { for (i = 1; i <= NR; i++) print name[i] ":1-" len[i]; exit }
            K = int((target + W - 1) / W)
            for (iter = 0; iter < 8; iter++) {
                total = sample(K)
                if (total >= 0.95 * target) break
                K = int(K * target / total * 1.05) + 1
            }
            for (i = 1; i <= n; i++) print out[i]
        }' ${reference_fai} > regions.txt

    angsd -i ${bam} -ref ${reference} -anc ${reference} -rf regions.txt \\
        -dosaf 1 -GL 2 \\
        -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \\
        -C 50 -baq 1 \\
        -P ${task.cpus} -out ${sample_id}

    # One diploid sample, reference as ancestral: 3 values = expected
    # number of sites with 0, 1, 2 non-reference alleles
    realSFS ${sample_id}.saf.idx -P ${task.cpus} > ${sample_id}.sfs
    """
}

process HETEROZYGOSITY_SUMMARY {
    tag "summary"
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: "heterozygosity_estimate.txt"

    input:
    path sfs_files          // <sample_id>.sfs, one per sampled sample
    val  indel_override     // --indel_heterozygosity, or '' to scale from SNPs

    output:
    path "heterozygosity_estimate.txt", emit: report
    path "gatk_prior_args.txt",         emit: args

    script:
    """
    python3 << 'PYEOF'
import glob, statistics

rows = []
for f in sorted(glob.glob('*.sfs')):
    vals = [float(x) for x in open(f).read().split()]
    if len(vals) != 3:
        raise SystemExit(f'ERROR: {f}: expected 3 SFS values for one diploid sample, got {len(vals)}')
    sites = sum(vals)
    rows.append((f[:-len('.sfs')], sites, vals[1] / sites if sites > 0 else 0.0))

median_sites = statistics.median(r[1] for r in rows)
used = [r for r in rows if r[1] > 0 and r[1] >= 0.5 * median_sites]
if not used:
    raise SystemExit('ERROR: no sample had enough covered sites to estimate heterozygosity '
                     '(single-end reads are dropped by the proper-pair filter)')

het = float(f'{statistics.median(r[2] for r in used):.3g}')
indel_override = "${indel_override}".strip()
indel = float(indel_override) if indel_override else float(f'{het * 0.125:.3g}')

with open('gatk_prior_args.txt', 'w') as f:
    f.write(f'--heterozygosity {het} --indel-heterozygosity {indel}')

with open('heterozygosity_estimate.txt', 'w') as f:
    f.write('# --heterozygosity auto (ANGSD per-sample SFS, ${params.angsd_container})\\n')
    f.write(f'heterozygosity (SNP prior)       {het}   median of {len(used)} sample(s)\\n')
    f.write(f'indel heterozygosity (prior)     {indel}   '
            f'{"from --indel_heterozygosity" if indel_override else "SNP x 0.125 (GATK default ratio)"}\\n')
    f.write(f'GenotypeGVCFs arguments          --heterozygosity {het} --indel-heterozygosity {indel}\\n\\n')
    f.write('sample\\tsites\\theterozygosity\\tused\\n')
    for name, sites, h in rows:
        f.write(f'{name}\\t{sites:.0f}\\t{h:.5f}\\t{"yes" if (name, sites, h) in used else "no (low coverage)"}\\n')
print(open('heterozygosity_estimate.txt').read())
PYEOF
    """
}
