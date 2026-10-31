/*
 * End-of-run cohort QC: one self-contained HTML report (qc_report.html).
 *
 *   QC_EXTRACT    bcftools, on the filtered cohort VCF (PASS sites):
 *                 the SNP sets for SNPRelate and the per-sample stats the
 *                 report reads
 *   QC_SNPRELATE  SNPRelate (container): nuclear PCA, KING kinship and,
 *                 with --mito_chrom, a mitochondrial PCA
 *   QC_REPORT     bin/qc_report.py: the HTML report
 *
 * --sex_chrom and --mito_chrom are optional. Set, their contigs are left out
 * of the nuclear SNP set and stats, and get their own report sections
 * (chrX heterozygosity vs the samplesheet's sex column; mito PCA and
 * heteroplasmy). Unset, every contig counts as nuclear.
 */

process QC_EXTRACT {
    tag "qc"
    label 'process_medium'

    input:
    tuple path(vcf), path(tbi)
    val   sex_chrom    // contig name or ''
    val   mito_chrom   // contig name or ''
    val   maf

    output:
    path "nuclear_snps.vcf.gz",     emit: nuclear
    path "mito_snps.*",             emit: mito      // mito_snps.vcf.gz, or mito_snps.NONE
    path "stats_*.txt",             emit: stats
    path "density.tsv",             emit: density

    script:
    def exclude = [sex_chrom, mito_chrom].findAll { it }
    def nuclear_targets = exclude ? "--targets '^${exclude.join(',')}'" : ''
    // PASS biallelic SNPs, common enough to carry structure. +prune keeps the
    // first SNP per kb, a cheap thinning before SNPRelate's LD pruning
    def snp_select = "--apply-filters 'PASS,.' --min-alleles 2 --max-alleles 2 --types snps --min-af ${maf}:minor"
    def sex_stats  = sex_chrom
        ? "bcftools stats --regions '${sex_chrom}' --apply-filters 'PASS,.' --samples - ${vcf} > stats_sex.txt"
        : ''
    """
    set -o pipefail

    # Full passes over the VCF run in parallel; each is waited on by PID so a
    # failure in any of them fails the task
    bcftools view ${snp_select} ${nuclear_targets} -Ou ${vcf} \\
        | bcftools +prune -w 1kb -n 1 -N 1st -Oz -o nuclear_snps.vcf.gz - &
    p1=\$!
    bcftools stats ${nuclear_targets} --apply-filters 'PASS,.' --samples - ${vcf} > stats_nuclear.txt &
    p2=\$!
    # PASS sites per Mb, contigs in VCF order: contig, mb, n
    bcftools view --drop-genotypes --no-header --apply-filters 'PASS,.' ${vcf} \\
        | awk -F'\\t' -v OFS='\\t' '
            BEGIN { print "contig", "mb", "n" }
            { k = \$1 OFS int((\$2 - 1) / 1000000); if (!(k in n)) order[++m] = k; n[k]++ }
            END { for (i = 1; i <= m; i++) print order[i], n[order[i]] }' > density.tsv &
    p3=\$!

    # Single contigs: region queries through the index, quick
    ${sex_stats}
    if [ -n '${mito_chrom}' ]; then
        bcftools stats --regions '${mito_chrom}' --apply-filters 'PASS,.' --samples - ${vcf} > stats_mito.txt
        # Every polymorphic mito SNP (no MAF cut-off: lineages can be rare)
        bcftools view --regions '${mito_chrom}' --apply-filters 'PASS,.' --min-alleles 2 --max-alleles 2 \\
            --types snps -Oz -o mito_snps.vcf.gz ${vcf}
    else
        touch mito_snps.NONE
    fi

    wait \$p1
    wait \$p2
    wait \$p3
    echo "Nuclear QC SNPs: \$(bcftools view --no-header nuclear_snps.vcf.gz | grep -c . || true)"
    """

    stub:
    """
    touch nuclear_snps.vcf.gz mito_snps.NONE stats_nuclear.txt density.tsv
    """
}

process QC_SNPRELATE {
    // Container (params.snprelate_container) is set in nextflow.config
    tag "qc"

    input:
    path nuclear
    path mito
    val  ld_threshold

    output:
    path "pca_*.tsv",       emit: pca
    path "kinship.tsv",     emit: kinship
    path "snp_counts.tsv",  emit: counts
    path "notes.txt",       emit: notes

    script:
    """
    qc_snprelate.R ${nuclear} ${mito} ${ld_threshold} ${task.cpus}
    """

    stub:
    """
    touch pca_nuclear.tsv pca_nuclear_eigen.tsv kinship.tsv snp_counts.tsv notes.txt
    """
}

process QC_REPORT {
    tag "qc"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path files, stageAs: 'inputs/*'
    val  samples_tsv   // sample, sex, group (from the samplesheet)
    val  meta_json     // run metadata and parameters

    output:
    path "qc_report.html", emit: report

    script:
    """
    cat > samples.tsv << 'SAMPLES'
${samples_tsv}
SAMPLES
    cat > meta.json << 'META'
${meta_json}
META
    qc_report.py --inputs inputs --samples samples.tsv --meta meta.json --out qc_report.html
    """

    stub:
    """
    touch qc_report.html
    """
}
