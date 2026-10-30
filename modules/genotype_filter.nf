process GENOTYPE_FILTER {
    /*
     * Per-genotype filter on the whole cohort VCF: genotypes with low depth
     * (DP) or low genotype quality (GQ / RGQ) are set to missing (./.).
     *
     * bcftools +setGT does this in one pass with multithreaded compression.
     * GATK VariantFiltration evaluated a JEXL expression for every sample at
     * every site, single-threaded: 7.6 h for the 200-sample F2 cohort.
     *
     * Same rules as before:
     *   DP  < --genotype_dp_threshold    always
     *   GQ  < --genotype_gq_threshold    not with --pooled (GQ is meaningless
     *                                    for pools) or --use_gnarly (no GQ)
     *   RGQ < --genotype_rgq_threshold   --use_gnarly only, if set
     * A missing value (e.g. DP '.') never fails a filter, as in GATK.
     *
     * Unlike GATK, failed genotypes get no FORMAT/FT tag; SAMPLE_MISSINGNESS
     * counts low-depth genotypes from FORMAT/DP directly. The site FILTER
     * column is left as it was ('.'), and the site filters set PASS later.
     */
    tag "genotype_filter"
    label 'process_low'

    publishDir "${params.outdir}", mode: params.publish_mode, pattern: "cohort.genotype_filtered.vcf.gz*"

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("cohort.genotype_filtered.vcf.gz"), path("cohort.genotype_filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    def conditions = []
    if (params.genotype_dp_threshold)
        conditions << "FMT/DP<${params.genotype_dp_threshold}"
    if (!params.pooled && !params.use_gnarly && params.genotype_gq_threshold)
        conditions << "FMT/GQ<${params.genotype_gq_threshold}"
    if (params.use_gnarly && !params.pooled && params.genotype_rgq_threshold)
        conditions << "FMT/RGQ<${params.genotype_rgq_threshold}"
    def set_gt = conditions
        ? "bcftools +setGT --threads ${task.cpus} -Oz -o cohort.genotype_filtered.vcf.gz ${vcf} -- -t q -n . -i '${conditions.join(' | ')}'"
        : "bcftools view --threads ${task.cpus} -Oz -o cohort.genotype_filtered.vcf.gz ${vcf}"   // no thresholds set
    """
    ${set_gt}
    bcftools index --tbi --threads ${task.cpus} cohort.genotype_filtered.vcf.gz
    """
}
