process GENOTYPE_FILTER {
    tag "genotype_filter"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy', pattern: "cohort.genotype_filtered.vcf.gz*"

    input:
    tuple path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple path("cohort.genotype_filtered.vcf.gz"), path("cohort.genotype_filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    // GQ filter — skipped for pooled samples (GQ not meaningful for pools)
    //           — skipped for GnarlyGenotyper output (GQ not emitted; use RGQ instead)
    def gq_filter = (!params.pooled && !params.use_gnarly && params.genotype_gq_threshold)
        ? "--genotype-filter-expression \"GQ < ${params.genotype_gq_threshold}\" --genotype-filter-name \"LowGQ\""
        : ''

    // RGQ filter — only applied when using GnarlyGenotyper (--use_gnarly)
    // RGQ = Reference Genotype Quality: confidence a site is hom-ref
    // Analogous to GQ but specific to GnarlyGenotyper output
    def rgq_filter = (params.use_gnarly && !params.pooled && params.genotype_rgq_threshold)
        ? "--genotype-filter-expression \"RGQ < ${params.genotype_rgq_threshold}\" --genotype-filter-name \"LowRGQ\""
        : ''

    // DP filter — applied regardless of mode
    def dp_filter = params.genotype_dp_threshold
        ? "--genotype-filter-expression \"DP < ${params.genotype_dp_threshold}\" --genotype-filter-name \"LowDP\""
        : ''
    """
    gatk VariantFiltration \\
        --reference ${reference} \\
        --variant ${vcf} \\
        ${gq_filter} \\
        ${rgq_filter} \\
        ${dp_filter} \\
        --set-filtered-genotype-to-no-call true \\
        --output cohort.genotype_filtered.vcf.gz
    """
}
