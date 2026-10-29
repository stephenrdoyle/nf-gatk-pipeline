process BCFTOOLS_STATS {
    tag "${label}"
    label 'process_low'

    publishDir "${params.outdir}/bcftools_stats", mode: 'copy'

    conda 'bioconda::bcftools=1.19'
    container 'biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple path(vcf), path(tbi)
    val   label   // e.g. "before_filter", "after_filter", "after_genotype_filter"

    output:
    path "${label}_bcftools_stats.txt", emit: stats

    script:
    """
    bcftools stats -s - ${vcf} > ${label}_bcftools_stats.txt
    """
}
