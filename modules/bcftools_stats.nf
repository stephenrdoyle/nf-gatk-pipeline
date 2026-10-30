process BCFTOOLS_STATS {
    /*
     * bcftools stats for one stage of the cohort VCF.
     *
     * Site filters only flag records in the FILTER column; they never remove
     * them, and bcftools stats counts every record unless told otherwise. So
     * the after-filter stats use pass_only; without it they describe failed
     * sites too and look almost identical to the before-filter stats.
     */
    tag "${label}"
    label 'process_low'

    publishDir "${params.outdir}/bcftools_stats", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    val   label       // "before_filter" or "after_filter"
    val   pass_only   // true: count only sites whose FILTER is PASS (or '.')

    output:
    path "${label}_bcftools_stats.txt", emit: stats

    script:
    def filter_arg = pass_only ? "--apply-filters 'PASS,.'" : ''
    """
    bcftools stats --threads ${task.cpus} ${filter_arg} --samples - ${vcf} > ${label}_bcftools_stats.txt
    """
}
