process FILTER_SUMMARY {
    /*
     * Sites per filter for the filtered SNP and indel VCFs, read from their
     * FILTER columns and ##FILTER headers (see bin/count_filter_stats.py), so
     * every filter mode gets real per-filter counts.
     */
    tag "summary"
    label 'process_low'

    publishDir "${params.outdir}/filter_reports", mode: 'copy'

    input:
    tuple path(snp_vcf),   path(snp_tbi)
    tuple path(indel_vcf), path(indel_tbi)

    output:
    path "filter_summary.tsv", emit: summary

    script:
    """
    count_filter_stats.py snps=${snp_vcf} indels=${indel_vcf} > filter_summary.tsv
    """
}
