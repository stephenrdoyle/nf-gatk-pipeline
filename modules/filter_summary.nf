process FILTER_SUMMARY {
    tag "summary"
    label 'process_low'

    publishDir "${params.outdir}/filter_reports", mode: 'copy'

    conda 'bioconda::bcftools=1.19 conda-forge::python=3.11'
    container 'biocontainers/bcftools:1.19--h8b25389_0'

    input:
    // Raw (pre-filter) SNP and indel VCFs
    tuple path(raw_snp_vcf),      path(raw_snp_tbi)
    tuple path(raw_indel_vcf),    path(raw_indel_tbi)
    // Filtered SNP and indel VCFs
    tuple path(filt_snp_vcf),     path(filt_snp_tbi)
    tuple path(filt_indel_vcf),   path(filt_indel_tbi)
    // Optional per-annotation threshold TSVs (empty list [] if not optimised mode)
    path snp_summary_tsv
    path indel_summary_tsv

    output:
    path "filter_summary.tsv", emit: summary

    script:
    def snp_tsv_arg   = snp_summary_tsv   instanceof List ? '' : snp_summary_tsv.name
    def indel_tsv_arg = indel_summary_tsv instanceof List ? '' : indel_summary_tsv.name
    """
    count_filter_stats.py ${raw_snp_vcf}   ${filt_snp_vcf}   snps   ${snp_tsv_arg}
    count_filter_stats.py ${raw_indel_vcf} ${filt_indel_vcf} indels ${indel_tsv_arg}

    # Combine SNP and indel counts into one TSV
    head -1 snps_counts.tsv > filter_summary.tsv
    tail -n +2 snps_counts.tsv   >> filter_summary.tsv
    tail -n +2 indels_counts.tsv >> filter_summary.tsv

    echo "Combined filter summary written: filter_summary.tsv"
    """
}
