process MERGE_FILTERED {
    tag "merge"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    conda 'bioconda::gatk4=4.6.1.0'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple path(snp_vcf),   path(snp_tbi)
    tuple path(indel_vcf), path(indel_tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    """
    gatk MergeVcfs \\
        --INPUT ${snp_vcf} \\
        --INPUT ${indel_vcf} \\
        --OUTPUT cohort.filtered.vcf.gz
    """
}
