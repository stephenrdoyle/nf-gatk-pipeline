process MERGE_FILTERED {
    tag "merge"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(snp_vcf),   path(snp_tbi)
    tuple path(indel_vcf), path(indel_tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    // Heap tied to the LSF request: without -Xmx the JVM sizes its heap from
    // the node's RAM, not task.memory, and LSF kills the job (exit 130)
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    """
    gatk MergeVcfs \\
        --java-options "${java_opts}" \\
        --INPUT ${snp_vcf} \\
        --INPUT ${indel_vcf} \\
        --OUTPUT cohort.filtered.vcf.gz
    """
}
