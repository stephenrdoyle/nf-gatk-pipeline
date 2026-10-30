process HAPLOTYPE_CALLER {
    /*
     * Per-sample, per-interval HaplotypeCaller.
     * --interval-padding 200 ensures variants near interval boundaries
     * are called with full haplotype context; the padded bases are
     * trimmed back to primary interval boundaries by GenotypeGVCFs.
     *
     * Fixes applied:
     *   - Intel GKL deflater disabled (SIGBUS on farm22)
     *   - JVM heap headroom: reserve 4 GB for native allocations
     *   - Explicit -Xms4g to avoid slow JVM startup on large heaps
     */

    tag { "${sample_id}:${interval_list.baseName}" }
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), path(bai)
    each path(interval_list)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(sample_id), val("${interval_list.baseName}"), path("${sample_id}.${interval_list.baseName}.g.vcf.gz"), path("${sample_id}.${interval_list.baseName}.g.vcf.gz.tbi"), emit: gvcfs

    script:
    def avail_mem = task.memory.toGiga() - 4
    def ref_abs = reference.toRealPath()
    """
    gatk HaplotypeCaller \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --reference ${ref_abs} \\
        --input ${bam} \\
        --output ${sample_id}.${interval_list.baseName}.g.vcf.gz \\
        --intervals ${interval_list} \\
        --interval-padding 200 \\
        --emit-ref-confidence GVCF \\
        --tmp-dir .
    """

    stub:
    """
    touch ${sample_id}.${interval_list.baseName}.g.vcf.gz
    touch ${sample_id}.${interval_list.baseName}.g.vcf.gz.tbi
    """
}

process HAPLOTYPE_CALLER_DRAGEN {
    /*
     * DRAGEN-GATK variant of HAPLOTYPE_CALLER (--dragen_mode): identical
     * apart from --dragen-mode and the sample's calibrated DRAGstr model.
     *
     * Kept as a separate process, rather than a flag on HAPLOTYPE_CALLER,
     * because Nextflow's cache key includes the script block's source and
     * the params it reads: any edit to HAPLOTYPE_CALLER would make -resume
     * rerun every per-sample, per-window job of existing runs.
     */

    tag { "${sample_id}:${interval_list.baseName}" }
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), path(bai), path(dragstr_params)
    each path(interval_list)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(sample_id), val("${interval_list.baseName}"), path("${sample_id}.${interval_list.baseName}.g.vcf.gz"), path("${sample_id}.${interval_list.baseName}.g.vcf.gz.tbi"), emit: gvcfs

    script:
    def avail_mem = task.memory.toGiga() - 4
    def ref_abs = reference.toRealPath()
    """
    gatk HaplotypeCaller \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --reference ${ref_abs} \\
        --input ${bam} \\
        --output ${sample_id}.${interval_list.baseName}.g.vcf.gz \\
        --intervals ${interval_list} \\
        --interval-padding 200 \\
        --emit-ref-confidence GVCF \\
        --dragen-mode true \\
        --dragstr-params-path ${dragstr_params} \\
        --tmp-dir .
    """

    stub:
    """
    touch ${sample_id}.${interval_list.baseName}.g.vcf.gz
    touch ${sample_id}.${interval_list.baseName}.g.vcf.gz.tbi
    """
}
