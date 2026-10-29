process GNARLY_GENOTYPER {
    /*
     * GnarlyGenotyper: faster alternative to GenotypeGVCFs for large cohorts.
     * Incompatible with --dragen_mode (guard enforced in main.nf).
     *
     * Updated to use path(interval_list) to match the interval-scatter design,
     * replacing the old val(interval_arg) pattern.
     * Emit name changed to 'vcf' to match GENOTYPE_GVCFS — both branches
     * in main.nf feed into the same genotyped_vcfs channel.
     */

    tag { group_id }
    label 'process_high'

    input:
    tuple val(group_id), path(genomicsdb), path(interval_list)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(group_id), path("${group_id}.genotyped.vcf.gz"), path("${group_id}.genotyped.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = task.memory.toGiga() - 8
    def ref_abs = reference.toRealPath()
    """
    gatk GnarlyGenotyper \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --reference ${ref_abs} \\
        --variant gendb://${genomicsdb} \\
        --output ${group_id}.genotyped.vcf.gz \\
        --intervals ${interval_list} \\
        --only-output-calls-starting-in-intervals \\
        --tmp-dir .
    """

    stub:
    """
    touch ${group_id}.genotyped.vcf.gz
    touch ${group_id}.genotyped.vcf.gz.tbi
    """
}
