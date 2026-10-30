process GNARLY_GENOTYPER {
    /*
     * GnarlyGenotyper: faster alternative to GenotypeGVCFs for large cohorts.
     * Incompatible with --dragen_mode (guard enforced in main.nf).
     *
     * Same inputs and outputs as GENOTYPE_GVCFS (one joint-genotyping piece,
     * "contig:start-end"), so both feed the same genotyped_vcfs channel.
     */

    tag { part_id }
    label 'process_high'

    input:
    tuple val(part_id), path(genomicsdb), val(interval)   // interval: "contig:start-end"
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(part_id), path("${part_id}.genotyped.vcf.gz"), path("${part_id}.genotyped.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = task.memory.toGiga() - 8
    def ref_abs = reference.toRealPath()
    """
    gatk GnarlyGenotyper \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --reference ${ref_abs} \\
        --variant gendb://${genomicsdb} \\
        --output ${part_id}.genotyped.vcf.gz \\
        --intervals ${interval} \\
        --only-output-calls-starting-in-intervals \\
        --tmp-dir .
    """

    stub:
    """
    touch ${part_id}.genotyped.vcf.gz
    touch ${part_id}.genotyped.vcf.gz.tbi
    """
}
