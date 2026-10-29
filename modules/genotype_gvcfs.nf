process GENOTYPE_GVCFS {
    /*
     * Genotype a single GenomicsDB workspace.
     * Outputs one VCF per interval group — these are concatenated
     * downstream by GATHER_VCFS into per-chromosome VCFs.
     *
     * The --intervals flag restricts output to the primary (unpadded)
     * interval, discarding any boundary-padded calls that bled in
     * from HaplotypeCaller's --interval-padding 200.
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
    gatk GenotypeGVCFs \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --reference ${ref_abs} \\
        --variant gendb://${genomicsdb} \\
        --output ${group_id}.genotyped.vcf.gz \\
        --intervals ${interval_list} \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --tmp-dir .
    """

    stub:
    """
    touch ${group_id}.genotyped.vcf.gz
    touch ${group_id}.genotyped.vcf.gz.tbi
    """
}
