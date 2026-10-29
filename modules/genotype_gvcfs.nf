process GENOTYPE_GVCFS {
    /*
     * Genotype a single GenomicsDB workspace.
     * Outputs one VCF per interval group — these are concatenated
     * downstream by GATHER_VCFS into the cohort VCF.
     *
     * --intervals restricts genotyping to the primary (unpadded) interval.
     * On its own that still emits any record that merely overlaps the
     * interval, so a deletion starting just before a window boundary is
     * written by both neighbouring groups and duplicated after GATHER_VCFS
     * (bcftools concat does not deduplicate).
     * --only-output-calls-starting-in-intervals keeps each record only in
     * the group where it starts, so every site is emitted exactly once.
     * GNARLY_GENOTYPER uses the same flag.
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
        --only-output-calls-starting-in-intervals \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --tmp-dir .
    """

    stub:
    """
    touch ${group_id}.genotyped.vcf.gz
    touch ${group_id}.genotyped.vcf.gz.tbi
    """
}
