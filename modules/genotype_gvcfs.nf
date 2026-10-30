process GENOTYPE_GVCFS {
    /*
     * Genotype one joint-genotyping piece's GenomicsDB workspace.
     * Outputs one VCF per piece; GATHER_VCFS concatenates them, in
     * part_id order, into the cohort VCF.
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
    gatk GenotypeGVCFs \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --reference ${ref_abs} \\
        --variant gendb://${genomicsdb} \\
        --output ${part_id}.genotyped.vcf.gz \\
        --intervals ${interval} \\
        --only-output-calls-starting-in-intervals \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --tmp-dir .
    """

    stub:
    """
    touch ${part_id}.genotyped.vcf.gz
    touch ${part_id}.genotyped.vcf.gz.tbi
    """
}
