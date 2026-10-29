process SAMTOOLS_INDEX {
    tag "${sample_id}"
    label 'process_low'

    conda 'bioconda::samtools=1.19'
    container 'biocontainers/samtools:1.19--h50ea8bc_1'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}
