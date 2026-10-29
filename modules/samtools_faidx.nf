process SAMTOOLS_FAIDX {
    tag "faidx"
    label 'process_low'

    // Publish permanently next to the reference so GATK can find it
    publishDir "${file(params.reference).parent}", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}.fai", emit: fai

    script:
    """
    samtools faidx ${reference}
    """
}
