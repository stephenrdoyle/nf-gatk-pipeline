process SAMTOOLS_FAIDX {
    tag "faidx"
    label 'process_low'

    // Publish next to the reference so GATK can find it. overwrite: false:
    // every new launch directory reruns this, and rewriting the shared file
    // could briefly disturb other runs reading it. A stale .fai (older than
    // the FASTA) is caught at launch in main.nf instead.
    publishDir "${file(params.reference).parent}", mode: 'copy', overwrite: false

    input:
    path reference

    output:
    path "${reference}.fai", emit: fai

    script:
    """
    samtools faidx ${reference}
    """
}
