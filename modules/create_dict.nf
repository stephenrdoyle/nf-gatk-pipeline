process CREATE_DICT {
    /*
     * GATK requires the .dict to exist next to the reference .fa.
     * publishDir writes it there permanently so GATK's engine-level
     * validation always finds it, regardless of which tool is called.
     *
     * The cache only spans one launch directory's runs, so every new
     * project reruns this. overwrite: false stops it rewriting a .dict that
     * other runs may be reading; a stale .dict (older than the FASTA) is
     * caught at launch in main.nf instead.
     */
    tag "reference"
    label 'process_low'

    publishDir "${file(params.reference).parent}", mode: 'copy', overwrite: false

    input:
    path reference

    output:
    path "${reference.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${reference}
    """
}
