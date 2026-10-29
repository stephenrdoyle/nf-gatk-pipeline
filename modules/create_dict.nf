process CREATE_DICT {
    /*
     * GATK requires the .dict to exist next to the reference .fa.
     * publishDir writes it there permanently so GATK's engine-level
     * validation always finds it, regardless of which tool is called.
     * Once present, Nextflow's cache means this process never reruns.
     */
    tag "reference"
    label 'process_low'

    publishDir "${file(params.reference).parent}", mode: 'copy'

    input:
    path reference

    output:
    path "${reference.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${reference}
    """
}
