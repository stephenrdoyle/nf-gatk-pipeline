/*
 * DRAGstr model calibration for DRAGEN-GATK mode (--dragen_mode).
 *
 * DRAGEN-mode HaplotypeCaller models indel errors in short tandem repeats
 * (STRs) with per-sample parameters estimated from the sample's own reads:
 *   COMPOSE_STR_TABLE   once per reference: catalogue of STR sites to sample
 *   CALIBRATE_DRAGSTR   once per sample: fit the DRAGstr model on those sites
 * The resulting <sample>.dragstr file is passed to HAPLOTYPE_CALLER_DRAGEN.
 */

process COMPOSE_STR_TABLE {
    tag "reference"
    label 'process_low'

    input:
    path reference
    path reference_fai
    path reference_dict

    output:
    path "${reference.baseName}.str.zip", emit: str_table

    script:
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    def ref_abs   = reference.toRealPath()   // same reference form as every other GATK step
    """
    gatk ComposeSTRTableFile \\
        --java-options "${java_opts}" \\
        --reference ${ref_abs} \\
        --output ${reference.baseName}.str.zip
    """
}

process CALIBRATE_DRAGSTR {
    tag "${sample_id}"
    label 'process_low'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict
    path str_table

    output:
    tuple val(sample_id), path("${sample_id}.dragstr"), emit: params

    script:
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    def ref_abs   = reference.toRealPath()   // same reference form as every other GATK step
    """
    gatk CalibrateDragstrModel \\
        --java-options "${java_opts}" \\
        --reference ${ref_abs} \\
        --input ${bam} \\
        --str-table-path ${str_table} \\
        --output ${sample_id}.dragstr
    """
}
