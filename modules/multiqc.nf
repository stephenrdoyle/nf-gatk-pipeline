process MULTIQC {
    tag "multiqc"
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path collected_files

    output:
    path "multiqc_report.html",    emit: report
    path "multiqc_report_data/",   emit: data

    script:
    def config_arg = file("${projectDir}/assets/multiqc_config.yml").exists()
        ? "--config ${projectDir}/assets/multiqc_config.yml"
        : ''
    """
    # multiqc comes from its own module (see MULTIQC in nextflow.config).
    # Unset PYTHONPATH in case a python module's numpy is still inherited
    # from the submitting environment and clashes with MultiQC's
    unset PYTHONPATH

    multiqc . \\
        --filename multiqc_report.html \\
        ${config_arg} \\
        --force
    """
}
