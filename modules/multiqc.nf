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
    # Unset PYTHONPATH to prevent bcftools-1.19/python-3.11.6 numpy
    # from conflicting with MultiQC's Python 3.12 numpy
    unset PYTHONPATH

    multiqc . \\
        --filename multiqc_report.html \\
        ${config_arg} \\
        --force
    """
}
