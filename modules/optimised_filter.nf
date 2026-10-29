process OPTIMISED_FILTER {
    tag "${variant_type}"
    label 'process_medium'

    publishDir "${params.outdir}/filter_reports", mode: 'copy', pattern: "*.html"

    conda 'bioconda::gatk4=4.6.1.0 conda-forge::python=3.11'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(variant_type), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(variant_type), path("${variant_type}.optimised_filtered.vcf.gz"), path("${variant_type}.optimised_filtered.vcf.gz.tbi"), emit: filtered_vcf
    path "${variant_type}_filter_report.html",       emit: report
    path "${variant_type}_filter_summary.tsv",        emit: summary

    script:
    // Same reference form as HAPLOTYPE_CALLER / GENOTYPE_GVCFS: the real path,
    // with symlinks (e.g. lustre_link) resolved. GATK runs in a Singularity
    // wrapper, so a resolved /lustre path doesn't depend on the container
    // also seeing a staged symlink's target; .fai/.dict are published next
    // to the real FASTA by SAMTOOLS_FAIDX / CREATE_DICT
    def ref_abs = reference.toRealPath()
    // Heap tied to the LSF request: without -Xmx the JVM sizes its heap from
    // the node's RAM, not task.memory, and LSF kills the job (exit 130)
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    def lo_pct  = params.filter_percentile_lo as int
    def ann_str = variant_type == 'snps'
                  ? 'QD,MQ,FS,SOR,MQRankSum,ReadPosRankSum'
                  : 'QD,FS,SOR,ReadPosRankSum'
    def lo_bad  = 'QD,MQ,MQRankSum,ReadPosRankSum'
    def hi_bad  = 'FS,SOR'
    def gatk_f  = ann_str.split(',').collect { "-F ${it}" }.join(' ')
    """
    # Extract annotation values from VCF
    gatk VariantsToTable \\
        --java-options "${java_opts}" \\
        --reference ${ref_abs} \\
        --variant ${vcf} \\
        ${gatk_f} \\
        --output annotations.tsv

    # Compute thresholds, write filter_args.txt and HTML report
    optimised_filter.py \\
        annotations.tsv \\
        ${variant_type} \\
        ${lo_pct} \\
        ${ann_str} \\
        ${lo_bad} \\
        ${hi_bad}

    # Apply the derived filters — source the array script to safely handle
    # negative numbers and spaces in filter expressions
    source filter_args.sh
    gatk VariantFiltration \\
        --java-options "${java_opts}" \\
        --reference ${ref_abs} \\
        --variant ${vcf} \\
        "\${FILTER_ARGS[@]}" \\
        --output ${variant_type}.optimised_filtered.vcf.gz
    """
}
