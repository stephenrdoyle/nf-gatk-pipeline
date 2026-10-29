process CUSTOM_FILTER {
    tag "${variant_type}"
    label 'process_low'

    conda 'bioconda::gatk4=4.6.1.0'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(variant_type), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(variant_type), path("${variant_type}.custom_filtered.vcf.gz"), path("${variant_type}.custom_filtered.vcf.gz.tbi"), emit: filtered_vcf
    path "${variant_type}_filter_summary.tsv", emit: summary

    script:
    // Heap tied to the LSF request: without -Xmx the JVM sizes its heap from
    // the node's RAM, not task.memory, and LSF kills the job (exit 130)
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    // params.snp_filters / params.indel_filters are comma-separated expressions
    // e.g. "QD < 2.0,FS > 60.0"
    def raw_filters = variant_type == 'snps' ? params.snp_filters : params.indel_filters
    def exprs = raw_filters.split(',').collect { it.trim() }.findAll { it }
    def filter_args = exprs.collect { expr ->
        def name = expr.replaceAll(/[^A-Za-z0-9_]/, '_').replaceAll(/_+/, '_')
        "--filter-expression \"${expr}\" --filter-name \"${name}\""
    }.join(' \\\n        ')
    // Threshold TSV rows: "ANN OP VALUE" -> annotation, operator, threshold
    def summary_rows = exprs.collect { expr ->
        def parts = expr.split(/\s+/)
        def ann = parts.size() > 0 ? parts[0] : 'NA'
        def op  = parts.size() > 1 ? parts[1] : 'NA'
        def thr = parts.size() > 2 ? parts[2] : 'NA'
        [ variant_type, ann, op, thr, 'NA', 'NA', 'NA', 'NA', 'NA' ].join('\t')
    }.join('\n')
    """
    gatk VariantFiltration \\
        --java-options "${java_opts}" \\
        --reference ${reference} \\
        --variant ${vcf} \\
        ${filter_args} \\
        --output ${variant_type}.custom_filtered.vcf.gz

    # Write threshold TSV for downstream summary
    cat > ${variant_type}_filter_summary.tsv << 'TSV'
variant_type\tannotation\toperator\tthreshold\tn_total\tn_flagged\tpct_flagged\tmean\tmedian
${summary_rows}
TSV
    """
}
