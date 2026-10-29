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
    // params.snp_filters / params.indel_filters are comma-separated expressions
    // e.g. "QD < 2.0,FS > 60.0"
    def raw_filters = variant_type == 'snps' ? params.snp_filters : params.indel_filters
    def filter_args = raw_filters.split(',').collect { expr ->
        def name = expr.trim().replaceAll(/[^A-Za-z0-9_]/, '_').replaceAll(/_+/, '_')
        "--filter-expression \"${expr.trim()}\" --filter-name \"${name}\""
    }.join(' \\\n        ')
    """
    gatk VariantFiltration \\
        --reference ${reference} \\
        --variant ${vcf} \\
        ${filter_args} \\
        --output ${variant_type}.custom_filtered.vcf.gz

    # Write threshold TSV for downstream summary
    python3 -c "
import sys
vtype    = '${variant_type}'
raw      = params_filters = '${variant_type == "snps" ? params.snp_filters : params.indel_filters}'
filters  = [e.strip() for e in raw.split(',')]
with open(f'{vtype}_filter_summary.tsv', 'w') as f:
    f.write('variant_type	annotation	operator	threshold	n_total	n_flagged	pct_flagged	mean	median
')
    for expr in filters:
        parts = expr.split()
        ann   = parts[0] if len(parts) > 0 else 'NA'
        op    = parts[1] if len(parts) > 1 else 'NA'
        val   = parts[2] if len(parts) > 2 else 'NA'
        f.write(f'{vtype}	{ann}	{op}	{val}	NA	NA	NA	NA	NA
')
"
    """
}
