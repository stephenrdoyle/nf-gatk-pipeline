process HARD_FILTER {
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
    tuple val(variant_type), path("${variant_type}.hard_filtered.vcf.gz"), path("${variant_type}.hard_filtered.vcf.gz.tbi"), emit: filtered_vcf
    path "${variant_type}_filter_summary.tsv", emit: summary

    script:
    def filters = variant_type == 'snps' ? """
        --filter-expression "QD < 2.0"               --filter-name "SNP_QD2" \\
        --filter-expression "MQ < 40.0"              --filter-name "SNP_MQ40" \\
        --filter-expression "FS > 60.0"              --filter-name "SNP_FS60" \\
        --filter-expression "SOR > 3.0"              --filter-name "SNP_SOR3" \\
        --filter-expression "MQRankSum < -12.5"      --filter-name "SNP_MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0"  --filter-name "SNP_ReadPosRankSum-8"
    """ : """
        --filter-expression "QD < 2.0"               --filter-name "INDEL_QD2" \\
        --filter-expression "FS > 200.0"             --filter-name "INDEL_FS200" \\
        --filter-expression "SOR > 10.0"             --filter-name "INDEL_SOR10" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "INDEL_ReadPosRankSum-20"
    """
    """
    gatk VariantFiltration \\
        --reference ${reference} \\
        --variant ${vcf} \\
        ${filters} \\
        --output ${variant_type}.hard_filtered.vcf.gz

    # Write threshold TSV for downstream summary
    python3 -c "
import sys
vtype = '${variant_type}'
snp_filters = [
    ('QD','lt',2.0),('MQ','lt',40.0),('FS','gt',60.0),
    ('SOR','gt',3.0),('MQRankSum','lt',-12.5),('ReadPosRankSum','lt',-8.0)
]
indel_filters = [
    ('QD','lt',2.0),('FS','gt',200.0),
    ('SOR','gt',10.0),('ReadPosRankSum','lt',-20.0)
]
filters = snp_filters if vtype == 'snps' else indel_filters
with open(f'{vtype}_filter_summary.tsv', 'w') as f:
    f.write('variant_type	annotation	operator	threshold	n_total	n_flagged	pct_flagged	mean	median
')
    for ann, direction, val in filters:
        op = '<' if direction == 'lt' else '>'
        f.write(f'{vtype}	{ann}	{op}	{val}	NA	NA	NA	NA	NA
')
"
    """
}
