process POOL_ANNOTATE {
    tag "pool_annotate"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("cohort.pooled_annotated.vcf.gz"), path("cohort.pooled_annotated.vcf.gz.tbi"), emit: annotated_vcf

    script:
    def het_lo     = params.pool_het_lo
    def het_hi     = params.pool_het_hi
    def hom_alt_lo = params.pool_hom_alt_lo
    def min_dp     = params.genotype_dp_threshold ?: 0
    """
    # Get all sample IDs from the VCF header
    sample_id=\$(bcftools query -l ${vcf} | tr '\\n' ',' | sed 's/,\$//')

    # Write uncompressed VCF — Python gzip is not BGZF and cannot be tabix indexed
    pool_annotate.py \\
        ${vcf} \\
        cohort.pooled_annotated.vcf \\
        "\${sample_id}" \\
        ${het_lo} \\
        ${het_hi} \\
        ${hom_alt_lo} \\
        ${min_dp}

    # Compress as BGZF and index
    bcftools view -O z -o cohort.pooled_annotated.vcf.gz cohort.pooled_annotated.vcf
    bcftools index --tbi cohort.pooled_annotated.vcf.gz
    rm cohort.pooled_annotated.vcf
    """
}
