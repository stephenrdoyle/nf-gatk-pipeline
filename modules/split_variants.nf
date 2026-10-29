process SPLIT_VARIANTS {
    tag "split"
    label 'process_low'

    input:
    tuple path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple path("snps.vcf.gz"),   path("snps.vcf.gz.tbi"),   emit: snps
    tuple path("indels.vcf.gz"), path("indels.vcf.gz.tbi"), emit: indels

    script:
    """
    gatk SelectVariants \\
        --reference ${reference} \\
        --variant ${vcf} \\
        --select-type-to-include SNP \\
        --output snps.vcf.gz

    gatk SelectVariants \\
        --reference ${reference} \\
        --variant ${vcf} \\
        --select-type-to-include INDEL \\
        --output indels.vcf.gz
    """
}
