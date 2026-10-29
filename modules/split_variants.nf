process SPLIT_VARIANTS {
    tag "split"
    label 'process_low'

    conda 'bioconda::gatk4=4.6.1.0'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple path("snps.vcf.gz"),   path("snps.vcf.gz.tbi"),   emit: snps
    tuple path("indels.vcf.gz"), path("indels.vcf.gz.tbi"), emit: indels

    script:
    // Heap tied to the LSF request: without -Xmx the JVM sizes its heap from
    // the node's RAM, not task.memory, and LSF kills the job (exit 130)
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    // MIXED records (a site with both SNP and indel alternate alleles) go
    // with the indels, as in GATK's hard-filtering guidance: they are kept
    // intact as multiallelic records and filtered with the indel thresholds.
    // Selecting only SNP and INDEL dropped them from the final VCF entirely.
    """
    gatk SelectVariants \\
        --java-options "${java_opts}" \\
        --reference ${reference} \\
        --variant ${vcf} \\
        --select-type-to-include SNP \\
        --output snps.vcf.gz

    gatk SelectVariants \\
        --java-options "${java_opts}" \\
        --reference ${reference} \\
        --variant ${vcf} \\
        --select-type-to-include INDEL \\
        --select-type-to-include MIXED \\
        --output indels.vcf.gz
    """
}
