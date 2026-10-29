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
    // Same reference form as HAPLOTYPE_CALLER / GENOTYPE_GVCFS: the real path,
    // with symlinks (e.g. lustre_link) resolved. GATK runs in a Singularity
    // wrapper, so a resolved /lustre path doesn't depend on the container
    // also seeing a staged symlink's target; .fai/.dict are published next
    // to the real FASTA by SAMTOOLS_FAIDX / CREATE_DICT
    def ref_abs = reference.toRealPath()
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
        --reference ${ref_abs} \\
        --variant ${vcf} \\
        --select-type-to-include SNP \\
        --output snps.vcf.gz

    gatk SelectVariants \\
        --java-options "${java_opts}" \\
        --reference ${ref_abs} \\
        --variant ${vcf} \\
        --select-type-to-include INDEL \\
        --select-type-to-include MIXED \\
        --output indels.vcf.gz
    """
}
