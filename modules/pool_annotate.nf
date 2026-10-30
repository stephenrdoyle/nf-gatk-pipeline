process POOL_ANNOTATE {
    tag "pool_annotate"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("cohort.pooled_annotated.vcf.gz"), path("cohort.pooled_annotated.vcf.gz.tbi"), emit: annotated_vcf

    script:
    // GT is set to missing for every sample (diploid calls are meaningless for
    // pools); FORMAT/AF per ALT allele is added from AD, '.' below this depth
    def min_dp = params.genotype_dp_threshold ?: 0
    """
    # Stream straight into BGZF: no uncompressed whole-cohort intermediate.
    # pipefail: Nextflow runs bash -ue only, so without it a crash in the
    # Python step would leave a truncated VCF and a task that "succeeded"
    set -o pipefail
    pool_annotate.py ${vcf} - ${min_dp} \\
        | bcftools view --threads ${task.cpus} -O z -o cohort.pooled_annotated.vcf.gz -
    bcftools index --tbi cohort.pooled_annotated.vcf.gz
    """
}
