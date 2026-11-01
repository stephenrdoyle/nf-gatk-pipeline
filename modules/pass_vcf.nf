process PASS_VCF {
    /*
     * --write_pass_vcf: the final VCF with failing sites removed, for tools
     * or collaborators that don't read the FILTER column. Off by default:
     * it is another copy of most of the VCF, and bcftools view -f PASS on the
     * final VCF gives the same sites.
     */
    tag "pass_vcf"

    publishDir "${params.outdir}/vcf", mode: params.publish_mode,
        saveAs: { fn -> fn.replace('final_PASS', "${params.cohort}.${params.pooled ? '3_pooled_final' : '2_final'}_PASS") }

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("final_PASS.vcf.gz"), path("final_PASS.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view --threads ${task.cpus} --apply-filters 'PASS,.' -Oz -o final_PASS.vcf.gz ${vcf}
    bcftools index --tbi --threads ${task.cpus} final_PASS.vcf.gz
    """

    stub:
    """
    touch final_PASS.vcf.gz final_PASS.vcf.gz.tbi
    """
}
