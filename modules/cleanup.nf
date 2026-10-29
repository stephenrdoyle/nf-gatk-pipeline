process CLEANUP {
    /*
     * Remove intermediate VCF work directories once the pipeline is complete.
     *
     * Receives staged symlinks to the intermediate VCFs. Because Nextflow
     * stages files as symlinks pointing into the work directory, we resolve
     * each symlink to find the actual work directory and delete it.
     *
     * Targets:
     *   - GATHER_VCFS work dir       (cohort.genotyped.vcf.gz       ~2-5 GB)
     *   - GENOTYPE_FILTER work dir   (cohort.genotype_filtered.vcf.gz)
     *   - SPLIT_VARIANTS work dir    (snps.vcf.gz, indels.vcf.gz)
     *   - variant filter work dirs   (snps/indels.*.filtered.vcf.gz)
     *
     * Runs locally — no LSF submission needed.
     */

    tag "cleanup"
    label 'process_low'
    executor 'local'

    input:
    path multiqc_report              // sentinel — only runs after MultiQC completes
    tuple path(raw_vcf), path(raw_tbi)
    tuple path(geno_vcf), path(geno_tbi)
    tuple path(snp_vcf), path(snp_tbi)
    tuple path(indel_vcf), path(indel_tbi)
    tuple path(filt_snp_vcf), path(filt_snp_tbi)
    tuple path(filt_indel_vcf), path(filt_indel_tbi)

    output:
    val true, emit: done

    script:
    """
    echo "Pipeline complete — cleaning intermediate VCF work directories"

    # Resolve each symlink to its work directory and delete it
    for f in ${raw_vcf} ${geno_vcf} ${snp_vcf} ${indel_vcf} ${filt_snp_vcf} ${filt_indel_vcf}; do
        [ -z "\$f" ] && continue
        if [ -L "\$f" ]; then
            workdir=\$(dirname \$(readlink -f "\$f"))
            if [ -d "\$workdir" ] && echo "\$workdir" | grep -qE "/(work|scratch)/"; then
                echo "Removing work dir: \$workdir"
                rm -rf "\$workdir"
            fi
        fi
    done

    echo "Cleanup complete"
    """
}
