process GATHER_VCFS {
    /*
     * Concatenate all per-interval VCFs into a single cohort VCF,
     * then sort to enforce canonical contig order.
     *
     * Uses a sed-based header replacement rather than bcftools reheader
     * to avoid bcftools version-dependent merge vs replace behaviour.
     */

    label 'process_medium'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path vcf_list    // all group_NNNN.genotyped.vcf.gz, name-sorted
    path vcfi_list   // corresponding .tbi files
    path fai         // reference .fai — must match the reference used for calling

    output:
    path "cohort.genotyped.vcf.gz",     emit: vcf
    path "cohort.genotyped.vcf.gz.tbi", emit: tbi

    script:
    def sorted_vcfs = vcf_list instanceof List
        ? vcf_list.sort { it.name }.collect { it.toString() }.join(' ')
        : vcf_list.toString()
    """
    # 1. Concatenate in group order
    bcftools concat \\
        --naive \\
        -O z \\
        -o concat.vcf.gz \\
        ${sorted_vcfs}
    bcftools index --tbi concat.vcf.gz

    # 2. Build replacement header:
    #    - All ## meta lines except ##contig
    #    - Fresh ##contig lines from the .fai (correct lengths + order)
    #    - #CHROM column header line last
    bcftools view -h concat.vcf.gz | grep "^##" | grep -v "^##contig" > new_header.txt
    awk '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' ${fai}  >> new_header.txt
    bcftools view -h concat.vcf.gz | grep "^#CHROM"                  >> new_header.txt

    # 3. Replace header by stripping the old one and prepending the new one.
    #    This avoids bcftools reheader version-dependent merge behaviour.
    #    - Decompress, strip all header lines, reattach new header, recompress
    bgzip -d -c concat.vcf.gz | grep -v "^#" > body.vcf
    cat new_header.txt body.vcf | bgzip -c > reheadered.vcf.gz
    bcftools index --tbi reheadered.vcf.gz

    # 4. Sort to canonical contig order
    bcftools sort \\
        -m ${task.memory.toGiga() - 4}G \\
        -O z \\
        -o cohort.genotyped.vcf.gz \\
        reheadered.vcf.gz
    bcftools index --tbi cohort.genotyped.vcf.gz

    # Cleanup intermediates
    rm -f concat.vcf.gz concat.vcf.gz.tbi \\
          reheadered.vcf.gz reheadered.vcf.gz.tbi \\
          new_header.txt body.vcf
    """
}
