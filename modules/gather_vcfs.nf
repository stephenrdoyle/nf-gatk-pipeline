process GATHER_VCFS {
    /*
     * Concatenate all per-interval VCFs into a single cohort VCF and
     * replace its ##contig header lines with ones built from the .fai.
     *
     * No decompression and no sort:
     *   - groups are numbered in .fai order and never overlap, and
     *     GenotypeGVCFs emits each variant only from the group it starts
     *     in, so concatenating in group-name order is already sorted
     *   - concat --naive copies compressed blocks as-is
     *   - reheader -h replaces the whole header (no merging) and copies the
     *     record blocks unchanged
     * Indexing the output fails on any out-of-order record, so the final
     * `bcftools index` also checks that the concatenation is sorted.
     */

    label 'process_medium'
    publishDir "${params.outdir}/variants", mode: params.publish_mode   // see publish_mode in nextflow.config

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
    # 1. Concatenate in group order (compressed blocks copied, not recompressed)
    bcftools concat \\
        --naive \\
        -O z \\
        -o concat.vcf.gz \\
        ${sorted_vcfs}

    # 2. Build replacement header:
    #    - All ## meta lines except ##contig
    #    - Fresh ##contig lines from the .fai (correct lengths + order)
    #    - #CHROM column header line last
    bcftools view -h concat.vcf.gz | grep "^##" | grep -v "^##contig" > new_header.txt
    awk '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' ${fai}  >> new_header.txt
    bcftools view -h concat.vcf.gz | grep "^#CHROM"                  >> new_header.txt

    # 3. Swap in the new header; the record blocks are copied unchanged
    bcftools reheader \\
        -h new_header.txt \\
        -o cohort.genotyped.vcf.gz \\
        concat.vcf.gz

    # 4. Index — errors on unsorted positions, so this also verifies order
    bcftools index --tbi cohort.genotyped.vcf.gz

    # Cleanup intermediates
    rm -f concat.vcf.gz new_header.txt
    """
}
