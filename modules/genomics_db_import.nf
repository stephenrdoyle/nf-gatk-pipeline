process GENOMICS_DB_IMPORT {
    /*
     * Consolidate all per-sample GVCFs for a single interval group
     * into a GenomicsDB workspace.
     *
     * Receives:
     *   group_id      - e.g. "group_0003"
     *   gvcf_list     - list of all .g.vcf.gz paths for this group
     *   gvcfi_list    - list of all .g.vcf.gz.tbi paths for this group
     *   interval_list - the primary (unpadded) interval .list file
     *
     * The workspace is named genomicsdb_${group_id} and is emitted
     * as a directory path for downstream GenotypeGVCFs.
     *
     * Note: uses primary interval (no padding) so workspace boundaries
     * are clean and non-overlapping across groups.
     */

    tag { group_id }
    label 'process_high'

    input:
    tuple val(group_id), path(gvcf_list), path(gvcfi_list), path(interval_list)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(group_id), path("genomicsdb_${group_id}"), path(interval_list), emit: genomicsdb

    script:
    def avail_mem = task.memory.toGiga() - 8
    def sample_map = gvcf_list.collect { gvcf ->
        // derive sample name from filename: sample_id.group_NNNN.g.vcf.gz
        def sname = gvcf.name.replaceAll(/\.group_\d+\.g\.vcf\.gz$/, '')
        "${sname}\t${gvcf}"
    }.join('\n')

    """
    # Write the sample map
    cat > sample_map.txt << 'SAMPLEMAP'
${sample_map}
SAMPLEMAP

    gatk GenomicsDBImport \\
        --java-options "-Xmx${avail_mem}g -Xms4g -Dsamjdk.use_intel_deflater=false" \\
        --sample-name-map sample_map.txt \\
        --genomicsdb-workspace-path genomicsdb_${group_id} \\
        --intervals ${interval_list} \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --reader-threads 8 \\
        --batch-size 50 \\
        --tmp-dir .
    """

    stub:
    """
    mkdir -p genomicsdb_${group_id}
    """
}
