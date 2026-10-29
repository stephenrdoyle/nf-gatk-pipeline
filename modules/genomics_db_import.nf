process GENOMICS_DB_IMPORT {
    /*
     * Import all samples' GVCFs for one joint-genotyping piece into a
     * GenomicsDB workspace.
     *
     * Receives:
     *   group_id   - HaplotypeCaller window, e.g. "group_0003"
     *   gvcf_list  - that window's GVCFs, one per sample (name-sorted)
     *   gvcfi_list - their .tbi indexes
     *   part_id    - the piece, e.g. "group_0003.part_002"
     *   interval   - the piece as "contig:start-end": a --genotype_window_size
     *                slice of the window. The GVCF indexes let GenomicsDB read
     *                only this slice of each GVCF.
     *
     * Pieces don't overlap, so workspace boundaries are clean.
     *
     * --consolidate (--genomicsdb_consolidate, default true): samples are
     * imported in batches of 50, each stored as a separate fragment that
     * GenotypeGVCFs would otherwise have to merge on every read. Consolidating
     * merges them once at the end of the import.
     */

    tag { part_id }
    label 'process_high'

    input:
    tuple val(group_id), path(gvcf_list), path(gvcfi_list), val(part_id), val(interval)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(part_id), path("genomicsdb_${part_id}"), val(interval), emit: genomicsdb

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
        --genomicsdb-workspace-path genomicsdb_${part_id} \\
        --intervals ${interval} \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --reader-threads ${task.cpus} \\
        --batch-size 50 \\
        --consolidate ${params.genomicsdb_consolidate} \\
        --tmp-dir .
    """

    stub:
    """
    mkdir -p genomicsdb_${part_id}
    """
}
