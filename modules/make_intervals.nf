process MAKE_INTERVALS {
    /*
     * Generate fixed-width interval list files from the reference .fai.
     * Produces one .list file per window, named group_NNNN.list.
     * Chromosome boundaries are respected — windows never span contigs.
     * Mitochondrion (any contig <= window_size) gets emitted as one interval.
     *
     * Output: a flat list of Path objects, one per interval group.
     */

    label 'process_low'

    input:
    path fai
    val  window_size

    output:
    path "intervals/group_*.list", emit: interval_lists

    script:
    """
    mkdir -p intervals

    # Build BED windows from .fai, respecting contig boundaries.
    # Contigs shorter than window_size are emitted as a single interval.
    awk -v W=${window_size} '
    {
        contig = \$1
        len    = \$2
        if (len <= W) {
            print contig "\t" 0 "\t" len
        } else {
            start = 0
            while (start < len) {
                end = start + W
                if (end > len) end = len
                print contig "\t" start "\t" end
                start = end
            }
        }
    }' ${fai} > windows.bed

    # Convert BED (0-based half-open) to GATK interval format (1-based closed).
    # Write one group_NNNN.list file per window.
    awk '{
        outfile = sprintf("intervals/group_%04d.list", NR-1)
        printf "%s:%d-%d\\n", \$1, \$2+1, \$3 > outfile
        close(outfile)
    }' windows.bed

    echo "Generated \$(ls intervals/group_*.list | wc -l) interval files"
    """
}
