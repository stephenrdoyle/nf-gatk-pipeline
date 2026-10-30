process PREFLIGHT {
    /*
     * Checks that need the farm's tools, run once before any BAM is processed
     * (main.nf holds VALIDATE_BAM, and so HaplotypeCaller, until this passes):
     *   - the process.beforeScript modules provide gatk, samtools, bcftools
     *     (with its setGT plugin) and python3
     *   - every BAM's @SQ contigs (names, lengths, order) match the
     *     reference .fai; GATK refuses mismatched dictionaries, and would
     *     otherwise fail once per window for each bad BAM
     * All problems are reported together. The cheaper checks (samplesheet,
     * missing or stale indexes, writable reference dir) run on the head node
     * at launch; see preflight_launch_checks() in main.nf.
     *
     * cache false: BAMs or modules can change between runs without changing
     * this task's inputs, so it re-checks every time (a few minutes).
     */
    tag "preflight"
    label 'process_low'
    cache false

    input:
    val  bam_list   // "sample_id<TAB>/abs/path.bam" lines, one per sample
    path fai

    output:
    val true, emit: ok

    script:
    """
    errors=0
    fail() { echo "PREFLIGHT ERROR: \$*" >&2; errors=\$((errors + 1)); }

    # ── tools from process.beforeScript ───────────────────────────────────────
    for tool in gatk samtools bcftools python3; do
        command -v "\$tool" >/dev/null 2>&1 \\
            || fail "'\$tool' not found: check the module load lines in process.beforeScript (nextflow.config)"
    done
    # GENOTYPE_FILTER uses the bcftools setGT plugin
    if command -v bcftools >/dev/null 2>&1; then
        bcftools plugin -l 2>/dev/null | grep -qw setGT \\
            || fail "bcftools plugin 'setGT' not found (needed by GENOTYPE_FILTER): check BCFTOOLS_PLUGINS / the bcftools module"
    fi

    # ── BAM contigs vs reference ──────────────────────────────────────────────
    cut -f1,2 ${fai} > reference_contigs.tsv
    cat > bams.tsv << 'BAMS'
${bam_list}
BAMS

    if command -v samtools >/dev/null 2>&1; then
        while IFS=\$'\\t' read -r sample bam; do
            [ -n "\$sample" ] || continue
            if ! samtools view -H "\$bam" > header.sam 2> header.err; then
                fail "\$sample: cannot read the BAM header of \$bam (\$(head -1 header.err))"
                continue
            fi
            awk -F'\\t' '\$1 == "@SQ" {
                    name = ""; len = ""
                    for (i = 2; i <= NF; i++) {
                        if (\$i ~ /^SN:/) name = substr(\$i, 4)
                        if (\$i ~ /^LN:/) len  = substr(\$i, 4)
                    }
                    print name "\\t" len
                }' header.sam > bam_contigs.tsv
            if ! cmp -s reference_contigs.tsv bam_contigs.tsv; then
                if cmp -s <(sort reference_contigs.tsv) <(sort bam_contigs.tsv); then
                    detail="same contigs and lengths, but in a different order"
                else
                    detail="\$(grep -c . bam_contigs.tsv) BAM vs \$(grep -c . reference_contigs.tsv) reference contigs; first difference, < reference > BAM: \$(diff reference_contigs.tsv bam_contigs.tsv | grep '^[<>]' | head -2 | tr '\\t\\n' ': ')"
                fi
                fail "\$sample: BAM contigs don't match the reference .fai (\$detail): \$bam"
            fi
        done < bams.tsv
    fi

    if [ "\$errors" -gt 0 ]; then
        echo "PREFLIGHT FAILED: \$errors problem(s), listed above. Nothing has been submitted for the BAMs." >&2
        exit 1
    fi
    echo "Preflight passed: tools found, \$(grep -c . bams.tsv) BAMs match the reference"
    """
}
