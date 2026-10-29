process VALIDATE_BAM {
    /*
     * Make sure each BAM has usable read groups for GATK: one sample name,
     * and PL/LB/PU on every @RG line.
     *
     * Most BAMs already do, so they are passed through untouched as
     * symlinks — no copy, no re-index. Only BAMs whose header needs a
     * change are rewritten:
     *   - no @RG lines      -> addreplacerg (reads need an RG tag; full rewrite)
     *   - header needs edit -> samtools reheader (copies the BAM) + index
     *
     * SM handling: a BAM with exactly one sample name on all its @RG lines
     * keeps it (GenomicsDBImport's sample map names it ${sample_id} anyway).
     * Otherwise — several sample names, or some @RG lines without SM — SM is
     * set to ${sample_id} on every @RG line, so GATK sees a single sample.
     * Existing read groups (IDs, and RG tags on reads) are kept.
     */
    tag "${sample_id}"
    label 'process_low'

    conda 'bioconda::samtools=1.19'
    container 'biocontainers/samtools:1.19--h50ea8bc_1'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(sample_id), path("${sample_id}.validated.bam"), path("${sample_id}.validated.bam.bai"), emit: validated_bam

    script:
    """
    # --no-PG: without it samtools adds its own @PG line, and the header
    # would never compare equal to the patched one below
    samtools view -H --no-PG ${bam} > header.txt

    # ── 1. No @RG lines: reads carry no RG tag, so they must be rewritten ───
    if ! grep -q '^@RG' header.txt; then
        echo "WARN: BAM ${bam} has no @RG lines -- adding read group for ${sample_id}"
        samtools addreplacerg \\
            -r "ID:1\tSM:${sample_id}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \\
            -o ${sample_id}.validated.bam \\
            ${bam}
        samtools index -@ ${task.cpus} ${sample_id}.validated.bam
        exit 0
    fi

    # ── 2. Decide whether SM must be forced to ${sample_id} ─────────────────
    # 0 = exactly one sample name, present on every @RG line (keep it)
    # 1 = several sample names, or some @RG lines without SM (force it)
    force_sm=\$(awk -F'\\t' '
        /^@RG/ {
            s = ""
            for (i = 2; i <= NF; i++) if (\$i ~ /^SM:/) s = \$i
            if (s == "") missing++; else names[s] = 1
        }
        END { n = 0; for (k in names) n++; print (n == 1 && missing == 0) ? 0 : 1 }
    ' header.txt)

    # ── 3. Build the patched header: set SM if forced, add missing PL/LB/PU ─
    awk -F'\\t' -v OFS='\\t' -v id='${sample_id}' -v force="\${force_sm}" '
        /^@RG/ {
            line = \$1; has_sm = has_pl = has_lb = has_pu = 0
            for (i = 2; i <= NF; i++) {
                if (\$i ~ /^SM:/) { if (force) continue; has_sm = 1 }
                if (\$i ~ /^PL:/) has_pl = 1
                if (\$i ~ /^LB:/) has_lb = 1
                if (\$i ~ /^PU:/) has_pu = 1
                line = line OFS \$i
            }
            if (!has_sm) line = line OFS "SM:" id
            if (!has_pl) line = line OFS "PL:ILLUMINA"
            if (!has_lb) line = line OFS "LB:lib1"
            if (!has_pu) line = line OFS "PU:unit1"
            print line; next
        }
        { print }
    ' header.txt > patched_header.txt

    # ── 4. Pass through untouched, or rewrite the header ────────────────────
    if cmp -s header.txt patched_header.txt; then
        # Header already valid: link to the original BAM and index (no copy)
        ln -s "\$(readlink -f ${bam})" ${sample_id}.validated.bam
        ln -s "\$(readlink -f ${bai})" ${sample_id}.validated.bam.bai
    else
        echo "WARN: BAM ${bam} read-group header patched (SM forced to ${sample_id}: \${force_sm}):"
        diff header.txt patched_header.txt || true
        samtools reheader patched_header.txt ${bam} > ${sample_id}.validated.bam
        samtools index -@ ${task.cpus} ${sample_id}.validated.bam
    fi
    """
}
