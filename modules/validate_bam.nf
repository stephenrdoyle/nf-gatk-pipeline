process VALIDATE_BAM {
    tag "${sample_id}"
    label 'process_low'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(sample_id), path("${sample_id}.validated.bam"), path("${sample_id}.validated.bam.bai"), emit: validated_bam

    script:
    """
    # ── 1. Check the number of unique SM tags ─────────────────────────────────
    n_samples=\$(samtools view -H ${bam} | grep '^@RG' | grep -o 'SM:[^\t]*' | sort -u | wc -l)

    if [ "\${n_samples}" -eq 0 ]; then
        # No @RG lines at all — add a complete read group from scratch
        echo "WARN: BAM ${bam} has no @RG lines -- adding read group for ${sample_id}"
        samtools addreplacerg \
            -r "ID:1\tSM:${sample_id}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
            -o normalised.bam \
            ${bam}
        input_bam=normalised.bam

    elif [ "\${n_samples}" -gt 1 ]; then
        # Multiple SM tags — normalise to single sample name
        echo "WARN: BAM ${bam} contains \${n_samples} sample names -- normalising to ${sample_id} via addreplacerg"
        samtools addreplacerg \
            -r "ID:1\tSM:${sample_id}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
            -w \
            -o normalised.bam \
            ${bam}
        input_bam=normalised.bam

    else
        # Exactly one SM tag — use as-is
        input_bam=${bam}
    fi

    # ── 2. Patch any @RG lines missing PL/SM/LB/PU tags ─────────────────────
    samtools view -H \${input_bam} | awk '
        /^@RG/ {
            if (index(\$0, "\\tPL:") == 0) \$0 = \$0 "\\tPL:ILLUMINA"
            if (index(\$0, "\\tSM:") == 0) \$0 = \$0 "\\tSM:${sample_id}"
            if (index(\$0, "\\tLB:") == 0) \$0 = \$0 "\\tLB:lib1"
            if (index(\$0, "\\tPU:") == 0) \$0 = \$0 "\\tPU:unit1"
        }
        { print }
    ' > patched_header.txt

    samtools reheader patched_header.txt \${input_bam} > ${sample_id}.validated.bam

    samtools index -@ ${task.cpus} ${sample_id}.validated.bam
    """
}
