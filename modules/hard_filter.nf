process HARD_FILTER {
    tag "${variant_type}"
    label 'process_low'

    conda 'bioconda::gatk4=4.6.1.0'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(variant_type), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(variant_type), path("${variant_type}.hard_filtered.vcf.gz"), path("${variant_type}.hard_filtered.vcf.gz.tbi"), emit: filtered_vcf
    path "${variant_type}_filter_summary.tsv", emit: summary

    script:
    // Same reference form as HAPLOTYPE_CALLER / GENOTYPE_GVCFS: the real path,
    // with symlinks (e.g. lustre_link) resolved. GATK runs in a Singularity
    // wrapper, so a resolved /lustre path doesn't depend on the container
    // also seeing a staged symlink's target; .fai/.dict are published next
    // to the real FASTA by SAMTOOLS_FAIDX / CREATE_DICT
    def ref_abs = reference.toRealPath()
    // Heap tied to the LSF request: without -Xmx the JVM sizes its heap from
    // the node's RAM, not task.memory, and LSF kills the job (exit 130)
    def java_opts = "-Xmx${task.memory.toGiga() - 4}g -Dsamjdk.use_intel_deflater=false"
    // [ annotation, operator, threshold, filter name ] — single source for both
    // the VariantFiltration arguments and the threshold TSV
    def filters = variant_type == 'snps' ? [
        [ 'QD',             '<', '2.0',   'SNP_QD2'                ],
        [ 'MQ',             '<', '40.0',  'SNP_MQ40'               ],
        [ 'FS',             '>', '60.0',  'SNP_FS60'               ],
        [ 'SOR',            '>', '3.0',   'SNP_SOR3'               ],
        [ 'MQRankSum',      '<', '-12.5', 'SNP_MQRankSum-12.5'     ],
        [ 'ReadPosRankSum', '<', '-8.0',  'SNP_ReadPosRankSum-8'   ],
    ] : [
        [ 'QD',             '<', '2.0',   'INDEL_QD2'              ],
        [ 'FS',             '>', '200.0', 'INDEL_FS200'            ],
        [ 'SOR',            '>', '10.0',  'INDEL_SOR10'            ],
        [ 'ReadPosRankSum', '<', '-20.0', 'INDEL_ReadPosRankSum-20'],
    ]
    def filter_args = filters.collect { ann, op, thr, name ->
        "--filter-expression \"${ann} ${op} ${thr}\" --filter-name \"${name}\""
    }.join(' \\\n        ')
    def summary_rows = filters.collect { ann, op, thr, name ->
        [ variant_type, ann, op, thr, 'NA', 'NA', 'NA', 'NA', 'NA' ].join('\t')
    }.join('\n')
    """
    gatk VariantFiltration \\
        --java-options "${java_opts}" \\
        --reference ${ref_abs} \\
        --variant ${vcf} \\
        ${filter_args} \\
        --output ${variant_type}.hard_filtered.vcf.gz

    # Write threshold TSV for downstream summary
    cat > ${variant_type}_filter_summary.tsv << 'TSV'
variant_type\tannotation\toperator\tthreshold\tn_total\tn_flagged\tpct_flagged\tmean\tmedian
${summary_rows}
TSV
    """
}
