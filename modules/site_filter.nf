process SITE_FILTER {
    /*
     * Site filters, in one pass over the genotype-filtered cohort VCF (this
     * replaces SPLIT_VARIANTS -> SNP/indel VariantFiltration -> MERGE_FILTERED
     * and FILTER_SUMMARY).
     *
     *   1. bcftools: one line per record with the annotations the rules need,
     *      plus F_MISSING (fraction of samples with no call, after the
     *      genotype filter)
     *   2. bin/site_filter.py apply: decides each record's FILTER (SNP rules
     *      for SNPs, indel rules for indels and mixed sites; a missing
     *      annotation passes, as in GATK), and writes the filter summary and
     *      the diagnostics report
     *   3. bcftools annotate copies the decisions into the VCF
     *   4. bin/site_filter.py verify: the FILTER values actually written must
     *      match the decisions, record for record by count, or the task fails
     *
     * Rules: --filter_mode hard (GATK's thresholds) or custom (--snp_filters,
     * --indel_filters); optional --max_depth, --max_missing, --excess_het.
     * Sites are flagged, never removed.
     */
    tag "site_filter"
    label 'process_medium'

    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'cohort.filtered.vcf.gz*'
    publishDir "${params.outdir}/filter_reports", mode: 'copy', pattern: '{filter_summary.tsv,filter_diagnostics.*}'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi"), emit: filtered_vcf
    path "filter_summary.tsv",                                                emit: summary
    path "filter_diagnostics.{html,tsv}",                                     emit: diagnostics

    script:
    // INFO tags in the site table: the hard-filter annotations, DP and
    // ExcessHet (diagnostics, --max_depth, --excess_het), and any tag a
    // custom expression uses
    def tags = ['QD', 'FS', 'SOR', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'DP', 'ExcessHet']
    if (params.filter_mode == 'custom') {
        "${params.snp_filters},${params.indel_filters}".tokenize(',').each { expr ->
            def m = expr =~ /^\s*(?:INFO\/)?([A-Za-z_][A-Za-z0-9_.]*)/
            if (m.find() && !(m.group(1) in tags)) tags << m.group(1)
        }
    }
    def query_fmt = (['%CHROM', '%POS', '%REF', '%ALT', '%INFO/F_MISSING'] + tags.collect { "%INFO/${it}" }).join('\\t') + '\\n'
    def opt = { flag, value -> value != null && value != false && value != '' ? "${flag} '${value}'" : '' }
    def rule_args = params.filter_mode == 'custom'
        ? "--mode custom ${opt('--snp-filters', params.snp_filters)} ${opt('--indel-filters', params.indel_filters)}"
        : '--mode hard'
    def excess_het = params.excess_het == true ? 'true' : params.excess_het
    """
    set -o pipefail

    # ── 1. site table (F_MISSING recomputed from the genotype-filtered calls) ──
    bcftools view --header-only ${vcf} | grep '^##contig' > contigs.txt || true
    bcftools +fill-tags --threads ${task.cpus} ${vcf} -Ou -- -t F_MISSING \\
        | bcftools query -f '${query_fmt}' \\
        | gzip -1 > sites.tsv.gz

    # ── 2. decisions (sites-only VCF) + summary + diagnostics ──────────────────
    site_filter.py apply \\
        --sites sites.tsv.gz \\
        --columns ${tags.join(',')} \\
        --contigs contigs.txt \\
        ${rule_args} \\
        ${opt('--max-depth', params.max_depth)} \\
        ${opt('--max-missing', params.max_missing)} \\
        ${opt('--excess-het', excess_het)} \\
        | bcftools view -Oz -o decisions.vcf.gz
    bcftools index --tbi decisions.vcf.gz
    rm sites.tsv.gz

    # ── 3. copy FILTER into the cohort VCF (matched on CHROM, POS, REF, ALT) ──
    bcftools annotate --threads ${task.cpus} \\
        --annotations decisions.vcf.gz \\
        --columns FILTER \\
        -Oz -o cohort.filtered.vcf.gz \\
        ${vcf}
    bcftools index --tbi --threads ${task.cpus} cohort.filtered.vcf.gz

    # ── 4. every record must carry the FILTER that was decided for it ──────────
    bcftools query -f '%FILTER\\n' cohort.filtered.vcf.gz > observed_filters.txt
    site_filter.py verify expected_filters.tsv observed_filters.txt
    rm observed_filters.txt
    """

    stub:
    """
    touch cohort.filtered.vcf.gz cohort.filtered.vcf.gz.tbi filter_summary.tsv filter_diagnostics.html filter_diagnostics.tsv
    """
}
