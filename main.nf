/*
 * nf-gatk-pipeline: main workflow
 * WGS variant calling pipeline
 * HaplotypeCaller -> GenomicsDBImport -> GenotypeGVCFs -> Filter -> QC
 *
 * Interval-scatter design:
 *   - Fixed 10 Mb windows generated from .fai at runtime
 *   - HaplotypeCaller runs per-sample x per-interval (with 200 bp padding)
 *   - GenomicsDBImport runs per-interval (all samples, primary interval)
 *   - GenotypeGVCFs runs per-interval
 *   - GatherVcfs concatenates all intervals into final cohort VCF
 *   - Post-calling: genotype filter -> site filter -> QC
 */

nextflow.enable.dsl = 2

include { VALIDATE_BAM        } from './modules/validate_bam'
include { SAMTOOLS_FAIDX      } from './modules/samtools_faidx'
include { CREATE_DICT         } from './modules/create_dict'
include { MAKE_INTERVALS      } from './modules/make_intervals'
include { HAPLOTYPE_CALLER    } from './modules/haplotype_caller'
include { GENOMICS_DB_IMPORT  } from './modules/genomics_db_import'
include { GENOTYPE_GVCFS      } from './modules/genotype_gvcfs'
include { GNARLY_GENOTYPER    } from './modules/gnarly_genotyper'
include { GATHER_VCFS         } from './modules/gather_vcfs'
include { GENOTYPE_FILTER     } from './modules/genotype_filter'
include { SITE_FILTER         } from './modules/site_filter'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_RAW      } from './modules/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_FILTERED } from './modules/bcftools_stats'
include { SAMPLE_MISSINGNESS as SAMPLE_MISSINGNESS_RAW      } from './modules/sample_missingness'
include { SAMPLE_MISSINGNESS as SAMPLE_MISSINGNESS_FILTERED } from './modules/sample_missingness'
include { MULTIQC             } from './modules/multiqc'
include { CLEANUP             } from './modules/cleanup'
include { POOL_ANNOTATE       } from './modules/pool_annotate'
include { ESTIMATE_HETEROZYGOSITY; HETEROZYGOSITY_SUMMARY } from './modules/heterozygosity'
include { PREFLIGHT           } from './modules/preflight'
include { SOFTWARE_VERSIONS   } from './modules/software_versions'

// ---------------------------------------------------------------------------
// Parameter defaults: defined once, in the params block of nextflow.config
// ---------------------------------------------------------------------------

// Guard: required inputs
if (!params.samplesheet || !params.reference) {
    error "ERROR: --samplesheet and --reference are required"
}

// Guard: joint-genotyping piece size
if (!(params.genotype_window_size instanceof Number) || params.genotype_window_size < 1) {
    error "ERROR: --genotype_window_size must be a positive number of bp (got ${params.genotype_window_size})"
}

// Guard: --heterozygosity is unset, 'auto', or a fraction
if (params.heterozygosity != null && params.heterozygosity != 'auto' &&
        !(params.heterozygosity instanceof Number && params.heterozygosity > 0 && params.heterozygosity < 1)) {
    error "ERROR: --heterozygosity must be 'auto' or a number between 0 and 1 (got ${params.heterozygosity})"
}

// Guard: site filters (see modules/site_filter.nf)
if (params.filter_mode == 'optimised') {
    error "ERROR: --filter_mode optimised has been removed: its percentile thresholds flagged a fixed share of " +
          "sites whatever the data, and on the F2 cohort were no better than GATK's (Ts/Tv 2.00 vs 1.98) while " +
          "losing 2.3 M more SNPs. Use the default (hard), and check filter_reports/filter_diagnostics.html; " +
          "tune with --filter_mode custom if a threshold doesn't suit the cohort."
}
if (!(params.filter_mode in ['hard', 'custom'])) {
    error "ERROR: --filter_mode must be 'hard' or 'custom' (got '${params.filter_mode}')"
}
if (params.filter_mode == 'custom') {
    if (!params.snp_filters || !params.indel_filters) {
        error "ERROR: --filter_mode custom requires --snp_filters and --indel_filters"
    }
    // same grammar as bin/site_filter.py: '<INFO tag> <op> <number>', comma-separated
    def bad = "${params.snp_filters},${params.indel_filters}".tokenize(',')*.trim().findAll { e ->
        e && !(e ==~ /(?:INFO\/)?[A-Za-z_][A-Za-z0-9_.]*\s*(<=|>=|==|!=|<|>)\s*-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?/)
    }
    if (bad) error "ERROR: can't parse filter expression(s) ${bad.join('; ')}: expected '<INFO tag> <op> <number>', e.g. 'QD < 2.0'"
}
if (params.max_depth != null && params.max_depth != false &&
        !(params.max_depth.toString() ==~ /\d+(\.\d+)?x?/)) {
    error "ERROR: --max_depth must be a multiple of the median site depth (e.g. 2x) or a number (got ${params.max_depth})"
}
if (params.max_missing != null && params.max_missing != false &&
        !(params.max_missing instanceof Number && params.max_missing >= 0 && params.max_missing < 1)) {
    error "ERROR: --max_missing must be a fraction from 0 to 1 (got ${params.max_missing})"
}
if (!(params.excess_het in [true, false, null] || (params.excess_het instanceof Number && params.excess_het > 0))) {
    error "ERROR: --excess_het must be true (GATK's 54.69), a number, or false (got ${params.excess_het})"
}

// ---------------------------------------------------------------------------
// Launch preflight: instant checks on the launch host, before anything is
// submitted. Checks that need the farm's tools (modules, BAM contigs vs the
// reference) run in the PREFLIGHT process. All problems are reported together.
// ---------------------------------------------------------------------------
def preflight_launch_checks() {
    def problems = []
    def sheet = file(params.samplesheet)
    def rows  = sheet.splitCsv(header: true)
    def cols  = rows ? rows[0].keySet() : []
    ['sample_id', 'bam'].each { c ->
        if (!(c in cols)) problems << "samplesheet ${sheet} has no '${c}' column (header: ${cols.join(',') ?: 'empty file'})"
    }

    if (!problems) {
        def ids = rows.collect { (it.sample_id ?: '').trim() }
        if (ids.any { !it }) problems << "samplesheet has rows with an empty sample_id"
        // sample_id becomes part of file names (<sample_id>.group_NNNN.g.vcf.gz)
        def unsafe = ids.findAll { it && !(it ==~ /[A-Za-z0-9._-]+/) }
        if (unsafe) problems << "sample_ids may only use letters, digits, '.', '_' and '-': ${unsafe.take(5).join(', ')}"
        // duplicates would overwrite each other's GVCFs, silently losing a sample
        def dups = ids.countBy { it }.findAll { id, n -> id && n > 1 }.keySet()
        if (dups) problems << "duplicate sample_ids: ${dups.take(5).join(', ')}"

        rows.each { row ->
            def bam = file(row.bam)
            if (!bam.exists()) { problems << "${row.sample_id}: BAM not found: ${row.bam}"; return }
            // same index lookup as parse_samplesheet()
            def bai = row.bai ? file(row.bai)
                    : file("${row.bam}.bai").exists() ? file("${row.bam}.bai")
                    : file(row.bam.replaceAll(/\.bam$/, '.bai'))
            if (!bai.exists()) problems << "${row.sample_id}: BAM index not found (tried ${row.bai ?: row.bam + '.bai and .bai'})"
            else if (bai.lastModified() < bam.lastModified()) problems << "${row.sample_id}: index is older than the BAM, so it may be stale: ${bai}"
        }
    }

    // .fai/.dict are published beside the reference (GATK reads them there);
    // that needs write access unless both already exist
    def ref = file(params.reference)
    if (!ref.exists()) problems << "reference not found: ${params.reference}"
    else {
        def have_index = file("${ref}.fai").exists() && file("${ref.parent}/${ref.baseName}.dict").exists()
        if (!have_index && !java.nio.file.Files.isWritable(ref.parent))
            problems << "reference directory ${ref.parent} is not writable, and .fai/.dict must be created there"
    }

    if (problems) error "Preflight failed (${problems.size()} problem(s)):\n  - " + problems.join('\n  - ')
}
preflight_launch_checks()

// ---------------------------------------------------------------------------
// Helper: parse samplesheet CSV -> [ sample_id, bam_path, bai_path ]
// ---------------------------------------------------------------------------
def parse_samplesheet(csv_path) {
    Channel.fromPath(csv_path)
        .splitCsv(header: true)
        .map { row ->
            def bam = file(row.bam, checkIfExists: true)
            // Accept explicit bai column, or derive from bam path
            def bai = row.bai
                ? file(row.bai, checkIfExists: true)
                : file("${row.bam}.bai").exists()
                    ? file("${row.bam}.bai")
                    : file(row.bam.replaceAll(/\.bam$/, '.bai'), checkIfExists: true)
            [ row.sample_id, bam, bai ]
        }
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    // --- Software versions (results/pipeline_info/software_versions.txt) ---
    // Pipeline commit is read on the launch host; '-dirty' marks uncommitted
    // changes, i.e. a run whose code isn't exactly any commit
    def git_commit = 'unknown (not a git checkout)'
    try {
        def sha   = "git -C ${projectDir} rev-parse --short HEAD".execute().text.trim()
        def dirty = "git -C ${projectDir} status --porcelain --untracked-files=no".execute().text.trim()
        if (sha) git_commit = sha + (dirty ? '-dirty' : '')
    } catch (Exception e) { }
    SOFTWARE_VERSIONS([
        "pipeline: ${workflow.manifest.name} ${workflow.manifest.version}",
        "git commit: ${git_commit}",
        "nextflow: ${workflow.nextflow.version}",
        "run name: ${workflow.runName}",
        "command: ${workflow.commandLine}",
    ].join('\n'))

    // --- Publishing: result VCFs are hard-linked by default (publish_mode) ---
    if (!(params.publish_mode in ['link', 'copy'])) {
        error "ERROR: --publish_mode must be 'link' or 'copy' (got '${params.publish_mode}')"
    }
    if (params.publish_mode == 'link') {
        // Hard links can't cross filesystems: compare the filesystem of the
        // nearest existing directory of --outdir with that of the work dir
        // (null, so the check is skipped, if the filesystem can't be determined)
        def store_of = { java.nio.file.Path p ->
            try {
                p = p.toAbsolutePath()
                while (p != null && !java.nio.file.Files.exists(p)) p = p.parent
                return p != null ? java.nio.file.Files.getFileStore(p) : null
            } catch (Exception e) {
                return null
            }
        }
        def out_store  = store_of(file(params.outdir))
        def work_store = store_of(workflow.workDir)
        if (out_store && work_store && out_store != work_store) {
            error "ERROR: --outdir ${params.outdir} is on a different filesystem from the work directory " +
                  "(${workflow.workDir}), so result VCFs can't be hard-linked: use --publish_mode copy"
        }
    }

    // --- 0. Reference files ---
    reference = file(params.reference, checkIfExists: true)

    // .fai/.dict beside the FASTA are never overwritten (see SAMTOOLS_FAIDX),
    // so an index left over from an older FASTA would be used as-is: refuse
    // one that is older than the FASTA
    [ file("${reference}.fai"), file("${reference.parent}/${reference.baseName}.dict") ].each { idx ->
        if (idx.exists() && idx.lastModified() < reference.lastModified()) {
            error "ERROR: ${idx} is older than the reference FASTA, so it may describe a previous " +
                  "version of it. Delete it and relaunch; the pipeline will recreate it."
        }
    }

    SAMTOOLS_FAIDX(reference)
    reference_fai = SAMTOOLS_FAIDX.out.fai

    CREATE_DICT(reference)
    reference_dict = CREATE_DICT.out.dict

    // --- 1. Validate / normalise BAM read groups ---
    // PREFLIGHT checks tools and BAM contigs; each BAM is held until it
    // passes. The combine/map only delays the channel: VALIDATE_BAM still gets
    // exactly [ sample_id, bam, bai ], so cached tasks stay cached.
    PREFLIGHT(
        parse_samplesheet(params.samplesheet)
            .map { sample_id, bam, bai -> "${sample_id}\t${bam}" }
            .collect()
            .map { it.join('\n') },
        reference_fai
    )
    raw_bams = parse_samplesheet(params.samplesheet)
        .combine(PREFLIGHT.out.ok)
        .map { sample_id, bam, bai, ok -> [ sample_id, bam, bai ] }
    VALIDATE_BAM(raw_bams, reference, reference_fai, reference_dict)
    validated_bams = VALIDATE_BAM.out.validated_bam

    // --- Heterozygosity prior for GenotypeGVCFs (--heterozygosity) ---
    //   unset   GATK's defaults (0.001 SNP, 1.25e-4 indel)
    //   number  that value; indels = --indel_heterozygosity or SNP x 0.125;
    //           --heterozygosity-stdev = --heterozygosity_stdev or sqrt(het / 10)
    //   'auto'  estimated from the BAMs with ANGSD (modules/heterozygosity.nf),
    //           alongside HaplotypeCaller, from the first
    //           --heterozygosity_max_samples samples by sample_id
    // Not applied in pooled mode (ANGSD's per-sample model assumes diploid
    // individuals) or with --use_gnarly (GnarlyGenotyper takes no such prior).
    def use_prior = params.heterozygosity != null && !params.pooled && !params.use_gnarly
    if (params.heterozygosity != null && !use_prior) {
        log.warn "--heterozygosity is ignored with ${params.pooled ? '--pooled' : '--use_gnarly'}: GATK's default priors are used"
    }
    if (use_prior && params.heterozygosity == 'auto') {
        ESTIMATE_HETEROZYGOSITY(
            validated_bams
                .toSortedList { a, b -> a[0] <=> b[0] }
                .flatMap { it.take(params.heterozygosity_max_samples as int) },
            reference,
            reference_fai,
            params.heterozygosity_sample_bp
        )
        HETEROZYGOSITY_SUMMARY(
            ESTIMATE_HETEROZYGOSITY.out.sfs.collect(),
            params.indel_heterozygosity ?: '',
            params.heterozygosity_stdev ?: ''
        )
        gatk_prior_args = HETEROZYGOSITY_SUMMARY.out.args.map { it.text.trim() }
    } else if (use_prior) {
        def het   = params.heterozygosity as double
        def indel = params.indel_heterozygosity ? params.indel_heterozygosity as double : het * 0.125
        // stdev = sqrt(het / 10) keeps the prior's strength at GATK's default;
        // see modules/heterozygosity.nf
        def stdev = params.heterozygosity_stdev ? params.heterozygosity_stdev as double : Math.sqrt(het / 10)
        gatk_prior_args = Channel.value(("--heterozygosity ${het} --indel-heterozygosity ${String.format('%.3g', indel)} " +
                                         "--heterozygosity-stdev ${String.format('%.3g', stdev)}").toString())
    } else {
        gatk_prior_args = Channel.value('')
    }

    // --- 2. Generate interval lists from .fai ---
    MAKE_INTERVALS(reference_fai, params.window_size)
    interval_ch = MAKE_INTERVALS.out.interval_lists.flatten()

    // --- 3. HaplotypeCaller: per-sample x per-interval ---
    HAPLOTYPE_CALLER(
        validated_bams,
        interval_ch,
        reference,
        reference_fai,
        reference_dict
    )

    // --- 4. Reshape for GenomicsDBImport ---
    // groupKey carries the expected group size (one GVCF per sample), so each
    // group is released as soon as all its samples finish HaplotypeCaller,
    // rather than when every HaplotypeCaller job in the run has finished.
    // Sorting the grouped files by name makes the GVCF list (and so the
    // sample map and task hash) identical between runs, so -resume reuses
    // completed GenomicsDBImport tasks instead of re-running them.
    n_samples = raw_bams.count()

    gvcfs_by_group = HAPLOTYPE_CALLER.out.gvcfs
        .combine(n_samples)
        .map { sample_id, group_id, gvcf, tbi, n -> [ groupKey(group_id, n), gvcf, tbi ] }
        .groupTuple(by: 0, sort: { a, b -> a.name <=> b.name })
        .map { key, gvcfs, tbis -> [ key.getGroupTarget().toString(), gvcfs, tbis ] }

    // Joint-genotyping pieces: each HaplotypeCaller window cut into slices of
    // --genotype_window_size for GenomicsDBImport and GenotypeGVCFs, so those
    // run as many shorter jobs. Computed here, not in MAKE_INTERVALS, so the
    // HaplotypeCaller windows (and their cache) are unchanged.
    // Emits [ group_id, part_id, "contig:start-end" ]; part_ids sort in
    // genomic order (group_0000.part_000, group_0000.part_001, ...).
    genotype_pieces = interval_ch.flatMap { list_file ->
        def group_id = list_file.baseName
        def size     = params.genotype_window_size as long
        def pieces   = []
        list_file.readLines().findAll { it.trim() }.each { line ->
            def colon  = line.lastIndexOf(':')
            def contig = line.substring(0, colon)
            def range  = line.substring(colon + 1).tokenize('-')*.toLong()
            for (long s = range[0]; s <= range[1]; s += size) {
                pieces << "${contig}:${s}-${Math.min(s + size - 1, range[1])}".toString()
            }
        }
        pieces.withIndex().collect { iv, i ->
            [ group_id, "${group_id}.part_${String.format('%03d', i)}".toString(), iv ]
        }
    }

    // Each window's GVCFs go to every piece of that window
    genomicsdb_input = gvcfs_by_group
        .combine(genotype_pieces, by: 0)
        // [ group_id, gvcfs, tbis, part_id, interval ]

    // --- 5. GenomicsDBImport: per piece, all samples ---
    GENOMICS_DB_IMPORT(
        genomicsdb_input,
        reference,
        reference_fai,
        reference_dict
    )

    // --- 6. Genotype: GenotypeGVCFs or GnarlyGenotyper, per piece ---
    if (params.use_gnarly) {
        GNARLY_GENOTYPER(
            GENOMICS_DB_IMPORT.out.genomicsdb,
            reference,
            reference_fai,
            reference_dict
        )
        genotyped_vcfs = GNARLY_GENOTYPER.out.vcf
    } else {
        GENOTYPE_GVCFS(
            GENOMICS_DB_IMPORT.out.genomicsdb,
            gatk_prior_args,
            reference,
            reference_fai,
            reference_dict
        )
        genotyped_vcfs = GENOTYPE_GVCFS.out.vcf
    }

    // --- 7. Gather all piece VCFs, in part_id (= genomic) order, into the cohort VCF ---
    all_vcfs = genotyped_vcfs
        .map { part_id, vcf, tbi -> [ part_id, vcf, tbi ] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { sorted ->
            [ sorted.collect { it[1] }, sorted.collect { it[2] } ]
        }

    GATHER_VCFS(
        all_vcfs.map { vcfs, tbis -> vcfs },
        all_vcfs.map { vcfs, tbis -> tbis },
        reference_fai
    )
    // Pair the two single-item GATHER_VCFS outputs into one tuple channel
    cohort_vcf = GATHER_VCFS.out.vcf
        .combine(GATHER_VCFS.out.tbi)
    // cohort_vcf: [ cohort.genotyped.vcf.gz, cohort.genotyped.vcf.gz.tbi ]

    // --- 8. Stats + per-sample missingness on the raw cohort VCF (all sites) ---
    BCFTOOLS_STATS_RAW(cohort_vcf, 'before_filter', false)
    SAMPLE_MISSINGNESS_RAW(cohort_vcf, BCFTOOLS_STATS_RAW.out.stats, 'before_filter', false)

    // --- 9. Per-genotype filter (DP / GQ / RGQ) ---
    GENOTYPE_FILTER(cohort_vcf)
    geno_filtered_vcf = GENOTYPE_FILTER.out.filtered_vcf

    // --- 10. Site filters (hard or custom, plus optional depth, missingness
    //     and ExcessHet), one pass; also the filter summary and diagnostics ---
    SITE_FILTER(geno_filtered_vcf)
    merged_vcf = SITE_FILTER.out.filtered_vcf

    // --- 11. Pool annotation (pooled mode only) ---
    if (params.pooled) {
        POOL_ANNOTATE(merged_vcf)
        final_vcf = POOL_ANNOTATE.out.annotated_vcf
    } else {
        final_vcf = merged_vcf
    }

    // --- 12. Stats + per-sample missingness on the filtered VCF, PASS sites ---
    // Filters only flag sites, so without restricting to PASS these would
    // count failed sites too. Run on merged_vcf, before pool annotation: in
    // pooled mode every GT is then masked, which would report every sample
    // as 100% missing. Without --pooled merged_vcf is the final VCF.
    BCFTOOLS_STATS_FILTERED(merged_vcf, 'after_filter', true)
    SAMPLE_MISSINGNESS_FILTERED(merged_vcf, BCFTOOLS_STATS_FILTERED.out.stats, 'after_filter', true)

    // --- 13. MultiQC ---
    multiqc_files = BCFTOOLS_STATS_RAW.out.stats
        .mix(BCFTOOLS_STATS_FILTERED.out.stats)
        .mix(SAMPLE_MISSINGNESS_RAW.out.missingness)
        .mix(SAMPLE_MISSINGNESS_FILTERED.out.missingness)
        .mix(SITE_FILTER.out.summary)
        .collect()
    MULTIQC(multiqc_files)

    // --- 14. Cleanup intermediate work directories (opt-in) ---
    // One output per task is enough: CLEANUP empties the task directory
    // that produced it. Validated BAMs are deliberately not included.
    if (params.cleanup) {
        cleanup_targets = Channel.empty()
            .mix(
                HAPLOTYPE_CALLER.out.gvcfs.map        { sample_id, group_id, gvcf, tbi -> gvcf },
                GENOMICS_DB_IMPORT.out.genomicsdb.map { group_id, db, interval_list -> db },
                genotyped_vcfs.map                    { group_id, vcf, tbi -> vcf },
                cohort_vcf.map                        { vcf, tbi -> vcf },
                GENOTYPE_FILTER.out.filtered_vcf.map  { vcf, tbi -> vcf }
            )
            .collect()

        CLEANUP(
            MULTIQC.out.report,
            cleanup_targets,
            workflow.workDir.toString(),
            file(params.outdir).toAbsolutePath().toString()
        )
    }
}
