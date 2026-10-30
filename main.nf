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
 *   - Post-calling: genotype filter -> split -> variant filter -> merge -> QC
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
include { SPLIT_VARIANTS      } from './modules/split_variants'
include { HARD_FILTER         } from './modules/hard_filter'
include { OPTIMISED_FILTER    } from './modules/optimised_filter'
include { CUSTOM_FILTER       } from './modules/custom_filter'
include { MERGE_FILTERED      } from './modules/merge_filtered'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_RAW      } from './modules/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_FILTERED } from './modules/bcftools_stats'
include { SAMPLE_MISSINGNESS as SAMPLE_MISSINGNESS_RAW      } from './modules/sample_missingness'
include { SAMPLE_MISSINGNESS as SAMPLE_MISSINGNESS_FILTERED } from './modules/sample_missingness'
include { FILTER_SUMMARY      } from './modules/filter_summary'
include { MULTIQC             } from './modules/multiqc'
include { CLEANUP             } from './modules/cleanup'
include { POOL_ANNOTATE       } from './modules/pool_annotate'
include { PREFLIGHT           } from './modules/preflight'
include { SOFTWARE_VERSIONS   } from './modules/software_versions'

// ---------------------------------------------------------------------------
// Parameter defaults: defined once, in the params block of nextflow.config
// ---------------------------------------------------------------------------

// Guard: required inputs
if (!params.samplesheet || !params.reference) {
    error "ERROR: --samplesheet and --reference are required"
}

// Guard: --use_gnarly and --dragen_mode are incompatible
if (params.use_gnarly && params.dragen_mode) {
    error "ERROR: --use_gnarly and --dragen_mode cannot be used together"
}

// Guard: custom filter mode requires filter expressions
if (params.filter_mode == 'custom' && (!params.snp_filters || !params.indel_filters)) {
    error "ERROR: --filter_mode custom requires --snp_filters and --indel_filters"
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

    // --- 0. Reference files ---
    reference = file(params.reference, checkIfExists: true)

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
        .map { key, gvcfs, tbis -> [ key.getGroupTarget(), gvcfs, tbis ] }

    interval_by_group = interval_ch
        .map { list_file -> [ list_file.baseName, list_file ] }

    genomicsdb_input = gvcfs_by_group
        .join(interval_by_group, by: 0)

    // --- 5. GenomicsDBImport: per-interval, all samples ---
    GENOMICS_DB_IMPORT(
        genomicsdb_input,
        reference,
        reference_fai,
        reference_dict
    )

    // --- 6. Genotype: GenotypeGVCFs or GnarlyGenotyper per-interval ---
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
            reference,
            reference_fai,
            reference_dict
        )
        genotyped_vcfs = GENOTYPE_GVCFS.out.vcf
    }

    // --- 7. Gather all interval VCFs into final cohort VCF ---
    all_vcfs = genotyped_vcfs
        .map { group_id, vcf, tbi -> [ group_id, vcf, tbi ] }
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
    GENOTYPE_FILTER(cohort_vcf, reference, reference_fai, reference_dict)
    geno_filtered_vcf = GENOTYPE_FILTER.out.filtered_vcf

    // --- 10. Split into SNPs and INDELs ---
    SPLIT_VARIANTS(geno_filtered_vcf, reference, reference_fai, reference_dict)
    // Emits: snps [ vcf, tbi ] and indels [ vcf, tbi ]

    // Tag with variant type for downstream channel handling
    snps_ch   = SPLIT_VARIANTS.out.snps.map   { vcf, tbi -> [ 'snps',   vcf, tbi ] }
    indels_ch = SPLIT_VARIANTS.out.indels.map { vcf, tbi -> [ 'indels', vcf, tbi ] }
    split_ch  = snps_ch.mix(indels_ch)

    // --- 11. Variant-level filter (hard / optimised / custom) ---
    if (params.filter_mode == 'optimised') {
        OPTIMISED_FILTER(split_ch, reference, reference_fai, reference_dict)
        filtered_ch = OPTIMISED_FILTER.out.filtered_vcf
    } else if (params.filter_mode == 'custom') {
        CUSTOM_FILTER(split_ch, reference, reference_fai, reference_dict)
        filtered_ch   = CUSTOM_FILTER.out.filtered_vcf
    } else {
        // default: hard filter
        HARD_FILTER(split_ch, reference, reference_fai, reference_dict)
        filtered_ch   = HARD_FILTER.out.filtered_vcf
    }

    // Separate filtered SNP and INDEL channels for merging
    filt_snps_ch   = filtered_ch.filter { it[0] == 'snps'   }.map { vtype, vcf, tbi -> [ vcf, tbi ] }
    filt_indels_ch = filtered_ch.filter { it[0] == 'indels' }.map { vtype, vcf, tbi -> [ vcf, tbi ] }

    // --- 12. Merge filtered SNPs + INDELs ---
    MERGE_FILTERED(
        filt_snps_ch,
        filt_indels_ch,
        reference,
        reference_fai,
        reference_dict
    )
    merged_vcf = MERGE_FILTERED.out.filtered_vcf

    // --- 13. Pool annotation (pooled mode only) ---
    if (params.pooled) {
        POOL_ANNOTATE(merged_vcf)
        final_vcf = POOL_ANNOTATE.out.annotated_vcf
    } else {
        final_vcf = merged_vcf
    }

    // --- 14. Stats + per-sample missingness on the filtered VCF, PASS sites ---
    // Filters only flag sites, so without restricting to PASS these would
    // count failed sites too. Run on merged_vcf, before pool annotation: in
    // pooled mode every GT is then masked, which would report every sample
    // as 100% missing. Without --pooled merged_vcf is the final VCF.
    BCFTOOLS_STATS_FILTERED(merged_vcf, 'after_filter', true)
    SAMPLE_MISSINGNESS_FILTERED(merged_vcf, BCFTOOLS_STATS_FILTERED.out.stats, 'after_filter', true)

    // --- 15. Filter summary: sites per filter, from the filtered VCFs ---
    FILTER_SUMMARY(filt_snps_ch, filt_indels_ch)

    // --- 16. MultiQC ---
    multiqc_files = BCFTOOLS_STATS_RAW.out.stats
        .mix(BCFTOOLS_STATS_FILTERED.out.stats)
        .mix(SAMPLE_MISSINGNESS_RAW.out.missingness)
        .mix(SAMPLE_MISSINGNESS_FILTERED.out.missingness)
        .mix(FILTER_SUMMARY.out.summary)
        .collect()
    MULTIQC(multiqc_files)

    // --- 17. Cleanup intermediate work directories (opt-in) ---
    // One output per task is enough: CLEANUP empties the task directory
    // that produced it. Validated BAMs are deliberately not included.
    if (params.cleanup) {
        cleanup_targets = Channel.empty()
            .mix(
                HAPLOTYPE_CALLER.out.gvcfs.map        { sample_id, group_id, gvcf, tbi -> gvcf },
                GENOMICS_DB_IMPORT.out.genomicsdb.map { group_id, db, interval_list -> db },
                genotyped_vcfs.map                    { group_id, vcf, tbi -> vcf },
                cohort_vcf.map                        { vcf, tbi -> vcf },
                GENOTYPE_FILTER.out.filtered_vcf.map  { vcf, tbi -> vcf },
                SPLIT_VARIANTS.out.snps.map           { vcf, tbi -> vcf },
                SPLIT_VARIANTS.out.indels.map         { vcf, tbi -> vcf },
                filt_snps_ch.map                      { vcf, tbi -> vcf },
                filt_indels_ch.map                    { vcf, tbi -> vcf }
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
