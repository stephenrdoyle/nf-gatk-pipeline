# nf-gatk-pipeline

A Nextflow (DSL2) pipeline for cohort germline variant calling with GATK: aligned BAMs in, one filtered multi-sample VCF out. It was built for *Haemonchus contortus* whole-genome data on the Sanger farm22 LSF cluster, but nothing in the workflow is species-specific.

The genome is split into fixed-size windows (10 Mb by default). HaplotypeCaller runs per sample per window, and joint genotyping runs on smaller pieces of each window (2 Mb by default). A failure then costs one window rather than a whole chromosome, and `-resume` reruns only what failed.

## Workflow

```
samplesheet (BAMs) ──► VALIDATE_BAM ──► HAPLOTYPE_CALLER ─────────► GENOMICS_DB_IMPORT ──► GENOTYPE_GVCFS
reference FASTA ──► .fai / .dict        (per sample × window,       (per window,           (or GNARLY_GENOTYPER;
                    ──► MAKE_INTERVALS   GVCF, 200 bp padding)        all samples)           per window)
                                                                                                  │
   ┌──────────────────────────────────────────────────────────────────────────────────────────────┘
   ▼
GATHER_VCFS ──► GENOTYPE_FILTER ──► SITE_FILTER ──────────────────► [POOL_ANNOTATE]
(cohort VCF)    (per-genotype       (hard / custom; SNP and           (--pooled only)
                 DP / GQ)            indel rules; filtered cohort VCF,
                                     summary, diagnostics)
                                          │
                                          bcftools stats, missingness, filter summary ──► MultiQC
```

Each step in more detail:

1. **`VALIDATE_BAM`** makes sure every BAM has usable read groups: one sample name, plus `PL`/`LB`/`PU`. BAMs that are already fine are linked, not copied.
2. **`MAKE_INTERVALS`** splits the reference into windows of `--window_size`. Windows never cross a contig boundary, and contigs shorter than the window (e.g. the mitochondrion) get one window each.
3. **`HAPLOTYPE_CALLER`** writes one GVCF per sample per window. It uses 200 bp of padding so variants at window edges are assembled with full context.
4. **`GENOMICS_DB_IMPORT`** builds a GenomicsDB workspace holding all samples for each joint-genotyping *piece*, a `--genotype_window_size` slice of a window (2 Mb by default). A window's pieces start as soon as all its samples have finished HaplotypeCaller, and fragments are consolidated at the end of the import (`--genomicsdb_consolidate`).
5. **`GENOTYPE_GVCFS`** joint-genotypes each piece. Each variant is emitted only by the piece it starts in, so there are no duplicates at piece edges.
6. **`GATHER_VCFS`** concatenates the pieces, in reference order, into `cohort.genotyped.vcf.gz`.
7. **`GENOTYPE_FILTER`** sets individual genotypes with low DP or GQ to no-call (`./.`).
8. **`SITE_FILTER`** flags failing sites in the `FILTER` column of `cohort.filtered.vcf.gz`, in one bcftools pass; nothing is removed. SNP rules apply to SNPs, indel rules to indels and mixed sites (SNP and indel alleles at one position). It also writes the filter summary and a diagnostics page (see [Variant filtering](#variant-filtering)).
9. **QC:** bcftools stats and per-sample missingness before filtering (all sites) and after (PASS sites only, since filters only flag sites), plus a per-filter summary, collected into a MultiQC report.

## Requirements

- Nextflow (DSL2). On farm22: `module load nextflow`.
- The software modules loaded by `process.beforeScript` in `nextflow.config`:
  - `gatk/4.6.1.0`
  - `samtools/1.21`
  - `bcftools-1.19/python-3.11.6` (also provides `python3` for the scripts in `bin/`)
  - `multiqc/1.17--pyhdfd78af_1` (loaded for the MultiQC step only)

  On farm22, GATK, samtools and MultiQC are Singularity images behind module wrappers. To run elsewhere, change these `module load` lines.
- An LSF cluster. LSF is the default executor; `-profile local` runs everything on the current machine, for small tests only.
- Write access to the reference FASTA's directory. The pipeline publishes `.fai` and `.dict` there, because GATK looks for them next to the FASTA.

## Input

**Samplesheet:** a CSV with a header row:

```csv
sample_id,bam
F2_001,/lustre/.../F2_001.bam
F2_002,/lustre/.../F2_002.bam
```

- An optional `bai` column gives the index path explicitly. Otherwise the pipeline looks for `<bam>.bai`, then `<bam without .bam>.bai`.
- `sample_id` names the sample in the output VCF, whatever the BAM's `SM` tag says.
- BAMs whose header is already valid are symlinked, not copied, so keep the original BAMs in place for as long as you might want to `-resume`.

**Reference:** a FASTA. `.fai` and `.dict` are created next to it if missing.

## Running

**On farm22, using the run script:** set `PIPELINE_DIR`, `SAMPLESHEET` and `REFERENCE` in the CONFIG section of `run_pipeline.sh`, then submit it from your project directory:

```bash
mkdir -p logs && bsub < /path/to/nf-gatk-pipeline/run_pipeline.sh
```

- **Where things go:** `work/`, `results/` and `logs/` are created in the directory you submit from, so each project keeps its own `-resume` history.
- **The submitted job runs Nextflow itself for the whole pipeline.** It goes on the `week` queue with 8 GB, and Nextflow's own Java heap is capped at 6 GB.
- **It runs with `-resume`,** so resubmitting from the same directory continues a stopped or failed run.

**Directly:**

```bash
nextflow run main.nf --samplesheet samples.csv --reference genome.fa --outdir results -resume
```

Run `nextflow` itself somewhere it can stay up for days: an LSF job, `screen`/`tmux` on a head node, or similar. The work directory defaults to `./work` in the launch directory.

**Preflight.** Before any BAM is processed, the pipeline checks its inputs and setup, and stops with one message listing every problem it found:
- **At launch:**
  - the samplesheet columns
  - unique, file-name-safe `sample_id`s
  - BAMs and indexes exist, and no index is older than its BAM
  - the reference directory is writable, if `.fai`/`.dict` need creating
- **In a short first job:**
  - the modules provide `gatk`, `samtools`, `bcftools` and `python3`
  - every BAM's contigs match the reference

## Parameters

Defaults are defined in one place: the `params` block of `nextflow.config`.

| Parameter | Default | Meaning |
|---|---|---|
| `--samplesheet` | *(required)* | Samplesheet CSV |
| `--reference` | *(required)* | Reference FASTA |
| `--outdir` | `results` | Output directory |
| `--window_size` | `10000000` | HaplotypeCaller window size in bp. Larger than the longest contig = one window per contig |
| `--genotype_window_size` | `2000000` | GenomicsDBImport/GenotypeGVCFs piece size in bp, within each window. ≥ `--window_size` = one piece per window. Can be changed without rerunning HaplotypeCaller |
| `--genomicsdb_consolidate` | `true` | Merge GenomicsDB's per-batch fragments after import, so GenotypeGVCFs reads one |
| `--genomicsdb_max_forks` | `50` | Most GenomicsDBImport jobs running at once; each reads every sample's GVCF from Lustre, so this limits the I/O burst |
| `--filter_mode` | `hard` | `hard` (GATK's thresholds) or `custom` (see [Variant filtering](#variant-filtering)) |
| `--snp_filters` / `--indel_filters` | — | Custom mode: comma-separated `<INFO tag> <op> <number>` expressions, e.g. `"QD < 2.0,FS > 60.0"` |
| `--max_depth` | — | Flag sites with INFO/DP above this: a multiple of the median site depth (`2x`) or a number. Off unless set |
| `--max_missing` | — | Flag sites where more than this fraction of samples has no call after the genotype filter, e.g. `0.2`. Off unless set |
| `--excess_het` | `false` | Flag sites with too many heterozygotes: `true` (ExcessHet > 54.69) or a number. Not for F1 crosses |
| `--genotype_dp_threshold` | `5` | Genotypes with `DP` below this are set to no-call |
| `--genotype_gq_threshold` | `20` | Genotypes with `GQ` below this are set to no-call (not in pooled or Gnarly mode) |
| `--heterozygosity` | — | GenotypeGVCFs' heterozygosity prior: unset = GATK default (0.001), a number, or `auto` (estimated from the BAMs, see [Modes](#modes)) |
| `--indel_heterozygosity` | — | Indel prior; unset = `--heterozygosity` × 0.125 |
| `--heterozygosity_stdev` | — | GATK's `--heterozygosity-stdev`; unset = √(heterozygosity ÷ 10), which keeps the prior's strength at GATK's default so only its mean moves |
| `--genotype_rgq_threshold` | — | Gnarly mode: genotypes with `RGQ` below this are set to no-call |
| `--use_gnarly` | `false` | Joint-genotype with GnarlyGenotyper instead of GenotypeGVCFs |
| `--pooled` | `false` | Samples are pools, not individuals (see [Modes](#modes)) |
| `--cleanup` | `false` | Empty intermediate work directories at the end (see [Cleanup](#cleanup)) |
| `--publish_mode` | `link` | How result VCFs are published: `link` (hard link: no extra space, index timestamps stay correct; needs `--outdir` on the same filesystem as `work/`, checked at launch) or `copy` |

## Variant filtering

Filtering never removes sites. Each site in `cohort.filtered.vcf.gz` gets `PASS` or the names of the filters it failed. To keep only passing sites:

```bash
bcftools view -f PASS results/cohort.filtered.vcf.gz
```

**`hard` (default):** GATK's recommended thresholds.

| SNPs | Indels and mixed sites |
|---|---|
| `QD < 2.0`, `MQ < 40.0`, `FS > 60.0`, `SOR > 3.0`, `MQRankSum < -12.5`, `ReadPosRankSum < -8.0` | `QD < 2.0`, `FS > 200.0`, `SOR > 10.0`, `ReadPosRankSum < -20.0` |


A site missing an annotation passes that filter, as in GATK.

**`custom`:** your own thresholds, via `--snp_filters` and `--indel_filters` (both required). Each is a comma-separated list of `<INFO tag> <op> <number>`, with `<`, `<=`, `>`, `>=`, `==` or `!=`, e.g. `--snp_filters "QD < 2.0,MQ < 40.0,FS > 60.0"`. Filter names are built from the expression (`SNP_QD_lt_2.0`).

**Optional site filters**, in either mode, all off by default:

| Option | Flags | Catches | Caution |
|---|---|---|---|
| `--max_depth 2x` | `HighDepth`: INFO/DP above 2 × the median site depth (SNPs and indels separately), or a fixed number | Collapsed repeats and paralogs | Depth varies along some genomes; check the DP histogram first |
| `--max_missing 0.2` | `HighMissing`: more than 20% of samples uncalled after the genotype filter | Sites genotyped in too few samples | How much to allow depends on the analysis. On the F2 cohort 0.2 flagged 30% of SNPs, so check per-sample missingness first: a few poor samples make every site look incomplete |
| `--excess_het` | `ExcessHet`: GATK's ExcessHet above 54.69 | Paralogs (nearly every sample heterozygous) | **Never for F1 crosses**: every site that differs between the parents is heterozygous in all offspring |

**Checking the thresholds suit your cohort.** The thresholds are fixed, so the same PASS definition holds across studies. `filter_reports/filter_diagnostics.html` shows whether they fit each cohort. For every annotation it shows a histogram with the thresholds marked, the share of sites flagged, and the sites flagged by that filter alone. It also shows depth, missingness and ExcessHet at several levels whether or not they are applied, so you can see what each option would do before turning it on. A warning appears when one filter flags more than 10% of sites. If a threshold sits inside the bulk of a distribution, tune it with `custom`.

**Why there is no `optimised` mode any more.** It set each threshold at the worst 5% of the cohort's own values. On the 200-sample F2 cohort that flagged 21% of SNPs. It flagged mostly real sites (Ts/Tv of the removed sites 1.6–2.1, against 2.00 for PASS), and on QD and SOR it was looser than GATK's thresholds. GATK's thresholds kept 2.3 M more SNPs at essentially the same Ts/Tv (1.98). The comparison script is in `benchmarks/filter_comparison/` to repeat this on another cohort.

## Modes

**GnarlyGenotyper (`--use_gnarly`):** a faster alternative to GenotypeGVCFs for large cohorts. GQ isn't emitted in this mode, so use `--genotype_rgq_threshold` instead of `--genotype_gq_threshold`.

**Estimated heterozygosity prior (`--heterozygosity auto`):** GATK's genotyping prior assumes human-like diversity (0.001). Many species, and crosses such as F2 cohorts, are 10–30× more diverse. With `auto`, the pipeline estimates heterozygosity from the BAMs while HaplotypeCaller runs:
- **How:** ANGSD per-sample site frequency spectra over a fixed ~20 Mb sample of the genome. There are no variant calls, so no caller prior is involved.
- **Which samples:** the first 20 samples by `sample_id`. Low-coverage samples are excluded.
- **The value:** the median, which is passed to GenotypeGVCFs and recorded with the per-sample values in `pipeline_info/heterozygosity_estimate.txt`.
- **Software:** ANGSD runs from a pinned Biocontainers image via Singularity (`module load ISG/singularity/3.11.4`). Images are cached in `$NXF_SINGULARITY_CACHEDIR`, which `run_pipeline.sh` sets.
- **Not used** with `--pooled` or `--use_gnarly`.
- **What reruns:** changing the prior reruns joint genotyping onwards, but not HaplotypeCaller.

**Pooled samples (`--pooled`):** for samples that are pools of individuals.
- **Genotypes are masked.** Variants are still called as if each sample were diploid, so every `GT` in the final VCF is set to missing. Use `FORMAT/AD`, or the added `FORMAT/AF` (alt reads ÷ total, one value per alt allele), instead.
- **The GQ filter is skipped.**
- **Low-frequency alleles may be missed,** because the diploid calling model is built for individuals rather than pools.

The output is `cohort.pooled_annotated.vcf.gz`. Stats and missingness are computed on the VCF before masking.

## Output

```
results/
├── variants/
│   └── cohort.genotyped.vcf.gz            joint-genotyped cohort, unfiltered
├── cohort.genotype_filtered.vcf.gz        after per-genotype DP/GQ filtering
├── cohort.filtered.vcf.gz                 final: sites flagged PASS or by filter name
├── cohort.pooled_annotated.vcf.gz         --pooled only: GT masked, FORMAT/AF added
├── filter_reports/
│   ├── filter_summary.tsv                 sites per filter: all, PASS, each filter, "only this filter"
│   └── filter_diagnostics.{html,tsv}      each threshold against this cohort's distributions
├── bcftools_stats/                        before_filter (all sites) / after_filter (PASS sites)
├── missingness/                           per-sample missing and low-DP genotypes:
│                                          before_filter (all sites) / after_filter (PASS sites)
├── multiqc/multiqc_report.html            QC summary
├── pipeline_info/software_versions.txt    tool and module versions, git commit, command line
├── pipeline_trace_<timestamp>.txt         one row per task attempt (see below)
├── pipeline_report.html
└── pipeline_timeline.html
```

Every VCF has a `.tbi` index alongside it. Result VCFs are hard links to the files in `work/` (`--publish_mode link`), so they take no extra space. `run_pipeline.sh` writes the Nextflow report and timeline (plus a DAG) to `results/reports/` instead of the top level.

## Resources, retries and failures

**Resources and queues.** Each step's CPUs, memory and time are set in `nextflow.config`. The LSF queue is chosen from the time a task asks for: `normal` up to 12h, `long` up to 48h, otherwise `week`. GATK's Java heap is set from the task's memory.

**Retries.** A task that fails with one of these exit codes is retried up to 3 times, with memory and time multiplied by the attempt number:

| Exit code | Usual cause |
|---|---|
| none (`-`) | lost exit status, e.g. a Lustre eviction |
| −1 | exit status couldn't be read |
| 130 | LSF memory limit |
| 140 | LSF run limit |
| 137, 143 | killed |
| 134, 135, 139 | crash (abort, SIGBUS, segfault) |
| 250 | GenomicsDB failure |
| 1 | GenomicsDBImport and GenotypeGVCFs only, where a SIGBUS can show up as exit 1 |

Any other failure, or a retryable one that keeps failing, stops the run with an error. Failures are never ignored, so a sample or window can't silently go missing.

**Finding what failed:**

```bash
awk -F'\t' 'NR==1 || ($5 != "COMPLETED" && $5 != "CACHED")' results/pipeline_trace_*.txt
```

```bash
nextflow log last -f name,status,exit,attempt,workdir -F 'status == "FAILED"'
```

Each task's own logs are in its work directory: `.command.err`, `.command.log` and `.exitcode`. After fixing the cause, rerun the same command with `-resume`.

**What `-resume` reruns.** Nextflow's cache key includes each step's script and inputs. Editing `HAPLOTYPE_CALLER`, or anything that runs before it (`VALIDATE_BAM`, `SAMTOOLS_FAIDX`, `CREATE_DICT`), reruns every HaplotypeCaller job. Changes to the software modules are *not* part of the cache key; `pipeline_info/software_versions.txt` records what each run actually used.

## Cleanup

`--cleanup` empties the work directories of large intermediates once MultiQC has finished:
- the per-sample GVCFs
- the GenomicsDB workspaces
- the per-window and intermediate VCFs

It keeps each directory's `.command.*` logs. It only touches task directories inside the work directory, and skips any whose published copy hasn't finished copying.

After `--cleanup`, `-resume` can't reuse those steps: changing a filter later means starting again from HaplotypeCaller. Use it once you're happy with the cohort VCF.

## Repository layout

```
main.nf               workflow
nextflow.config       parameters, resources, retries, LSF executor, reports
run_pipeline.sh       LSF submission script for farm22
modules/              one process per file
bin/                  scripts used by processes (filter thresholds and summary, pool annotation)
assets/               MultiQC configuration
samplesheet.csv       example samplesheet
MIGRATION.md          notes on the move from per-chromosome to 10 Mb interval scatter
```
