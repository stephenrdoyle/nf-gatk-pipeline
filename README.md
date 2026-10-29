# nf-gatk-pipeline

A Nextflow (DSL2) pipeline for cohort germline variant calling with GATK: aligned BAMs in, one filtered multi-sample VCF out. It was built for *Haemonchus contortus* whole-genome data on the Sanger farm22 LSF cluster, but nothing in the workflow is species-specific.

The genome is split into fixed-size windows (10 Mb by default). HaplotypeCaller runs per sample per window, and joint genotyping runs per window. A failure then costs one window rather than a whole chromosome, and `-resume` reruns only what failed.

## Workflow

```
samplesheet (BAMs) ──► VALIDATE_BAM ──► HAPLOTYPE_CALLER ─────────► GENOMICS_DB_IMPORT ──► GENOTYPE_GVCFS
reference FASTA ──► .fai / .dict        (per sample × window,       (per window,           (or GNARLY_GENOTYPER;
                    ──► MAKE_INTERVALS   GVCF, 200 bp padding)        all samples)           per window)
                                                                                                  │
   ┌──────────────────────────────────────────────────────────────────────────────────────────────┘
   ▼
GATHER_VCFS ──► GENOTYPE_FILTER ──► SPLIT_VARIANTS ──► SNP / indel filters ──► MERGE_FILTERED ──► [POOL_ANNOTATE]
(cohort VCF)    (per-genotype       (SNPs; indels       (hard / optimised /     (filtered cohort    (--pooled only)
                 DP / GQ)            + mixed sites)      custom)                  VCF)
                                                                                     │
                                          bcftools stats, missingness, filter summary ──► MultiQC
```

Each step in more detail:

1. **`VALIDATE_BAM`** makes sure every BAM has usable read groups: one sample name, plus `PL`/`LB`/`PU`. BAMs that are already fine are linked, not copied.
2. **`MAKE_INTERVALS`** splits the reference into windows of `--window_size`. Windows never cross a contig boundary, and contigs shorter than the window (e.g. the mitochondrion) get one window each.
3. **`HAPLOTYPE_CALLER`** writes one GVCF per sample per window. It uses 200 bp of padding so variants at window edges are assembled with full context.
4. **`GENOMICS_DB_IMPORT`** builds one GenomicsDB workspace per window, holding all samples. Each window starts as soon as all its samples have finished HaplotypeCaller.
5. **`GENOTYPE_GVCFS`** joint-genotypes each window. Each variant is emitted only by the window it starts in, so there are no duplicates at window edges.
6. **`GATHER_VCFS`** concatenates the windows, in reference order, into `cohort.genotyped.vcf.gz`.
7. **`GENOTYPE_FILTER`** sets individual genotypes with low DP or GQ to no-call (`./.`).
8. **`SPLIT_VARIANTS`** separates SNPs from indels. Mixed sites (SNP and indel alleles at one position) go with the indels.
9. **The SNP and indel filters** flag failing sites in the `FILTER` column; nothing is removed. **`MERGE_FILTERED`** recombines them into `cohort.filtered.vcf.gz`.
10. **QC:** bcftools stats and per-sample missingness before filtering (all sites) and after (PASS sites only, since filters only flag sites), plus a per-filter summary, collected into a MultiQC report.

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

## Parameters

Defaults are defined in one place: the `params` block of `nextflow.config`.

| Parameter | Default | Meaning |
|---|---|---|
| `--samplesheet` | *(required)* | Samplesheet CSV |
| `--reference` | *(required)* | Reference FASTA |
| `--outdir` | `results` | Output directory |
| `--window_size` | `10000000` | Window size in bp for the interval scatter. Larger than the longest contig = one window per contig |
| `--filter_mode` | `hard` | `hard`, `optimised` or `custom` (see [Variant filtering](#variant-filtering)) |
| `--filter_percentile_lo` | `5` | Optimised mode: percentile cut-off |
| `--snp_filters` / `--indel_filters` | — | Custom mode: comma-separated GATK expressions, e.g. `"QD < 2.0,FS > 60.0"` |
| `--genotype_dp_threshold` | `5` | Genotypes with `DP` below this are set to no-call |
| `--genotype_gq_threshold` | `20` | Genotypes with `GQ` below this are set to no-call (not in pooled or Gnarly mode) |
| `--genotype_rgq_threshold` | — | Gnarly mode: genotypes with `RGQ` below this are set to no-call |
| `--use_gnarly` | `false` | Joint-genotype with GnarlyGenotyper instead of GenotypeGVCFs |
| `--dragen_mode` | `false` | DRAGEN-GATK calling (see [Modes](#modes)) |
| `--pooled` | `false` | Samples are pools, not individuals (see [Modes](#modes)) |
| `--cleanup` | `false` | Empty intermediate work directories at the end (see [Cleanup](#cleanup)) |

## Variant filtering

Filtering never removes sites. Each site in `cohort.filtered.vcf.gz` gets `PASS` or the names of the filters it failed. To keep only passing sites:

```bash
bcftools view -f PASS results/cohort.filtered.vcf.gz
```

**`hard` (default):** GATK's recommended thresholds.

| SNPs | Indels and mixed sites |
|---|---|
| `QD < 2.0`, `MQ < 40.0`, `FS > 60.0`, `SOR > 3.0`, `MQRankSum < -12.5`, `ReadPosRankSum < -8.0` | `QD < 2.0`, `FS > 200.0`, `SOR > 10.0`, `ReadPosRankSum < -20.0` |

With `--dragen_mode`, the single DRAGEN filter `QUAL < 10.4139` replaces both columns.

**`optimised`:** thresholds are derived from your data. For each annotation, the worst `--filter_percentile_lo`% of sites are flagged: the bottom tail for `QD`, `MQ`, `MQRankSum` and `ReadPosRankSum`, and the top tail for `FS` and `SOR`. This flags at least that fraction per annotation by design, so check the HTML report in `filter_reports/` to see where each threshold landed.

**`custom`:** your own expressions, via `--snp_filters` and `--indel_filters`. Both are required in this mode.

## Modes

**GnarlyGenotyper (`--use_gnarly`):** a faster alternative to GenotypeGVCFs for large cohorts. GQ isn't emitted in this mode, so use `--genotype_rgq_threshold` instead of `--genotype_gq_threshold`. It can't be combined with `--dragen_mode`.

**DRAGEN-GATK (`--dragen_mode`):**
1. builds an STR table from the reference
2. calibrates a DRAGstr model for each sample
3. runs HaplotypeCaller with `--dragen-mode` and that model
4. in `hard` filter mode, applies DRAGEN's QUAL filter

The QUAL threshold was calibrated on human data, so compare it with the standard hard filters on a subset before relying on it.

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
│   └── {snps,indels}_filter_report.html   optimised mode: thresholds and histograms
├── bcftools_stats/                        before_filter (all sites) / after_filter (PASS sites)
├── missingness/                           per-sample missing and low-DP genotypes:
│                                          before_filter (all sites) / after_filter (PASS sites)
├── multiqc/multiqc_report.html            QC summary
├── pipeline_info/software_versions.txt    tool and module versions, git commit, command line
├── pipeline_trace_<timestamp>.txt         one row per task attempt (see below)
├── pipeline_report.html
└── pipeline_timeline.html
```

Every VCF has a `.tbi` index alongside it. `run_pipeline.sh` writes the Nextflow report and timeline (plus a DAG) to `results/reports/` instead of the top level.

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
modules/              one process per file (dragstr.nf holds the two DRAGEN calibration steps)
bin/                  scripts used by processes (filter thresholds and summary, pool annotation)
assets/               MultiQC configuration
samplesheet.csv       example samplesheet
MIGRATION.md          notes on the move from per-chromosome to 10 Mb interval scatter
```
