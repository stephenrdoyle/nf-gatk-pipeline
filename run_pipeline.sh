#!/usr/bin/env bash
# run_pipeline.sh — submit the GATK variant calling pipeline to LSF
#
# Usage:
#   bsub < run_pipeline.sh          # submit as an LSF job from anywhere
#   bash run_pipeline.sh            # run directly
#
# The script locates main.nf and nextflow.config relative to its own location,
# so it works correctly regardless of which directory you run it from.
#
# Edit the CONFIG section below before submitting.

# ── LSF directives ────────────────────────────────────────────────────────────
#BSUB -J gatk_pipeline
#BSUB -q normal
#BSUB -n 1
#BSUB -M 4096
#BSUB -R "rusage[mem=4096]"
#BSUB -W 72:00
#BSUB -o logs/pipeline_%J.out
#BSUB -e logs/pipeline_%J.err

set -euo pipefail

# ── resolve pipeline directory ────────────────────────────────────────────────
# Works whether the script is called by path, sourced, or submitted via bsub.

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_NF="${PIPELINE_DIR}/main.nf"
NF_CONFIG="${PIPELINE_DIR}/nextflow.config"

# ── CONFIG — edit these before running ────────────────────────────────────────

SAMPLESHEET="${PIPELINE_DIR}/samplesheet.csv"   # or an absolute path elsewhere
REFERENCE="/path/to/hg38.fa"
OUTDIR="${PIPELINE_DIR}/results"
# LSF queues are set per process in nextflow.config (normal/long by label)

# Nextflow work dir — set explicitly so it's always next to the pipeline,
# not in whatever directory you happened to submit from
WORK_DIR="${PIPELINE_DIR}/work"

# ── environment ───────────────────────────────────────────────────────────────

for f in "${SAMPLESHEET}" "${REFERENCE}"; do
    [[ -f "${f}" ]] || { echo "ERROR: ${f} not found -- edit the CONFIG section of $0" >&2; exit 1; }
done

mkdir -p logs

module load nextflow    # or: export PATH=/path/to/nextflow:$PATH

# ── run ───────────────────────────────────────────────────────────────────────

NF_OPTS=(
    --samplesheet  "${SAMPLESHEET}"
    --reference    "${REFERENCE}"
    --outdir       "${OUTDIR}"
)

nextflow run "${MAIN_NF}"        \
    -config  "${NF_CONFIG}"      \
    -work-dir "${WORK_DIR}"      \
    "${NF_OPTS[@]}"              \
    -resume                      \
    -with-report   "${OUTDIR}/reports/nextflow_report.html"  \
    -with-timeline "${OUTDIR}/reports/nextflow_timeline.html" \
    -with-dag      "${OUTDIR}/reports/pipeline_dag.html"
