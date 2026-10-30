#!/usr/bin/env bash
# run_pipeline.sh — submit the GATK variant calling pipeline to LSF
#
# Usage, from your project directory (where work/, results/ and logs/ go):
#   mkdir -p logs && bsub < /path/to/nf-gatk-pipeline/run_pipeline.sh
#   bash /path/to/nf-gatk-pipeline/run_pipeline.sh     # run without LSF
#
# Edit the CONFIG section below before submitting. logs/ must exist before
# bsub: LSF opens the -o/-e files itself, before this script runs.

# ── LSF directives ────────────────────────────────────────────────────────────
# This job runs Nextflow itself for the whole pipeline: days for a large
# cohort, so it goes on 'week'. farm22 run limits: normal 12h, long 48h,
# week 7d; bsub rejects a -W longer than the queue allows.
#BSUB -J gatk_pipeline
#BSUB -q week
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -W 168:00
#BSUB -o logs/pipeline_%J.out
#BSUB -e logs/pipeline_%J.err

set -euo pipefail

# ── CONFIG — edit these before running ────────────────────────────────────────

# Where the pipeline code lives. It has to be set here: `bsub < script`
# runs LSF's spooled copy of this file, so the script can't find its own
# location.
PIPELINE_DIR="/path/to/nf-gatk-pipeline"

SAMPLESHEET="${PWD}/samplesheet.csv"            # or an absolute path elsewhere
REFERENCE="/path/to/reference.fa"

# Run data goes in the directory you submit from (LSF starts the job there),
# so each project keeps its own work dir, results and -resume history
OUTDIR="${PWD}/results"
WORK_DIR="${PWD}/work"

# LSF queues for the pipeline's own tasks are chosen per task in
# nextflow.config from its requested time (normal <= 12h, long <= 48h, else week)

# ── checks ────────────────────────────────────────────────────────────────────

MAIN_NF="${PIPELINE_DIR}/main.nf"
for f in "${MAIN_NF}" "${SAMPLESHEET}" "${REFERENCE}"; do
    [[ -f "${f}" ]] || { echo "ERROR: ${f} not found -- edit the CONFIG section of run_pipeline.sh" >&2; exit 1; }
done

mkdir -p logs

# ── environment ───────────────────────────────────────────────────────────────

# A batch shell may not have the 'module' function, and module scripts often
# reference unset variables, so relax -u while loading
set +u
source /etc/profile.d/modules.sh
module load nextflow    # or: export PATH=/path/to/nextflow:$PATH
set -u

# Nextflow's own JVM: without a heap limit it sizes the heap from the node's
# RAM, not this job's -M, and LSF kills the whole run. Keep it below -M.
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx6g}"

# ── run ───────────────────────────────────────────────────────────────────────

NF_OPTS=(
    --samplesheet  "${SAMPLESHEET}"
    --reference    "${REFERENCE}"
    --outdir       "${OUTDIR}"
)

# nextflow.config in PIPELINE_DIR is loaded automatically
nextflow run "${MAIN_NF}"        \
    -work-dir "${WORK_DIR}"      \
    "${NF_OPTS[@]}"              \
    -resume                      \
    -with-report   "${OUTDIR}/reports/nextflow_report.html"  \
    -with-timeline "${OUTDIR}/reports/nextflow_timeline.html" \
    -with-dag      "${OUTDIR}/reports/pipeline_dag.html"
