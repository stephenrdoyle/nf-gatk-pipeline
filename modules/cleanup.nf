process CLEANUP {
    /*
     * Opt-in (--cleanup): free the space used by intermediate files once the
     * pipeline has finished.
     *
     * Targets — one staged output per task, collected in main.nf:
     *   - HAPLOTYPE_CALLER        per-sample, per-window GVCFs (the bulk)
     *   - GENOMICS_DB_IMPORT      GenomicsDB workspaces
     *   - GENOTYPE_GVCFS/GNARLY   per-window genotyped VCFs
     *   - GATHER_VCFS             cohort.genotyped.vcf.gz (published copy kept)
     *   - GENOTYPE_FILTER         cohort.genotype_filtered.vcf.gz (published copy kept)
     *   - SPLIT_VARIANTS / variant filters   snps/indels VCFs
     *
     * For each target's task directory, everything is deleted except the
     * .command.* logs, so failures stay debuggable.
     * NOTE: -resume can no longer reuse these tasks — a later rerun (e.g. with
     * different filters) starts again from HaplotypeCaller.
     *
     * Safety:
     *   - each task directory is invalidated BEFORE anything in it is deleted:
     *     .exitcode is renamed to .exitcode.cleaned, so Nextflow no longer sees
     *     the task as completed. If cleanup is interrupted part-way, -resume
     *     reruns every directory it touched instead of trusting a half-deleted
     *     one (a partly removed GenomicsDB workspace still exists as a
     *     directory, so it would otherwise be reused as a valid cached output)
     *   - each staged target is followed ONE link only, to the task directory
     *     that produced it; never through a chain (a passed-through BAM links
     *     to the original), and nothing is deleted unless that directory is
     *     exactly <workDir>/xx/<hash>
     *   - publishDir copies asynchronously, so a task directory holding a file
     *     whose published copy is still a different size is skipped
     *
     * Runs locally — no LSF submission needed.
     */

    tag "cleanup"
    label 'process_low'
    executor 'local'

    input:
    path multiqc_report                   // sentinel — only runs after MultiQC completes
    path targets, stageAs: 'targets/*'    // one output (file or dir) per intermediate task
    val  work_dir                         // workflow.workDir
    val  outdir                           // absolute params.outdir

    output:
    val true, emit: done

    script:
    """
    python3 << 'PYEOF'
import os, shutil

work_root = os.path.realpath('${work_dir}')
outdir    = os.path.realpath('${outdir}')

# Sizes of every published file, by name
published = {}
for root, _, files in os.walk(outdir):
    for name in files:
        published.setdefault(name, set()).add(os.path.getsize(os.path.join(root, name)))

# Task directory of each target: follow the staging link one level only
task_dirs = set()
for name in os.listdir('targets'):
    link = os.path.join('targets', name)
    target = os.path.join(os.path.dirname(link), os.readlink(link))
    task_dirs.add(os.path.realpath(os.path.dirname(target)))

keep = lambda e: e.startswith('.command.') or e.startswith('.exitcode')
cleaned, skipped = 0, 0
for d in sorted(task_dirs):
    if os.path.dirname(os.path.dirname(d)) != work_root:
        print(f'SKIP not a task directory under {work_root}: {d}')
        skipped += 1
        continue
    entries = [e for e in os.listdir(d) if not keep(e)]
    unpublished = [
        e for e in entries
        if e in published
        and not os.path.islink(os.path.join(d, e))
        and os.path.isfile(os.path.join(d, e))
        and os.path.getsize(os.path.join(d, e)) not in published[e]
    ]
    if unpublished:
        print(f'SKIP published copy incomplete for {unpublished}: {d}')
        skipped += 1
        continue
    # Invalidate first: a cleanup interrupted after this point leaves a task
    # that -resume reruns, never a half-deleted one that looks complete
    exitcode = os.path.join(d, '.exitcode')
    if os.path.exists(exitcode):
        os.rename(exitcode, exitcode + '.cleaned')
    for e in entries:
        p = os.path.join(d, e)
        if os.path.islink(p) or not os.path.isdir(p):
            os.remove(p)
        else:
            shutil.rmtree(p)
    cleaned += 1

print(f'Cleanup complete: {cleaned} task directories emptied (logs kept), {skipped} skipped')
PYEOF
    """
}
