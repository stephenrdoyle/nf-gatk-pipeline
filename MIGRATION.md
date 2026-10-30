# Migration notes: per-chromosome → 10 Mb interval-scatter

## What changed

### New process: MAKE_INTERVALS
Runs once at pipeline start. Reads the `.fai` and emits one `group_NNNN.list`
file per 10 Mb window. Groups are numbered in `.fai` contig order, so
`group_0000` is the first window of chr4 (first contig in the PRJEB506 .fai),
and so on. Mitochondrion gets its own group regardless of size.

Approximate group count for H. contortus (~283 Mb):
  ~28 autosomal/X groups + 1 mitochondrion = ~29 total groups

### HAPLOTYPE_CALLER: added `--interval-padding 200`
The primary interval is still the unpadded window. The padding flag tells
HaplotypeCaller to extend its active region search 200 bp either side of the
interval boundary, ensuring variants at boundaries are called with full
haplotype context. The output GVCF still only contains records within the
padded region. GenomicsDBImport imports only the primary interval, and
GenotypeGVCFs runs with `--only-output-calls-starting-in-intervals`, so each
variant is emitted only by the group in which it starts — no duplicates at
window boundaries after concatenation.

### Channel reshape in main.nf (step 4)
This is the most important change. Previously one channel item per chromosome
fed GenomicsDBImport. Now:

  HAPLOTYPE_CALLER emits: [ sample_id, group_id, gvcf, tbi ]
      ↓ .map to groupKey(group_id, n_samples) + .groupTuple(by: 0, sort by name)
  gvcfs_by_group:          [ group_id, [gvcf×N], [tbi×N] ]
      ↓ .join with interval_by_group
  genomicsdb_input:        [ group_id, [gvcf×N], [tbi×N], interval_list ]
      ↓ GENOMICS_DB_IMPORT
      ↓ GENOTYPE_GVCFS
      ↓ collect + toSortedList
  GATHER_VCFS:             single cohort VCF

Because the `groupKey` carries the sample count, `groupTuple` releases each
group as soon as ALL samples have completed HaplotypeCaller for that group,
without waiting for the other groups. (A plain `groupTuple` with no size would
hold every group until the whole HaplotypeCaller stage had finished.) Sorting
the grouped GVCFs by name keeps the sample map and task hash stable between
runs, so `-resume` reuses completed GenomicsDBImport tasks.

### GATHER_VCFS: --naive concat replaces --allow-overlaps
With non-overlapping primary intervals concatenated in genomic order,
`bcftools concat --naive` is safe and significantly faster than
`--allow-overlaps`. The concatenated file is already sorted (groups are in
`.fai` order and each variant is emitted only by the group it starts in), so
there is no `bcftools sort`; the `##contig` lines are replaced from the `.fai`
with `bcftools reheader -h`, which copies the compressed records unchanged.
`bcftools index` fails on any out-of-order record, so it doubles as a check.

### Resource reduction
| Process           | Before       | After (base) | Reason                        |
|-------------------|--------------|--------------|-------------------------------|
| GenomicsDBImport  | 64 GB / long | 32 GB / long | 5× smaller interval           |
| GenotypeGVCFs     | 64 GB / long | 32 GB / long | 5× fewer variants per job     |
| HaplotypeCaller   | 32 GB / long | 16 GB / normal | normal queue now feasible   |

## Recovery from failure

A single corrupt GenomicsDB workspace now costs ~1/5th of a chromosome:
  - Delete the failed GENOTYPE_GVCFS work dir
  - Delete the corresponding GENOMICS_DB_IMPORT work dir
  - `nextflow run ... -resume`
  - Only that one 10 Mb group rebuilds; the other ~4 groups on that
    chromosome are already cached

## What stays the same

- `##contig` header lines rebuilt from the `.fai` in GATHER_VCFS
- JVM heap headroom formula (task.memory.toGiga() - N)
- -Dsamjdk.use_intel_deflater=false on HAPLOTYPE_CALLER
- exit -1 retry for Lustre transport failures
- --use_gnarly / --dragen_mode incompatibility guard
- VALIDATE_BAM logic

## Switching back to per-chromosome (if needed)

Set `params.window_size` to a value larger than the longest contig (~50 Mb):

  nextflow run ... --window_size 100000000

MAKE_INTERVALS will emit one group per chromosome, reproducing the old
per-chromosome behaviour with all the new fixes still in place.
