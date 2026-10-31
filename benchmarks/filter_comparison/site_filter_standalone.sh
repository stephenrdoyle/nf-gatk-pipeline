#!/usr/bin/env bash
# Run the pipeline's SITE_FILTER step on an existing genotype-filtered cohort
# VCF, outside Nextflow: the same commands as modules/site_filter.nf, for a
# cohort run with older code (re-filtering it through -resume would rerun
# joint genotyping). Writes into OUTDIR; the input VCF is only read.
#
# Usage (farm22):
#   bsub -q long -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        -J site_filter -o site_filter_%J.out -e site_filter_%J.err \
#        bash site_filter_standalone.sh <cohort.genotype_filtered.vcf.gz> <outdir> [site_filter.py options]
#
# Options are passed to site_filter.py apply, e.g. --max-depth 2x --max-missing 0.2
# --excess-het true; default is hard mode with nothing optional applied.
set -euo pipefail

vcf=$(readlink -f "$1")
outdir=$2
shift 2
bin=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../bin" && pwd)
threads=${LSB_DJOB_NUMPROC:-4}

source /etc/profile.d/modules.sh
module load bcftools-1.19/python-3.11.6

mkdir -p "$outdir"
cd "$outdir"
tags=QD,FS,SOR,MQ,MQRankSum,ReadPosRankSum,DP,ExcessHet
fmt='%CHROM\t%POS\t%REF\t%ALT\t%INFO/F_MISSING'
for t in ${tags//,/ }; do fmt="$fmt\t%INFO/$t"; done

bcftools view --header-only "$vcf" | grep '^##contig' > contigs.txt || true
if [ ! -s sites.tsv.gz ]; then
    bcftools +fill-tags --threads "$threads" "$vcf" -Ou -- -t F_MISSING \
        | bcftools query -f "$fmt\n" | gzip -1 > sites.tsv.gz.tmp
    mv sites.tsv.gz.tmp sites.tsv.gz
fi

python3 "$bin/site_filter.py" apply --sites sites.tsv.gz --columns "$tags" --contigs contigs.txt --mode hard "$@" \
    | bcftools view -Oz -o decisions.vcf.gz
bcftools index -f --tbi decisions.vcf.gz

bcftools annotate --threads "$threads" --annotations decisions.vcf.gz --columns FILTER \
    -Oz -o cohort.filtered.vcf.gz "$vcf"
bcftools index -f --tbi --threads "$threads" cohort.filtered.vcf.gz

bcftools query -f '%FILTER\n' cohort.filtered.vcf.gz > observed_filters.txt
python3 "$bin/site_filter.py" verify expected_filters.tsv observed_filters.txt
rm observed_filters.txt
echo "done: $(pwd)/cohort.filtered.vcf.gz, filter_summary.tsv, filter_diagnostics.html"
