#!/usr/bin/env bash
# Compare site-filtering schemes on an existing cohort.filtered.vcf.gz.
# Read-only: writes into OUTDIR, never touches the VCF or the pipeline's work/.
#
# Usage (farm22):
#   bsub -q normal -n 4 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        -J filter_cmp -o filter_cmp_%J.out -e filter_cmp_%J.err \
#        bash run_comparison.sh <cohort.filtered.vcf.gz> <outdir>
#
# Writes <outdir>/sites.tsv.gz (one line per site; delete when done) and
# <outdir>/filter_comparison.md (the result).
set -euo pipefail

vcf=$(readlink -f "$1")
outdir=$2
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

source /etc/profile.d/modules.sh
module load bcftools-1.19/python-3.11.6

mkdir -p "$outdir"
cd "$outdir"

# AC/AN/F_MISSING recomputed from the genotypes as they are after the
# genotype filter (the INFO values GATK wrote predate it)
if [ ! -s sites.tsv.gz ]; then
    bcftools +fill-tags --threads 2 "$vcf" -Ou -- -t AC,AN,F_MISSING \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/DP\t%INFO/ExcessHet\t%INFO/F_MISSING\t%INFO/AC\t%INFO/AN\n' \
        | gzip -1 > sites.tsv.gz.tmp
    mv sites.tsv.gz.tmp sites.tsv.gz
fi
echo "sites: $(gzip -dc sites.tsv.gz | wc -l)"

python3 "$here/compare_filters.py" sites.tsv.gz > filter_comparison.md
echo "written: $(pwd)/filter_comparison.md"
