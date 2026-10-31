#!/usr/bin/env Rscript
# Population-structure QC with SNPRelate: nuclear PCA, KING-robust kinship
# and, if a mitochondrial SNP VCF is given, a mitochondrial PCA.
#
# Usage:
#   qc_snprelate.R <nuclear_snps.vcf.gz> <mito_snps.vcf.gz | *.NONE> <ld_threshold> <threads>
#
# Inputs are already reduced by QC_EXTRACT: PASS biallelic SNPs, minor allele
# frequency >= --qc_maf, thinned to one site per kb, with the sex and
# mitochondrial contigs excluded from the nuclear set.
#
# Writes (tab-separated):
#   pca_nuclear.tsv        sample, PC1..PCn (n <= 10)
#   pca_nuclear_eigen.tsv  pc, variance_explained
#   kinship.tsv            ID1, ID2, IBS0, kinship   (KING-robust, all pairs)
#   pca_mito.tsv, pca_mito_eigen.tsv   only with a mito VCF
#   snp_counts.tsv         set, n_snps   (input and after LD pruning)
#   notes.txt              anything skipped, e.g. no polymorphic mito sites
#
# Chromosome codes are arbitrary contig names, so every SNPRelate call uses
# autosome.only = FALSE; which contigs count as nuclear is decided upstream.

suppressPackageStartupMessages(library(SNPRelate))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: qc_snprelate.R <nuclear.vcf.gz> <mito.vcf.gz|*.NONE> <ld_threshold> <threads>")
nuclear_vcf <- args[1]
mito_vcf    <- args[2]
ld          <- as.numeric(args[3])
threads     <- as.integer(args[4])

# LD pruning visits SNPs in random order: fix the seed so a rerun on the same
# input gives the same SNP set, PCs and kinship
set.seed(1)

counts <- list()
notes  <- character()

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_pca <- function(pca, prefix) {
  n   <- min(10, ncol(pca$eigenvect))
  tab <- data.frame(sample = pca$sample.id, pca$eigenvect[, seq_len(n), drop = FALSE])
  names(tab)[-1] <- paste0("PC", seq_len(n))
  write_tsv(tab, paste0(prefix, ".tsv"))
  write_tsv(data.frame(pc = paste0("PC", seq_len(n)),
                       variance_explained = pca$varprop[seq_len(n)]),
            paste0(prefix, "_eigen.tsv"))
}

n_snps <- function(gds) length(read.gdsn(index.gdsn(gds, "snp.id")))

# Too few SNPs (a tiny test, or a near-clonal cohort) is reported, not fatal:
# the cohort VCF is already finished by the time this runs, and the report
# should still be written. Empty tables keep the outputs the same shape.
empty_outputs <- function() {
  write_tsv(data.frame(sample = character()), "pca_nuclear.tsv")
  write_tsv(data.frame(pc = character(), variance_explained = numeric()), "pca_nuclear_eigen.tsv")
  write_tsv(data.frame(ID1 = character(), ID2 = character(), IBS0 = numeric(), kinship = numeric()), "kinship.tsv")
}

# ── nuclear: LD pruning, PCA, kinship ───────────────────────────────────────
snpgdsVCF2GDS(nuclear_vcf, "nuclear.gds", method = "biallelic.only", verbose = FALSE)
g <- snpgdsOpen("nuclear.gds")
counts$nuclear_input <- n_snps(g)
if (counts$nuclear_input < 100) {
  notes <- c(notes, sprintf("Nuclear PCA and kinship skipped: only %d SNPs reached the QC set (100 needed).",
                            counts$nuclear_input))
  empty_outputs()
} else {
  pruned <- snpgdsLDpruning(g, ld.threshold = ld, method = "corr", slide.max.bp = 500000,
                            autosome.only = FALSE, num.thread = threads, verbose = FALSE)
  snps <- unlist(unname(pruned))
  counts$nuclear_ld_pruned <- length(snps)

  pca <- snpgdsPCA(g, snp.id = snps, autosome.only = FALSE, num.thread = threads, verbose = FALSE)
  write_pca(pca, "pca_nuclear")

  king <- snpgdsIBDKING(g, snp.id = snps, autosome.only = FALSE, num.thread = threads, verbose = FALSE)
  pairs <- snpgdsIBDSelection(king)
  write_tsv(pairs[, c("ID1", "ID2", "IBS0", "kinship")], "kinship.tsv")
}
snpgdsClose(g)

# ── mitochondrial PCA (optional) ─────────────────────────────────────────────
# The mitochondrion is inherited as one unit, so there is no LD pruning: all
# polymorphic sites are used. Genotypes are coded 0/2 (hom ref / hom alt);
# heterozygous calls are kept here and reported separately as heteroplasmy.
mito_pca <- function() {
  snpgdsVCF2GDS(mito_vcf, "mito.gds", method = "biallelic.only", verbose = FALSE)
  m <- snpgdsOpen("mito.gds")
  on.exit(snpgdsClose(m))
  counts$mito_input <<- n_snps(m)
  if (counts$mito_input < 2) stop(sprintf("%d mito SNPs", counts$mito_input))
  pm <- snpgdsPCA(m, autosome.only = FALSE, remove.monosnp = TRUE,
                  num.thread = threads, verbose = FALSE)
  write_pca(pm, "pca_mito")
}
if (!grepl("\\.NONE$", mito_vcf)) {
  # Any failure (no mito SNPs, all monomorphic) is a note, not a failed run
  tryCatch(mito_pca(), error = function(e) {
    notes <<- c(notes, paste("Mito PCA skipped:", conditionMessage(e)))
  })
}

write_tsv(data.frame(set = names(counts), n_snps = unlist(counts)), "snp_counts.tsv")
writeLines(notes, "notes.txt")
