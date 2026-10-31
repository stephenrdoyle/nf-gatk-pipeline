#!/usr/bin/env python3
"""
Tests for bin/qc_report.py on small synthetic inputs, laid out as QC_REPORT
stages them (inputs/ holding bcftools stats, SNPRelate tables, ...).

Each planted problem (a duplicate pair, a sex mismatch, a heteroplasmic
mito sample, a missing joint-genotyping piece, a missingness outlier) must
be flagged, and a clean cohort must flag nothing.

Usage: python3 -m unittest discover -s tests -p 'test_*.py'
"""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
import qc_report  # noqa: E402

SAMPLES = [f"S{i:02d}" for i in range(1, 11)]

# bcftools 1.19 PSC layout (columns are read by name, not position)
PSC_HEADER = ("# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions"
              "\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons"
              "\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing")


def stats_file(path, per_sample, records=1000, tstv=2.0):
    """per_sample: {sample: (refhom, nonrefhom, hets, depth, singletons)}"""
    with open(path, "w") as f:
        f.write(f"SN\t0\tnumber of records:\t{records}\n")
        f.write(f"SN\t0\tnumber of SNPs:\t{records - 100}\n")
        f.write(f"SN\t0\tnumber of indels:\t100\n")
        f.write(f"TSTV\t0\t600\t300\t{tstv}\t600\t300\t{tstv}\n")
        f.write(PSC_HEADER + "\n")
        for s, (rr, aa, het, dp, sing) in per_sample.items():
            f.write(f"PSC\t0\t{s}\t{rr}\t{aa}\t{het}\t0\t0\t0\t{dp}\t{sing}\t0\t0\t0\n")


def tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


class Cohort:
    """Writes a clean 10-sample cohort; tests then plant problems in it."""

    def __init__(self, tmp, sex_chrom="chrX", mito_chrom="mtDNA"):
        self.tmp, self.inputs = tmp, os.path.join(tmp, "inputs")
        os.makedirs(self.inputs)
        self.nuclear = {s: (500, 200, 300, 20 + i % 3, 10) for i, s in enumerate(SAMPLES)}
        # half male (no chrX hets), half female
        self.sex = {s: ("M" if i % 2 else "F") for i, s in enumerate(SAMPLES)}
        self.x = {s: (50, 20, 0 if self.sex[s] == "M" else 30, 10, 1) for s in SAMPLES}
        self.mito = {s: (10, 5, 0, 400, 0) for s in SAMPLES}
        self.missing = {s: 2.0 + (i % 3) * 0.1 for i, s in enumerate(SAMPLES)}
        self.kin = [(a, b, 0.0) for i, a in enumerate(SAMPLES) for b in SAMPLES[i + 1:]]
        self.meta = {"run_name": "test_run", "sex_chrom": sex_chrom, "mito_chrom": mito_chrom,
                     "qc_maf": 0.05, "qc_ld_threshold": 0.2, "n_samples_expected": len(SAMPLES),
                     "n_pieces_expected": 12, "n_pieces_gathered": 12}

    def build(self):
        p = lambda n: os.path.join(self.inputs, n)
        stats_file(p("stats_nuclear.txt"), self.nuclear)
        stats_file(p("before_filter_bcftools_stats.txt"), self.nuclear, records=2000, tstv=1.8)
        stats_file(p("after_filter_bcftools_stats.txt"), self.nuclear, records=1500, tstv=2.1)
        if self.meta["sex_chrom"]:
            stats_file(p("stats_sex.txt"), self.x)
        if self.meta["mito_chrom"]:
            stats_file(p("stats_mito.txt"), self.mito)
            tsv(p("pca_mito.tsv"), ["sample", "PC1", "PC2"], [(s, i * 0.1, -i * 0.05) for i, s in enumerate(SAMPLES)])
            tsv(p("pca_mito_eigen.tsv"), ["pc", "variance_explained"], [("PC1", 0.6), ("PC2", 0.2)])
        tsv(p("after_filter_missingness.tsv"), ["sample", "pct_missing"], self.missing.items())
        tsv(p("pca_nuclear.tsv"), ["sample", "PC1", "PC2", "PC3", "PC4"],
            [(s, (i % 5) * 0.1, (i % 3) * 0.1, (i % 4) * 0.1, (i % 2) * 0.1) for i, s in enumerate(SAMPLES)])
        tsv(p("pca_nuclear_eigen.tsv"), ["pc", "variance_explained"],
            [("PC1", 0.1), ("PC2", 0.05), ("PC3", 0.04), ("PC4", 0.03)])
        tsv(p("kinship.tsv"), ["ID1", "ID2", "IBS0", "kinship"], [(a, b, 0.05, k) for a, b, k in self.kin])
        tsv(p("snp_counts.tsv"), ["set", "n_snps"], [("nuclear_input", 5000), ("nuclear_ld_pruned", 2000), ("mito_input", 15)])
        tsv(p("filter_summary.tsv"), ["id", "variant_type", "filter", "expression", "n_sites", "pct_sites", "n_only"],
            [("SNPs: all sites", "SNPs", "(all sites)", "", 1800, "100.00", ""),
             ("SNPs: PASS", "SNPs", "PASS", "passes every filter", 1400, "77.78", ""),
             ("SNPs: QD", "SNPs", "QD", "QD < 2.0", 200, "11.11", 150)])
        tsv(p("density.tsv"), ["contig", "mb", "n"], [("chr1", 0, 400), ("chr1", 1, 380), ("chrX", 0, 120)])
        open(p("notes.txt"), "w").close()
        with open(p("software_versions.txt"), "w") as f:
            f.write("pipeline: nf-gatk-pipeline\ngatk: 4.6.1.0\n")
        samples = os.path.join(self.tmp, "samples.tsv")
        tsv(samples, ["sample", "sex", "group"],
            [(s, self.sex[s], "farmA" if i < 5 else "farmB") for i, s in enumerate(SAMPLES)])
        return qc_report.build(self.inputs, samples, self.meta)


class QcReportTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.c = Cohort(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def checks(self, html):
        """The 'Before the results' list only."""
        start = html.index("<ul class='checks'>")
        return html[start:html.index("</ul>", start)]

    def test_clean_cohort_flags_nothing(self):
        html = self.c.build()
        checks = self.checks(html)
        self.assertNotIn('class="bad"', checks)
        self.assertIn("0 sample(s) flagged", checks)
        self.assertIn("10 in the cohort VCF, 10 in the samplesheet", checks)

    def test_self_contained(self):
        html = self.c.build()
        self.assertNotIn("<script", html)
        self.assertNotIn("http://", html)
        self.assertNotIn("https://", html)
        self.assertIn("<svg", html)

    def test_duplicate_pair(self):
        self.c.kin[0] = ("S01", "S02", 0.49)
        html = self.c.build()
        self.assertIn('class="bad"><span class="mark">✗</span> <b>Duplicates:</b> 1 sample pair', html)
        self.assertIn("duplicate of S02", html)

    def test_sex_mismatch(self):
        self.c.sex["S01"] = "M"   # labelled male, but chrX heterozygous like a female
        html = self.c.build()
        self.assertIn("1 mismatch(es)", self.checks(html))
        self.assertIn("sex mismatch (labelled M, chrX says F)", html)

    def test_mito_heteroplasmy(self):
        self.c.mito["S03"] = (10, 2, 6, 400, 0)
        html = self.c.build()
        self.assertIn("1 sample(s) with heterozygous mito calls", self.checks(html))

    def test_missing_piece_and_sample(self):
        self.c.meta["n_pieces_gathered"] = 11
        self.c.meta["n_samples_expected"] = 11
        checks = self.checks(self.c.build())
        self.assertIn('class="bad"><span class="mark">✗</span> <b>Genome:</b> 11 of 12', checks)
        self.assertIn('class="bad"><span class="mark">✗</span> <b>Samples:</b> 10 in the cohort VCF, 11', checks)

    def test_missingness_outlier(self):
        self.c.missing["S05"] = 40.0
        html = self.c.build()
        self.assertIn("high missingness", html)

    def test_without_sex_and_mito(self):
        c = Cohort(os.path.join(self._tmp.name, "plain"), sex_chrom="", mito_chrom="")
        os.makedirs(c.tmp, exist_ok=True)
        html = c.build()
        self.assertNotIn("Sex chromosome", html)
        self.assertNotIn("Mitochondria (", html)
        self.assertNotIn("<b>Sex:</b>", html)

    def test_empty_snprelate_output(self):
        """Too few SNPs: SNPRelate writes header-only tables and a note."""
        p = lambda n: os.path.join(self.c.inputs, n)
        self.c.build()
        tsv(p("pca_nuclear.tsv"), ["sample"], [])
        tsv(p("kinship.tsv"), ["ID1", "ID2", "IBS0", "kinship"], [])
        with open(p("notes.txt"), "w") as f:
            f.write("Nuclear PCA and kinship skipped: only 12 SNPs reached the QC set (100 needed).\n")
        samples = os.path.join(self.c.tmp, "samples.tsv")
        html = qc_report.build(self.c.inputs, samples, self.c.meta)
        self.assertIn("only 12 SNPs", html)
        self.assertIn("No PCA output", html)

    def test_cli_writes_file(self):
        self.c.build()
        meta = os.path.join(self.c.tmp, "meta.json")
        with open(meta, "w") as f:
            json.dump(self.c.meta, f)
        out = os.path.join(self.c.tmp, "qc_report.html")
        rc = qc_report.main(["--inputs", self.c.inputs, "--samples", os.path.join(self.c.tmp, "samples.tsv"),
                             "--meta", meta, "--out", out])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.getsize(out) > 1000)

    def test_robust_outliers(self):
        vals = {f"s{i}": 10.0 + (i % 3) * 0.1 for i in range(20)}
        self.assertEqual(qc_report.robust_outliers(vals), {})
        vals["s0"] = 30.0
        self.assertEqual(list(qc_report.robust_outliers(vals)), ["s0"])
        self.assertEqual(qc_report.robust_outliers(vals, side="low"), {})

    def test_king_degree(self):
        self.assertEqual(qc_report.king_degree(0.5), "duplicate")
        self.assertEqual(qc_report.king_degree(0.25), "first degree")
        self.assertEqual(qc_report.king_degree(0.12), "second degree")
        self.assertEqual(qc_report.king_degree(0.06), "third degree")
        self.assertEqual(qc_report.king_degree(0.01), "unrelated")


if __name__ == "__main__":
    unittest.main()
