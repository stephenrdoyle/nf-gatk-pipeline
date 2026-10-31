#!/usr/bin/env python3
"""
Tests for bin/site_filter.py: the per-site decisions, the filter summary,
the diagnostics and the verify step, on small synthetic site tables laid out
as SITE_FILTER's bcftools query writes them.

Usage: python3 -m unittest discover -s tests -p 'test_*.py'
"""
import contextlib
import gzip
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
import site_filter as sf  # noqa: E402

TAGS = ["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "DP", "ExcessHet"]
GOOD = {"QD": 20, "FS": 1, "SOR": 0.7, "MQ": 60, "MQRankSum": 0, "ReadPosRankSum": 0, "DP": 1000, "ExcessHet": 3, "F_MISSING": 0.01}


def site(pos, ref="A", alt="G", **ann):
    v = dict(GOOD, **ann)
    return ["chr1", str(pos), ref, alt, str(v["F_MISSING"])] + [str(v[t]) for t in TAGS]


def read(path):
    with open(path) as fh:
        return fh.read()


class Run:
    def __init__(self, tmp, sites, **opts):
        self.tmp = tmp
        path = os.path.join(tmp, "sites.tsv.gz")
        with gzip.open(path, "wt") as fh:
            for s in sites:
                fh.write("\t".join(s) + "\n")
        contigs = os.path.join(tmp, "contigs.txt")
        with open(contigs, "w") as fh:
            fh.write("##contig=<ID=chr1,length=1000000>\n")
        p = lambda n: os.path.join(tmp, n)
        argv = ["apply", "--sites", path, "--columns", ",".join(TAGS), "--contigs", contigs,
                "--mode", opts.pop("mode", "hard"), "--summary", p("summary.tsv"), "--expected", p("expected.tsv"),
                "--diag-tsv", p("diag.tsv"), "--diag-html", p("diag.html")]
        for k, v in opts.items():
            argv += [f"--{k.replace('_', '-')}", str(v)]
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            sf.main(argv)
        self.vcf = out.getvalue()
        self.records = [l.split("\t") for l in self.vcf.splitlines() if not l.startswith("#")]
        self.filters = {int(r[1]): r[6] for r in self.records}
        lines = lambda n: read(p(n)).splitlines()
        self.summary = {line.split("\t")[0]: line.split("\t") for line in lines("summary.tsv")}
        diag = [l.split("\t") for l in lines("diag.tsv")]
        self.diag = [dict(zip(diag[0], row)) for row in diag[1:]]
        self.html = read(p("diag.html"))
        self.expected = p("expected.tsv")


class SiteFilterTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_classify(self):
        self.assertEqual(sf.classify("A", "G"), "snp")
        self.assertEqual(sf.classify("A", "G,T"), "snp")
        self.assertEqual(sf.classify("A", "G,*"), "snp")
        self.assertEqual(sf.classify("A", "AT"), "indel")
        self.assertEqual(sf.classify("A", "G,AT"), "indel")    # mixed goes with indels
        self.assertEqual(sf.classify("AC", "GT"), "indel")     # MNP: not a SNP

    def test_hard_rules_by_type(self):
        r = Run(self.tmp, [
            site(1),
            site(2, QD=1.5),
            site(3, ref="A", alt="AT", QD=1.5),
            site(4, ref="A", alt="AT", MQ=10),         # MQ is a SNP-only rule
            site(5, FS=61, SOR=3.5),
            site(6, ref="A", alt="AT", FS=150),        # indel FS limit is 200
            site(7, QD="."),                           # missing annotation passes
        ])
        self.assertEqual(r.filters[1], "PASS")
        self.assertEqual(r.filters[2], "SNP_QD2")
        self.assertEqual(r.filters[3], "INDEL_QD2")
        self.assertEqual(r.filters[4], "PASS")
        self.assertEqual(r.filters[5], "SNP_FS60;SNP_SOR3")
        self.assertEqual(r.filters[6], "PASS")
        self.assertEqual(r.filters[7], "PASS")

    def test_decisions_vcf_header(self):
        r = Run(self.tmp, [site(1), site(2, QD=1)])
        self.assertIn('##FILTER=<ID=SNP_QD2,Description="QD < 2">', r.vcf)
        self.assertIn("##contig=<ID=chr1,length=1000000>", r.vcf)
        self.assertEqual(len(r.records), 2)
        self.assertEqual(r.records[0][:5], ["chr1", "1", ".", "A", "G"])

    def test_optional_filters_off_by_default(self):
        r = Run(self.tmp, [site(1, DP=99999, F_MISSING=0.9, ExcessHet=500)] + [site(i) for i in range(2, 10)])
        self.assertEqual(r.filters[1], "PASS")
        self.assertNotIn("HighDepth", r.vcf)

    def test_max_depth_relative_per_type(self):
        sites = [site(i, DP=1000) for i in range(1, 10)] + [site(10, DP=2500)]
        sites += [site(20 + i, ref="A", alt="AT", DP=100) for i in range(9)] + [site(30, ref="A", alt="AT", DP=250)]
        r = Run(self.tmp, sites, max_depth="2x")
        self.assertEqual(r.filters[10], "HighDepth")   # > 2 x SNP median (1000)
        self.assertEqual(r.filters[30], "HighDepth")   # > 2 x indel median (100), not the SNP one
        self.assertEqual(r.filters[1], "PASS")

    def test_max_depth_absolute(self):
        r = Run(self.tmp, [site(1, DP=500), site(2, DP=1500)], max_depth="1000")
        self.assertEqual((r.filters[1], r.filters[2]), ("PASS", "HighDepth"))

    def test_missing_and_excess_het(self):
        r = Run(self.tmp, [site(1, F_MISSING=0.3), site(2, ExcessHet=60), site(3)], max_missing=0.2, excess_het="true")
        self.assertEqual(r.filters[1], "HighMissing")
        self.assertEqual(r.filters[2], "ExcessHet")
        self.assertEqual(r.filters[3], "PASS")

    def test_custom(self):
        r = Run(self.tmp, [site(1, QD=4), site(2, QD=6), site(3, ref="A", alt="AT", FS=100)],
                mode="custom", snp_filters="QD < 5.0, MQ<=30", indel_filters="INFO/FS > 50")
        self.assertEqual(r.filters[1], "SNP_QD_lt_5.0")
        self.assertEqual(r.filters[2], "PASS")
        self.assertEqual(r.filters[3], "INDEL_FS_gt_50")

    def test_custom_bad_expression(self):
        with self.assertRaises(SystemExit):
            sf.parse_expressions("QD < 2.0, vc.isSNP()", "snp")
        self.assertEqual(sf.custom_name("snp", "ReadPosRankSum", "<", "-8.0"), "SNP_ReadPosRankSum_lt_neg8.0")

    def test_summary(self):
        r = Run(self.tmp, [site(1), site(2, QD=1), site(3, QD=1, FS=70), site(4, ref="A", alt="AT")])
        self.assertEqual(r.summary["SNPs: all sites"][4], "3")
        self.assertEqual(r.summary["SNPs: PASS"][4], "1")
        self.assertEqual(r.summary["SNPs: SNP_QD2"][4:7], ["2", "66.67", "1"])   # 2 flagged, 1 by QD alone
        self.assertEqual(r.summary["SNPs: SNP_QD2"][3], "QD < 2")
        self.assertEqual(r.summary["Indels (incl. mixed): PASS"][4], "1")

    def test_diagnostics_warn_and_depth_rows(self):
        sites = [site(i, QD=1 if i <= 3 else 20) for i in range(1, 11)]   # QD < 2 flags 30%
        r = Run(self.tmp, sites)
        qd = next(d for d in r.diag if d["variant_type"] == "SNPs" and d["rule"] == "QD < 2")
        self.assertEqual(qd["pct_flagged"], "30.0")
        self.assertIn("flags 30.0%", qd["warning"])
        depth_rows = [d for d in r.diag if d["variant_type"] == "SNPs" and d["annotation"] == "DP"]
        self.assertEqual([d["applied"] for d in depth_rows], ["no"] * 4)   # shown, not applied
        self.assertIn("<svg", r.html)
        self.assertNotIn("<script", r.html)
        self.assertIn("QD &lt; 2", r.html)

    def test_verify(self):
        r = Run(self.tmp, [site(1), site(2, QD=1, FS=70)])
        good = os.path.join(self.tmp, "observed_good.txt")
        with open(good, "w") as fh:
            fh.write("PASS\nSNP_FS60;SNP_QD2\n")          # order within a value doesn't matter
        with contextlib.redirect_stderr(io.StringIO()):
            sf.main(["verify", r.expected, good])
        bad = os.path.join(self.tmp, "observed_bad.txt")
        with open(bad, "w") as fh:
            fh.write(".\n.\n")                           # annotate matched nothing
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            sf.main(["verify", r.expected, bad])


if __name__ == "__main__":
    unittest.main()
