"""Tests for workflows/methylation/bin/methylome_qc.py (the CPG_SUMMARY and QC_GATE steps)."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import sys

import pytest
from conftest import REPO

_spec = importlib.util.spec_from_file_location("methylome_qc", REPO / "workflows/methylation/bin/methylome_qc.py")
qc = importlib.util.module_from_spec(_spec)
sys.modules["methylome_qc"] = qc
_spec.loader.exec_module(qc)

# Trimmed from ROBERTS2022_PH_PLOIDY_WGBS SRR13207050 (klone, 2026-10-24).
PE_REPORT = """Final Alignment report
======================
Sequence pairs analysed in total:\t51972926
Number of paired-end alignments with a unique best hit:\t15706901
Mapping efficiency:\t30.2% 
Sequence pairs with no alignments under any condition:\t32523284

C methylated in CpG context:\t11.4%
"""
SE_REPORT = "Sequences analysed in total:\t1000\nMapping efficiency:\t97.5%\n"
DEDUP_REPORT = """Total number of alignments analysed in S.bam:\t15706897
Total number duplicated alignments removed:\t1677553 (10.68%)
Total count of deduplicated leftover sequences: 14029344 (89.32% of total)
"""
SPLITTING = """Total methylated C's in CpG context:\t5433044
Total methylated C's in CHG context:\t599585
Total methylated C's in CHH context:\t4114184

Total C to T conversions in CpG context:\t43475227
Total C to T conversions in CHG context:\t67070165
Total C to T conversions in CHH context:\t267870358
"""


def test_parse_reports():
    assert qc.parse_alignment_report(PE_REPORT) == {"pairs_or_reads": 51972926, "mapping_pct": 30.2}
    assert qc.parse_alignment_report(SE_REPORT) == {"pairs_or_reads": 1000, "mapping_pct": 97.5}
    assert qc.parse_dedup_report(DEDUP_REPORT) == {"duplication_pct": 10.68}
    s = qc.parse_splitting_report(SPLITTING)
    assert s["cpg_meth_pct"] == pytest.approx(11.109, abs=0.001)
    assert s["non_cpg_meth_pct"] == pytest.approx(1.388, abs=0.001)
    assert s["conversion_pct"] == pytest.approx(100 - s["non_cpg_meth_pct"])
    assert qc.parse_splitting_report("")["conversion_pct"] is None


def _cov(path, rows):
    with gzip.open(path, "wt") as fh:
        for i, (m, u) in enumerate(rows):
            fh.write(f"chr1\t{i * 10}\t{i * 10 + 2}\t0\t{m}\t{u}\n")
    return path


def test_cpg_summary(tmp_path):
    cov = _cov(tmp_path / "S1.cov.gz", [(0, 0), (1, 0), (4, 1), (2, 8), (10, 0)])
    s = qc.cpg_summary("S1", cov)
    assert (s["cpgs_1x"], s["cpgs_5x"], s["cpgs_10x"]) == (4, 3, 2)
    assert s["mean_depth"] == pytest.approx(26 / 4)
    assert s["cpg_meth_pooled_pct"] == pytest.approx(100 * 17 / 26, abs=0.001)
    assert s["cpg_meth_mean_5x_pct"] == pytest.approx(100 * (0.8 + 0.2 + 1.0) / 3, abs=0.001)


def test_qc_table_gate(tmp_path):
    sheet = tmp_path / "samplesheet.tsv"
    sheet.write_text("sample_id\tcondition\nS1\tA\nS10\tB\nS2\tB\n")
    (tmp_path / "S1_PE_report.txt").write_text(PE_REPORT)
    (tmp_path / "S1.dedup.deduplication_report.txt").write_text(DEDUP_REPORT)
    (tmp_path / "S1.dedup_splitting_report.txt").write_text(SPLITTING)
    qc.main(["cpg-summary", "S1", str(_cov(tmp_path / "c1.gz", [(5, 5)] * 3)), str(tmp_path / "S1.cpg_summary.tsv")])
    # S10 maps too little; S2 has no reports at all.
    (tmp_path / "S10_PE_report.txt").write_text(PE_REPORT.replace("30.2%", "12.0%"))
    (tmp_path / "S10.dedup_splitting_report.txt").write_text(SPLITTING)
    qc.main(["cpg-summary", "S10", str(_cov(tmp_path / "c10.gz", [(5, 5)] * 3)), str(tmp_path / "S10.cpg_summary.tsv")])
    out = tmp_path / "qc_summary.tsv"
    qc.main(["qc-table", "--samplesheet", str(sheet), "--out", str(out), "--min-cpgs-5x", "2",
             *map(str, sorted(tmp_path.glob("S*")))])
    with open(out) as fh:
        rows = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert list(rows) == ["S1", "S10", "S2"]
    assert rows["S1"]["qc_pass"] == "true" and rows["S1"]["qc_reasons"] == ""
    assert rows["S1"]["mapping_pct"] == "30.2" and rows["S1"]["cpgs_5x"] == "3" and rows["S1"]["condition"] == "A"
    assert rows["S10"]["qc_pass"] == "false" and rows["S10"]["qc_reasons"].startswith("mapping 12.0%")
    assert rows["S2"]["qc_pass"] == "false" and "no alignment report" in rows["S2"]["qc_reasons"]
