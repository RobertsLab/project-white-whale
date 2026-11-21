"""Tests for workflows/rnaseq/bin/rnaseq_quant.py (the GENE_MATRICES and RNASEQ_QC steps)."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys

import pytest
from conftest import REPO

_spec = importlib.util.spec_from_file_location("rnaseq_quant", REPO / "workflows/rnaseq/bin/rnaseq_quant.py")
rq = importlib.util.module_from_spec(_spec)
sys.modules["rnaseq_quant"] = rq
_spec.loader.exec_module(rq)

TX2GENE = {"T1a": "10", "T1b": "10", "T2": "9", "T3": "G3"}


def quant_sf(path, rows):
    path.mkdir(parents=True, exist_ok=True)
    with open(path / "quant.sf", "w") as fh:
        fh.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for name, length, eff, tpm, n in rows:
            fh.write(f"{name}\t{length}\t{eff}\t{tpm}\t{n}\n")


def matrix(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_aggregate_sums_tpm_and_counts_and_weights_length_like_tximport():
    q = [{"Name": "T1a", "Length": 1000, "EffectiveLength": 800.0, "TPM": 30.0, "NumReads": 24.0},
         {"Name": "T1b", "Length": 2000, "EffectiveLength": 1800.0, "TPM": 10.0, "NumReads": 18.0},
         {"Name": "T2", "Length": 500, "EffectiveLength": 300.0, "TPM": 0.0, "NumReads": 0.0},
         {"Name": "Tx", "Length": 500, "EffectiveLength": 300.0, "TPM": 5.0, "NumReads": 1.0}]
    out = rq.aggregate(q, TX2GENE)
    assert out["10"]["tpm"] == 40.0 and out["10"]["counts"] == 42.0
    assert out["10"]["length"] == pytest.approx((30 * 800 + 10 * 1800) / 40)
    # No expression: the plain mean of the effective lengths, as tximport's fallback.
    assert out["9"] == {"tpm": 0.0, "counts": 0.0, "length": 300.0}
    assert "Tx" not in out and set(out) == {"10", "9"}


def test_build_matrices_writes_gene_and_transcript_tables_in_sheet_order(tmp_path):
    quant_sf(tmp_path / "salmon" / "B", [("T1a", 1000, 800, 30, 24), ("T1b", 2000, 1800, 10, 18),
                                          ("T2", 500, 300, 60, 18), ("T3", 400, 200, 0, 0), ("Tz", 300, 100, 0, 0)])
    quant_sf(tmp_path / "salmon" / "A", [("T1a", 1000, 800, 0, 0), ("T1b", 2000, 1800, 50, 90),
                                          ("T2", 500, 300, 50, 15), ("T3", 400, 200, 0, 0), ("Tz", 300, 100, 0, 0)])
    dirs = {p.name: p for p in (tmp_path / "salmon").iterdir()}
    report = rq.build_matrices(["A", "B"], dirs, TX2GENE, tmp_path / "m")
    assert report["transcripts_quantified"] == 5 and report["transcripts_not_in_tx2gene"] == 1
    assert report["genes"] == 3 and report["tx2gene_genes_without_quantified_transcript"] == 0

    counts = matrix(tmp_path / "m" / "gene_counts.tsv.gz")
    # Numeric GeneIDs sort as numbers, then the rest.
    assert counts == [["gene_id", "A", "B"], ["9", "15", "18"], ["10", "90", "42"], ["G3", "0", "0"]]
    tpm = matrix(tmp_path / "m" / "gene_tpm.tsv.gz")
    assert tpm[2] == ["10", "50", "40"]
    length = matrix(tmp_path / "m" / "gene_length.tsv.gz")
    assert length[2] == ["10", "1800", "1050"] and length[3] == ["G3", "200", "200"]
    tx = matrix(tmp_path / "m" / "transcript_counts.tsv.gz")
    assert [r[0] for r in tx] == ["transcript_id", "T1a", "T1b", "T2", "T3", "Tz"]
    # Byte-identical on a re-run, so recording it twice commits nothing.
    before = (tmp_path / "m" / "gene_tpm.tsv.gz").read_bytes()
    rq.build_matrices(["A", "B"], dirs, TX2GENE, tmp_path / "m")
    assert (tmp_path / "m" / "gene_tpm.tsv.gz").read_bytes() == before


def test_build_matrices_refuses_mixed_indexes_and_foreign_tx2gene(tmp_path):
    quant_sf(tmp_path / "A", [("T1a", 1, 1, 1, 1)])
    quant_sf(tmp_path / "B", [("T2", 1, 1, 1, 1)])
    dirs = {"A": tmp_path / "A", "B": tmp_path / "B"}
    with pytest.raises(rq.QuantError, match="other transcripts"):
        rq.build_matrices(["A", "B"], dirs, TX2GENE, tmp_path / "m")
    with pytest.raises(rq.QuantError, match="no transcript in quant.sf"):
        rq.build_matrices(["A"], dirs, {"X": "1"}, tmp_path / "m")
    with pytest.raises(rq.QuantError, match="no Salmon output"):
        rq.build_matrices(["A", "C"], dirs, TX2GENE, tmp_path / "m")


def test_read_tx2gene_skips_a_header_and_refuses_two_genes(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("transcript_id\tgene_id\nT1\t1\nT2\t1\n")
    assert rq.read_tx2gene(p) == {"T1": "1", "T2": "1"}
    p.write_text("T1\t1\nT1\t2\n")
    with pytest.raises(rq.QuantError, match="both"):
        rq.read_tx2gene(p)


def test_parse_fastp_counts_pairs_for_paired_end():
    pe = {"summary": {"before_filtering": {"total_reads": 2000}, "after_filtering": {"total_reads": 1900, "q30_rate": 0.93}},
          "read2_before_filtering": {}}
    assert rq.parse_fastp(pe) == {"layout": "PE", "reads_in": 1000, "reads_after_trim_pct": 95.0, "q30_pct": 93.0}
    se = {"summary": {"before_filtering": {"total_reads": 500}, "after_filtering": {"total_reads": 500}}}
    assert rq.parse_fastp(se)["layout"] == "SE" and rq.parse_fastp(se)["reads_in"] == 500


@pytest.mark.parametrize("counts, call", [
    ({"ISR": 950, "ISF": 50}, ("reverse", 5.0)),
    ({"SF": 900, "SR": 100}, ("forward", 90.0)),
    ({"ISR": 500, "ISF": 480}, ("unstranded", 48.98)),
    ({"IU": 1000}, (None, None)),
])
def test_strandedness_comes_from_the_format_counts(counts, call):
    got = rq.strandedness(counts)
    assert got[0] == call[0] and (got[1] == pytest.approx(call[1], abs=0.01) if call[1] is not None else got[1] is None)


def _sample(tmp_path, sid, *, mapped, processed, decoys, lib):
    d = tmp_path / "reports" / sid
    (d / "aux_info").mkdir(parents=True)
    (d / "aux_info" / "meta_info.json").write_text(json.dumps(
        {"num_processed": processed, "num_mapped": mapped, "num_decoy_fragments": decoys,
         "percent_mapped": 100.0 * mapped / processed, "library_types": ["IU"]}))
    (d / "lib_format_counts.json").write_text(json.dumps(lib))
    f = tmp_path / "reports" / f"{sid}.fastp.json"
    f.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": processed},
                                         "after_filtering": {"total_reads": processed}}}))
    return [d, f]


def test_qc_table_applies_the_gate_and_keeps_failures(tmp_path):
    files = _sample(tmp_path, "A", mapped=800, processed=1000, decoys=50,
                    lib={"expected_format": "SR", "compatible_fragment_ratio": 0.99, "SR": 790, "SF": 10})
    files += _sample(tmp_path, "B", mapped=100, processed=1000, decoys=5,
                     lib={"expected_format": "U", "compatible_fragment_ratio": 1.0, "SR": 50, "SF": 50})
    m = tmp_path / "m"
    m.mkdir()
    with gzip.open(m / "gene_tpm.tsv.gz", "wt") as fh:
        fh.write("gene_id\tA\tB\ng1\t5\t0.5\ng2\t3\t2\n")
    sheet = [{"sample_id": "A", "condition": "c"}, {"sample_id": "B", "condition": "t"}]
    rows = rq.qc_table(sheet, files, m, min_mapping=30, min_mapped=500, min_genes=2)
    a, b = rows
    assert (a["qc_pass"], a["mapping_pct"], a["decoy_pct"], a["strandedness"], a["genes_tpm_ge1"]) == \
        ("true", 80.0, 5.0, "reverse", 2)
    assert b["qc_pass"] == "false" and b["strandedness"] == "unstranded"
    assert b["qc_reasons"] == "mapping 10.0% < 30%; mapped fragments 100 < 500; genes with TPM >= 1 1 < 2"
    rq.write_table(rows, tmp_path / "qc.tsv")
    with open(tmp_path / "qc.tsv") as fh:
        assert next(csv.reader(fh, delimiter="\t")) == rq.QC_COLUMNS
