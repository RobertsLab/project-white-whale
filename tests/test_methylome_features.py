"""Tests for workflows/methylation/bin/methylome_features.py (the FEATURE_SUMMARY step)."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import sys

import pytest
from conftest import REPO

_spec = importlib.util.spec_from_file_location("methylome_features",
                                               REPO / "workflows/methylation/bin/methylome_features.py")
mf = importlib.util.module_from_spec(_spec)
sys.modules["methylome_features"] = mf
_spec.loader.exec_module(mf)

FIXTURE_FEATURES = REPO / "tests/fixtures/reference/TEST/features"


@pytest.fixture
def features(tmp_path):
    d = tmp_path / "features"
    d.mkdir()
    # GA spans 100-500 with two overlapping exons; GB overlaps GA's end; GC has no CpG.
    (d / "genes.bed").write_text("c1\t100\t500\tGA\t.\t+\nc1\t400\t700\tGB\t.\t-\nc2\t0\t50\tGC\t.\t+\n")
    (d / "promoters.bed").write_text("c1\t0\t120\tGA\t.\t+\nc1\t690\t900\tGB\t.\t-\n")
    (d / "exons.bed").write_text("c1\t100\t200\tGA\t.\t+\nc1\t150\t250\tGA\t.\t+\nc1\t400\t700\tGB\t.\t-\n")
    (d / "introns.bed").write_text("c1\t250\t500\tGA\t.\t+\n")
    return d


def _cov(path, rows):
    with gzip.open(path, "wt") as fh:
        for chrom, pos, m, u in rows:
            fh.write(f"{chrom}\t{pos}\t{pos + 2}\t0\t{m}\t{u}\n")
    return path


def _read(path):
    with gzip.open(path, "rt") as fh:
        return {(r["feature"], r["gene_id"]): r for r in csv.DictReader(fh, delimiter="\t")}


def test_feature_summary(tmp_path, features):
    cov = _cov(tmp_path / "s.cov.gz", [
        ("c1", 50, 1, 1),      # GA promoter only
        ("c1", 160, 8, 2),     # GA body; in both GA exons, counted once
        ("c1", 300, 0, 3),     # GA body, intron
        ("c1", 450, 5, 0),     # GA body + intron, and GB body + exon
        ("c1", 800, 2, 0),     # GB promoter
        ("c2", 10, 0, 0),      # no calls: ignored
    ])
    out = tmp_path / "s.features.tsv.gz"
    assert mf.main(["s", str(cov), str(features), str(out)]) == 0
    r = _read(out)
    assert len(r) == 4 * 3                       # every gene, every class, covered or not
    ga = r[("gene_body", "GA")]
    assert (ga["cpgs"], ga["cpgs_5x"], ga["meth"], ga["unmeth"]) == ("3", "2", "13", "5")
    assert float(ga["meth_pct"]) == pytest.approx(100 * 13 / 18, abs=1e-3)
    assert float(ga["mean_meth_5x_pct"]) == pytest.approx(100 * (0.8 + 1.0) / 2, abs=1e-3)
    assert r[("exon", "GA")]["cpgs"] == "1" and r[("intron", "GA")]["cpgs"] == "2"
    assert r[("promoter", "GA")]["meth"] == "1" and r[("promoter", "GB")]["meth"] == "2"
    assert r[("gene_body", "GB")]["cpgs"] == "1" and r[("exon", "GB")]["cpgs"] == "1"
    gc = r[("gene_body", "GC")]
    assert (gc["cpgs"], gc["meth_pct"], gc["mean_meth_5x_pct"]) == ("0", "", "")
    # Byte-identical on a re-run.
    again = tmp_path / "again.tsv.gz"
    mf.main(["s", str(cov), str(features), str(again)])
    assert again.read_bytes() == out.read_bytes()


def test_snp_mask_leaves_cpgs_out(tmp_path, features):
    cov = _cov(tmp_path / "s.cov.gz", [("c1", 160, 8, 2), ("c1", 300, 0, 3), ("c1", 450, 5, 0)])
    mask = tmp_path / "cpg_mask.bed.gz"
    with gzip.open(mask, "wt") as fh:
        fh.write("c1\t300\t302\t301C>T\t9\t.\n")       # the unmethylated-looking CpG
    out = tmp_path / "m.tsv.gz"
    assert mf.main(["s", str(cov), str(features), str(out), str(mask)]) == 0
    ga = _read(out)[("gene_body", "GA")]
    assert (ga["cpgs"], ga["meth"], ga["unmeth"]) == ("2", "13", "2")
    stats = {}
    mf.summarize(cov, features, mf.read_mask(mask), stats)
    assert stats == {"covered": 3, "masked": 1}


def test_unsorted_input_is_refused(tmp_path, features):
    cov = _cov(tmp_path / "u.cov.gz", [("c1", 300, 1, 0), ("c1", 160, 1, 0)])
    with pytest.raises(ValueError, match="sorted"):
        mf.summarize(cov, features)
    cov = _cov(tmp_path / "v.cov.gz", [("c1", 1, 1, 0), ("c2", 1, 1, 0), ("c1", 5, 1, 0)])
    with pytest.raises(ValueError, match="contiguous"):
        mf.summarize(cov, features)


def test_fixture_features_cover_the_fixture_genome():
    """The workflow's test profile reads these; every gene sits inside its chromosome."""
    sizes = {"chr1": 20_000, "chr2": 10_000}
    genes = mf.gene_universe(FIXTURE_FEATURES)
    assert genes == ["G1", "G2", "G3", "G4"]
    for name in ("genes.bed", "promoters.bed", "exons.bed", "introns.bed"):
        for chrom, rows in mf.read_bed(FIXTURE_FEATURES / name).items():
            assert all(0 <= s < e <= sizes[chrom] and g in genes for s, e, g in rows), name
