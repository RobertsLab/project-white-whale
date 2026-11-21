"""Tests for curated sample tables (src/whale/curation.py) and their use in the atlas. No network."""
from __future__ import annotations

from click.testing import CliRunner
from test_atlas import repo  # noqa: F401 - fixture: studies A and B on main

from whale.atlas import pool
from whale.cli import main
from whale.curation import COLUMNS, validate_curation

BASE = {"animal_id": "", "tissue": "gill", "life_stage": "adult", "age": "", "sex": "unknown",
        "ploidy": "diploid", "sample_unit": "individual", "pool_size": "", "population": "",
        "treatment": "control", "confidence": "high", "sources": "tissue: BioSample S1", "notes": ""}


def _write(path, rows):
    path.write_text("\t".join(COLUMNS) + "\n" + "".join(
        "\t".join(str(r.get(c, "")) for c in COLUMNS) + "\n" for r in rows))


def _rows(ids, **over):
    return [{**BASE, "sample_id": i, **over} for i in ids]


def test_valid_table_passes(tmp_path, repo_root):
    path = tmp_path / "sample_curation.tsv"
    _write(path, _rows(["L1", "L2"]))
    assert validate_curation(path, ["L1", "L2"], repo_root) == []


def test_libraries_must_match_the_sample_sheet(tmp_path, repo_root):
    path = tmp_path / "sample_curation.tsv"
    _write(path, _rows(["L1", "L3"]))
    (err,) = validate_curation(path, ["L1", "L2"], repo_root)
    assert "not in the sample sheet: L3" in err and "not curated: L2" in err
    _write(path, _rows(["L2", "L1"]))
    assert validate_curation(path, ["L1", "L2"], repo_root) == ["rows must be in the sample sheet's order"]


def test_terms_pools_and_sources_are_checked(tmp_path, repo_root):
    path = tmp_path / "sample_curation.tsv"
    _write(path, [
        {**BASE, "sample_id": "L1", "tissue": "gills", "life_stage": "gastrula", "sex": "f", "ploidy": "2n"},
        {**BASE, "sample_id": "L2", "sample_unit": "pool"},
        {**BASE, "sample_id": "L3", "pool_size": "1", "sources": " ", "confidence": "sure"},
    ])
    errors = "\n".join(validate_curation(path, ["L1", "L2", "L3"], repo_root))
    for text in ("tissue 'gills'", "life_stage 'gastrula'", "sex 'f'", "ploidy '2n'",
                 "sample_unit is pool but pool_size is empty", "pool_size '1'",
                 "L3): no sources", "confidence 'sure'"):
        assert text in errors, text


def test_unknown_is_accepted_everywhere(tmp_path, repo_root):
    path = tmp_path / "sample_curation.tsv"
    _write(path, _rows(["L1"], tissue="unknown", life_stage="unknown", ploidy="unknown", sample_unit="unknown"))
    assert validate_curation(path, ["L1"], repo_root) == []


def test_atlas_prefers_curated_attributes(repo):  # noqa: F811
    d = repo / "data" / "studies" / "A2020_X_WGBS"
    _write(d / "sample_curation.tsv", [
        {**BASE, "sample_id": "A1", "animal_id": "oyster1", "tissue": "gonad", "sex": "male", "population": "wild"},
        {**BASE, "sample_id": "A2", "animal_id": "oyster1", "tissue": "gonad", "sex": "male", "population": "wild"},
        {**BASE, "sample_id": "A3", "tissue": "unknown", "life_stage": "unknown", "sample_unit": "pool",
         "pool_size": "5", "sex": "unknown"},
        {**BASE, "sample_id": "A4", "tissue": "unknown", "life_stage": "unknown", "sample_unit": "pool",
         "pool_size": "5", "sex": "unknown"},
    ])
    p = pool(repo)
    meta = {m["library_id"]: m for m in p["libraries"]}
    a1, a3, b1 = meta["A1"], meta["A3"], meta["B1"]
    assert (a1["tissue"], a1["sex"], a1["animal_id"], a1["population"]) == ("gonad", "male", "oyster1", "wild")
    assert (a3["tissue"], a3["life_stage"], a3["sample_unit"], a3["pool_size"], a3["sex"]) == ("", "", "pool", "5", "unknown")
    assert a1["metadata_source"] == "curated" and b1["metadata_source"] == "deposit"
    assert a1["unit_id"] == "S1"                     # units still come from BioSamples
    assert [s["curated"] for s in p["studies"]] == [True, False]


def test_cli_validates_every_table(repo, monkeypatch):  # noqa: F811
    d = repo / "data" / "studies" / "A2020_X_WGBS"
    _write(d / "sample_curation.tsv", _rows(["A1", "A2", "A3"]))
    monkeypatch.chdir(repo)
    res = CliRunner().invoke(main, ["validate-curation"])
    assert res.exit_code == 1 and "not curated: A4" in res.output
    _write(d / "sample_curation.tsv", _rows(["A1", "A2", "A3", "A4"]))
    res = CliRunner().invoke(main, ["validate-curation"])
    assert res.exit_code == 0 and "ok" in res.output
