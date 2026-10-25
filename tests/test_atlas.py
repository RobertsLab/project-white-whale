"""Tests for `whale build-atlas` (src/whale/atlas.py). No network access."""
from __future__ import annotations

import copy
import gzip
import json
import shutil
import subprocess

import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from whale.atlas import build, pool, spearman, units, write
from whale.cli import main

URL = "https://github.com/RobertsLab/project-white-whale"
MASK = {"mask_id": "GCF_X.snpmask.v1", "sha256": "0" * 64, "cpgs_masked": 1}
GENES = ["g1", "g2", "g3", "g4"]


def _gz(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def _study(repo, filled, sid, acc, libs, sheet, meth, calls, cpgs, *, flags=(), tissue="gill",
           mask=True, assembly="GCF_X", factors=None, by_factor=None):
    """Register a study on main with gene tables. meth/calls/cpgs: genes x libs."""
    s = copy.deepcopy(filled)
    s["study_id"], s["accessions"]["bioproject"] = sid, acc
    s["samples"]["tissues"] = [tissue]
    s["quality_flags"] = list(flags)
    if factors is not None:
        s["factors"] = factors
    (repo / "registry" / "studies" / f"{sid}.yaml").write_text(yaml.safe_dump(s))
    d = repo / "data" / "studies" / sid
    d.mkdir(parents=True)
    (d / "study_summary.yaml").write_text(yaml.safe_dump({
        "study_id": sid, "part": None, "genome_assembly": assembly, "chemistry": "wgbs",
        "snp_mask": MASK if mask else None, "nextflow_run": "run_x", "by_factor": by_factor or []}))
    cols = list(sheet[0])
    (d / "samplesheet.tsv").write_text("\t".join(cols) + "\n" + "".join(
        "\t".join(str(r[c]) for c in cols) + "\n" for r in sheet))
    (d / "qc_summary.tsv").write_text("sample_id\tcpg_meth_pct\tqc_pass\n" + "".join(
        f"{lib}\t{10 + k}\ttrue\n" for k, lib in enumerate(libs)))
    for cls in ("gene_body", "promoter"):
        for metric, m in (("meth", meth), ("calls", calls), ("cpgs", cpgs)):
            _gz(d / "genes" / f"{cls}.{metric}.tsv.gz",
                [["gene_id", *libs], *[[g, *row] for g, row in zip(GENES, m, strict=True)]])


@pytest.fixture
def repo(tmp_path, repo_root, filled_study):
    """Three studies on main: A (two animals, two technical-replicate libraries each),
    B (pools sharing one BioSample), C (not SNP-masked)."""
    for rel in ("registry", "schemas"):
        shutil.copytree(repo_root / rel, tmp_path / rel)
    for path in (tmp_path / "registry" / "studies").glob("[!_]*.yaml"):
        path.unlink()
    a_libs = ["A1", "A2", "A3", "A4"]
    _study(tmp_path, filled_study, "A2020_X_WGBS", "PRJNA1", a_libs,
           [{"sample_id": lib, "sample_accession": acc, "condition": cond, "pop": cond, "gender": sex}
            for lib, acc, cond, sex in (("A1", "S1", "w", "female"), ("A2", "S1", "w", "female"),
                                        ("A3", "S2", "i", ""), ("A4", "S2", "i", ""))],
           meth=[[10, 10, 0, 0], [5, 5, 50, 50], [80, 80, 90, 90], [0, 0, 0, 0]],
           calls=[[100] * 4, [100] * 4, [100] * 4, [0] * 4],
           cpgs=[[9, 7, 9, 9], [9] * 4, [9] * 4, [0] * 4],
           tissue="mantle",
           factors=[{"factor_id": "origin", "factor_original": "o", "factor_standardized": "genotype_or_population",
                     "levels": ["w", "i"], "reference_level": "w"}],
           by_factor=[{"factor_id": "origin", "column": "pop", "levels": []}])
    _study(tmp_path, filled_study, "B2021_Y_WGBS", "PRJNA2", ["B1", "B2"],
           [{"sample_id": b, "sample_accession": "SHARED", "condition": "c"} for b in ("B1", "B2")],
           meth=[[20, 30], [10, 0], [60, 70], [0, 0]], calls=[[100, 100], [100, 100], [100, 100], [0, 0]],
           cpgs=[[9, 9], [9, 9], [9, 9], [0, 0]], flags=["shared_biosample"])
    _study(tmp_path, filled_study, "C2022_Z_WGBS", "PRJNA3", ["C1"], [{"sample_id": "C1", "condition": "c"}],
           meth=[[1], [1], [1], [1]], calls=[[9], [9], [9], [9]], cpgs=[[9], [9], [9], [9]], mask=False)
    return tmp_path


def test_pool_takes_masked_merged_studies(repo):
    p = pool(repo)
    assert [s["study_id"] for s in p["studies"]] == ["A2020_X_WGBS", "B2021_Y_WGBS"]
    assert p["excluded"] == [{"study_id": "C2022_Z_WGBS", "part": None, "reason": "gene tables not SNP-masked"}]
    assert p["matrices"][("gene_body", "meth")].shape == (4, 6)
    meta = {m["library_id"]: m for m in p["libraries"]}
    assert (meta["A1"]["unit_id"], meta["A3"]["unit_id"]) == ("S1", "S2")      # replicates share a unit
    assert (meta["B1"]["unit_id"], meta["B2"]["unit_id"]) == ("B1", "B2")      # shared BioSample: per library
    assert (meta["A1"]["tissue"], meta["A1"]["sex"], meta["A3"]["sex"]) == ("mantle", "female", "unknown")
    assert meta["A1"]["factors"] == "origin=w"                                   # via by_factor's column
    assert meta["B1"]["snp_mask"] == "GCF_X.snpmask.v1"
    assert [s["units"] for s in p["studies"]] == [2, 2]


def test_other_assembly_is_excluded(repo, filled_study):
    _study(repo, filled_study, "D2023_W_WGBS", "PRJNA4", ["D1"], [{"sample_id": "D1", "condition": "c"}],
           meth=[[1]] * 4, calls=[[9]] * 4, cpgs=[[9]] * 4, assembly="GCF_OTHER")
    reasons = {e["study_id"]: e["reason"] for e in pool(repo)["excluded"]}
    assert reasons["D2023_W_WGBS"] == "on GCF_OTHER, not GCF_X"


def test_units_sum_calls_and_take_max_cpgs(repo):
    u = units(pool(repo))
    assert [m["unit_id"] for m in u["units"]] == ["S1", "S2", "B1", "B2"]
    assert u["units"][0]["libraries"] == ["A1", "A2"] and u["units"][0]["cpg_meth_pct"] == 10.5
    m = u["matrices"]
    assert m[("gene_body", "meth")][:, 0].tolist() == [20, 10, 160, 0]
    assert m[("gene_body", "calls")][:, 0].tolist() == [200, 200, 200, 0]
    assert m[("gene_body", "cpgs")][0, 0] == 9                                   # max of 9 and 7


def test_spearman_handles_ties_and_missing():
    a = np.array([0, 0, 1, 2, np.nan])
    b = np.array([0, 0, 2, 1, 5.0])
    rho, n = spearman(a, b)
    assert n == 4 and rho == pytest.approx(3.5 / 4.5)                            # ranks 1.5,1.5,3,4 vs 1.5,1.5,4,3
    assert spearman(np.array([1.0, np.nan]), np.array([1.0, 2.0])) == (None, 1)


def test_build_landscape_concordance_and_design(repo):
    data, _ = build(repo, repo_url=URL, now="2026-10-25T00:00:00Z")
    assert (data["genes"], data["libraries"], len(data["studies"])) == (4, 6, 2)
    s1 = data["landscape"]["units"][0]
    # S1 gene body: 10%, 5%, 80%; g4 has no CpGs
    assert s1["gene_body"]["genes"] == 3 and s1["gene_body"]["high_frac"] == pytest.approx(0.333)
    assert sum(s1["gene_body"]["bins"]) == 3
    c = data["concordance"]
    assert len(c["rho"]) == 4 and c["rho"][0][0] == 1.0
    within_a = next(b for b in c["blocks"] if b["a"] == b["b"] == "A2020_X_WGBS")
    assert within_a["pairs"] == 1
    between = next(b for b in c["blocks"] if b["a"] != b["b"])
    assert between["pairs"] == 4
    pair = c["study_pairs"][0]
    assert (pair["a"], pair["b"], pair["genes"], pair["shared_tissue"]) == ("A2020_X_WGBS", "B2021_Y_WGBS", 3, [])
    assert pair["both_low"] + pair["both_high"] + pair["opposite"] + pair["other"] == 3
    assert sum(map(sum, pair["density"])) == 3
    caveats = " ".join(data["design"]["caveats"])
    assert "Tissue 'mantle' comes from one study (A2020_X_WGBS)" in caveats
    assert "Every study has the same chemistry" in caveats
    assert any("2 studies (4 samples)" in h for h in data["headline"])


def test_write_tables_are_deterministic(repo, tmp_path):
    data, pooled = build(repo, repo_url=URL, now="2026-10-25T00:00:00Z")
    out = tmp_path / "atlas"
    write(data, pooled, out)
    first = (out / "tables" / "gene_body.meth_pct.tsv.gz").read_bytes()
    write(data, pooled, out)
    assert (out / "tables" / "gene_body.meth_pct.tsv.gz").read_bytes() == first
    rows = gzip.decompress(first).decode().splitlines()
    assert rows[0].split("\t") == ["gene_id", "A1", "A2", "A3", "A4", "B1", "B2"]
    assert rows[4].split("\t")[1] == ""                                         # g4: no calls
    for t in data["tables"]:
        assert (out / t).exists(), t
    meta = (out / "tables" / "sample_metadata.tsv").read_text().splitlines()
    assert meta[0].startswith("library_id\tunit_id\tstudy_id") and len(meta) == 7
    assert json.loads((out / "atlas.json").read_text())["studies"][0]["study_id"] == "A2020_X_WGBS"


def test_cli_writes_atlas(repo, tmp_path, monkeypatch):
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    monkeypatch.setenv("WHALE_REPO", str(repo))
    monkeypatch.setenv("GITHUB_REPOSITORY", "RobertsLab/project-white-whale")
    res = CliRunner().invoke(main, ["build-atlas", "--out", str(tmp_path / "out")])
    assert res.exit_code == 0, res.output
    assert "2 studies, 6 libraries, 4 samples, 4 genes" in res.output
    assert "excluded C2022_Z_WGBS: gene tables not SNP-masked" in res.output
    assert json.loads((tmp_path / "out" / "atlas.json").read_text())["repo_url"] == URL


def test_pages_read_every_key_they_use(repo, repo_root):
    """The pages and the builder agree on field names."""
    data, _ = build(repo, repo_url=URL)
    js = (repo_root / "site" / "atlas" / "atlas.js").read_text()
    land = (repo_root / "site" / "atlas" / "landscape.html").read_text()
    conc = (repo_root / "site" / "atlas" / "concordance.html").read_text()
    for key in ("generated_at", "commit", "repo_url", "genome_assembly", "tables", "excluded"):
        assert key in data and f"data.{key}" in js, key
    for key in ("studies", "libraries", "genes", "headline", "design", "landscape", "bin", "min_cpgs", "low", "high"):
        assert key in data and f"data.{key}" in land, key
    for key in ("unit_id", "study_id", "libraries", "condition", "cpg_meth_pct"):
        assert key in data["landscape"]["units"][0] and f"u.{key}" in land, key
    for key in ("median", "low_frac", "high_frac", "genes"):
        assert key in data["landscape"]["units"][0]["gene_body"], key
    for key in ("rho", "units", "blocks", "study_pairs", "min_genes"):
        assert key in data["concordance"] and key in conc, key
    for key in ("a", "b", "pairs", "mean_rho", "min_rho", "max_rho"):
        assert key in data["concordance"]["blocks"][0] and f"b.{key}" in conc, key
    for key in ("rho", "genes", "shared_tissue", "both_low", "both_high", "opposite", "other", "mean_abs_diff", "density"):
        assert key in data["concordance"]["study_pairs"][0] and f"p.{key}" in conc, key
