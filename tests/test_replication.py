"""Tests for replication reports (src/whale/replication.py). No network."""
from __future__ import annotations

import copy
import gzip
import shutil

import pytest
import yaml
from click.testing import CliRunner

from whale.cli import main
from whale.replication import evaluate, permutation_p, r_squared, render

SID = "T2026_REP_WGBS"
GENES = [("g1", "LOC1", "heat shock 70 kDa protein 12A-like"), ("g2", "LOC2", "heat shock 70 kDa protein 12A"),
         ("g3", "abc1", "syntenin-1"), ("g4", "LOC4", "uncharacterized LOC4")]
UNITS = ["U1", "U2", "U3", "U4", "U5", "U6"]
LEVELS = {"U1": "a", "U2": "a", "U3": "a", "U4": "b", "U5": "b", "U6": "b"}


def _gz(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def _tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


@pytest.fixture
def repo(tmp_path, repo_root, filled_study):
    for rel in ("registry", "schemas"):
        shutil.copytree(repo_root / rel, tmp_path / rel)
    s = copy.deepcopy(filled_study)
    s["study_id"] = SID
    s["factors"] = [{"factor_id": "grp", "factor_original": "g", "factor_standardized": "other",
                     "levels": ["a", "b"], "reference_level": "a"}]
    (tmp_path / "registry" / "studies" / f"{SID}.yaml").write_text(yaml.safe_dump(s))
    d = tmp_path / "data" / "studies" / SID
    (d / "study_summary.yaml").parent.mkdir(parents=True)
    (d / "study_summary.yaml").write_text(yaml.safe_dump({"study_id": SID, "nextflow_run": "run_x",
                                                         "genome_assembly": "GCF_X", "snp_mask": {"mask_id": "m1"}}))
    _tsv(d / "samplesheet.tsv", [["sample_id", "sample_accession", "grp"],
                                 *[[f"L{u}", u, LEVELS[u]] for u in UNITS]])
    # group b is ~2 points higher overall
    meth = {"U1": 10.0, "U2": 10.2, "U3": 9.9, "U4": 12.1, "U5": 12.0, "U6": 12.3}
    _tsv(d / "qc_summary.tsv", [["sample_id", "cpg_meth_pct", "qc_pass"], *[[f"L{u}", meth[u], "true"] for u in UNITS]])
    libs = [f"L{u}" for u in UNITS]
    # gene bodies 40% methylated, promoters 2%
    for region, frac in (("gene_body", 0.4), ("promoter", 0.02)):
        _gz(d / "genes" / f"{region}.meth.tsv.gz", [["gene_id", *libs], *[[g[0], *[int(100 * frac)] * 6] for g in GENES]])
        _gz(d / "genes" / f"{region}.calls.tsv.gz", [["gene_id", *libs], *[[g[0], *[100] * 6] for g in GENES]])
        _gz(d / "genes" / f"{region}.cpgs.tsv.gz", [["gene_id", *libs], *[[g[0], *[9] * 6] for g in GENES]])
    _gz(tmp_path / "data" / "reference" / "gene_info" / "GCF_X.tsv.gz",
        [["gene_id", "symbol", "description"], *[list(g) for g in GENES]])
    dd = d / "differential"
    dd.mkdir()
    (dd / "summary.yaml").write_text(yaml.safe_dump({"generated_on": "2026-11-21", "jobs": [{
        "job_id": "main", "all_units": {"units": 6, "regions": {
            "gene_body": {"genes_tested": 4, "tests": {"grpb": {"fdr_lt_0.05": 30, "lambda": 1.0, "up": 25, "down": 5}}},
            "promoter": {"genes_tested": 4, "tests": {"grpb": {"fdr_lt_0.05": 0, "lambda": 1.0, "up": 0, "down": 0}}},
        }}}]}))
    _tsv(dd / "jobs.tsv", [["job_id", "model_id", "stratum", "formula", "tests", "levels", "units"],
                           ["main", "main", "", "~ grp", "grp", "grp=a,b", ",".join(UNITS)]])
    _tsv(dd / "design.tsv", [["unit_id", "libraries", "grp"], *[[u, f"L{u}", LEVELS[u]] for u in UNITS]])
    _tsv(dd / "pca.tsv", [["job_id", "unit_id", "pc1", "pc2", "pc1_pct", "pc2_pct", "outlier"],
                          *[["main", u, (5 if LEVELS[u] == "b" else -5) + i * 0.1, (-1) ** i, 40, 10, "FALSE"]
                            for i, u in enumerate(UNITS)]])
    head = ["gene_id", "symbol", "description", "diff_grpb", "stat_grpb", "p_grpb", "fdr_grpb"]
    _gz(dd / "main.gene_body.tsv.gz", [head, ["g1", "LOC1", "", 1.0, 1, 0.2, 0.5], ["g2", "LOC2", "", 3.0, 4, 0.001, 0.01],
                                       ["g3", "abc1", "", -6.0, -5, 1e-6, 1e-4], ["g4", "LOC4", "", 0, 0, 0.9, 0.9]])
    _gz(dd / "main.promoter.tsv.gz", [head, ["g1", "LOC1", "", 0, 0, 0.5, 0.9]])
    return tmp_path


def _claims(repo, claims):
    (repo / "data" / "studies" / SID / "claims.yaml").write_text(yaml.safe_dump(
        {"study_id": SID, "sources": [{"id": "P1", "citation": "A (2020)."}], "claims": claims}))
    return {c["id"]: c for c in evaluate(repo, SID)["claims"]}


def _c(i, kind, **check):
    return {"id": i, "source": "P1", "where": "x", "statement": "s", "kind": kind, "check": check}


def test_levels_regions_and_group_differences(repo):
    r = _claims(repo, [
        _c("L1", "global_level", level_pct=11.5, groups="all"),
        _c("L2", "global_level", level_pct=20, groups={"grp": "a"}),
        _c("L3", "global_level", level_pct=40, groups="all", region="gene_body"),
        _c("R1", "region_contrast", higher="gene_body", lower="promoter"),
        _c("R2", "region_contrast", higher="promoter", lower="genome"),
        _c("G1", "group_difference", factor="grp", higher="b", lower="a"),
        _c("G2", "group_difference", factor="grp", higher="a", lower="b"),
        _c("G3", "group_difference", factor="grp", higher=None, lower=None, expect="no_difference"),
    ])
    assert r["L1"]["verdict"] == "replicated" and r["L2"]["verdict"] == "not_replicated"
    assert r["L3"]["verdict"] == "replicated"                    # per-gene level, 40%
    assert r["R1"]["verdict"] == "replicated" and r["R2"]["verdict"] == "not_replicated"
    assert r["G1"]["verdict"] == "replicated" and r["G1"]["p"] < 0.1
    assert r["G2"]["verdict"] == "not_replicated"
    assert r["G3"]["verdict"] == "inconclusive"                  # 3 vs 3 cannot show a difference


def test_fractions_without_a_cutoff_stay_open(repo):
    r = _claims(repo, [_c("F1", "methylated_fraction", unit="genes", region="gene_body", threshold_pct=None,
                          fraction_pct=50),
                       _c("F2", "methylated_fraction", unit="genes", region="gene_body", threshold_pct=30,
                          fraction_pct=100),
                       _c("F3", "methylated_fraction", unit="cpgs", fraction_pct=16)])
    assert r["F1"]["verdict"] == "inconclusive" and r["F1"]["ge10"] == 100.0
    assert r["F2"]["verdict"] == "replicated"
    assert r["F3"]["verdict"] == "not_testable"


def test_differential_counts_direction_genes_and_pca(repo):
    r = _claims(repo, [
        _c("D1", "dm_count", factor="grp", level="b", stratum=None, n=40, unit="DMGs"),
        _c("D2", "dm_count", factor="grp", level="b", stratum=None, n=400, unit="DMGs"),
        _c("D3", "dm_count", factor="grp", level="b", stratum=None, n=500, unit="DMRs"),
        _c("D4", "dm_count", factor="grp", level="b", stratum={"tissue": "gill"}, n=40, unit="DMGs"),
        _c("V1", "dm_direction", factor="grp", level="b", stratum=None, direction="hyper"),
        _c("N1", "dm_genes", factor="grp", level="b", stratum=None, genes=[
            {"name": "x / syntenin-1", "direction": "hypo"},
            {"name": "heat shock 70 kDa protein 12A (HSPA12A)", "direction": "hyper"},
            {"name": "no such protein anywhere", "direction": None}]),
        _c("P1", "clustering", factor="grp", component="PC1", statement="separate by grp"),
        _c("P2", "clustering", factor="grp", component="PC2", statement="no separation by grp"),
        _c("P3", "clustering", factor="grp", component="hierarchical", statement=""),
    ])
    assert r["D1"]["verdict"] == "replicated" and r["D2"]["verdict"] == "partly"
    assert r["D3"]["verdict"] == "partly"                        # other units: presence only
    assert r["D4"]["verdict"] == "partly" and r["D4"]["approximate"]   # pooled over tissue
    assert r["V1"]["verdict"] == "replicated"
    genes = {g["name"]: g for g in r["N1"]["genes"]}
    assert genes["x / syntenin-1"]["status"] == "replicated"
    # two genes match HSPA12A (one "-like"): the best p (0.001) is multiplied by 4 (2 genes x 2 regions)
    assert genes["heat shock 70 kDa protein 12A (HSPA12A)"]["status"] == "replicated"
    assert "best of 4 tests" in genes["heat shock 70 kDa protein 12A (HSPA12A)"]["detail"]
    assert genes["no such protein anywhere"]["status"] == "not_found"
    assert r["P1"]["verdict"] == "replicated" and r["P2"]["verdict"] == "replicated"
    assert r["P3"]["verdict"] == "not_testable"


def test_inflated_statistics_cap_counts(repo):
    p = repo / "data" / "studies" / SID / "differential" / "summary.yaml"
    doc = yaml.safe_load(p.read_text())
    doc["jobs"][0]["all_units"]["regions"]["gene_body"]["tests"]["grpb"]["lambda"] = 1.8
    p.write_text(yaml.safe_dump(doc))
    r = _claims(repo, [_c("D1", "dm_count", factor="grp", level="b", stratum=None, n=40, unit="DMGs")])
    assert r["D1"]["verdict"] == "partly" and r["D1"]["inflated"]


def test_without_differential_claims_are_pending(repo):
    shutil.rmtree(repo / "data" / "studies" / SID / "differential")
    r = _claims(repo, [_c("D1", "dm_count", factor="grp", level="b", stratum=None, n=40, unit="DMGs"),
                       {"id": "X", "kind": "not_testable", "statement": "s", "check": {"reason": "needs_rnaseq"}}])
    assert r["D1"]["verdict"] == "pending" and r["X"]["verdict"] == "not_testable"


def test_report_and_cli(repo, monkeypatch):
    _claims(repo, [_c("L1", "global_level", level_pct=11.5, groups="all")])
    text = render(evaluate(repo, SID))
    assert "| L1 |" in text and "**replicated**" in text and "run_x" in text
    monkeypatch.chdir(repo)
    res = CliRunner().invoke(main, ["replicate"])
    assert res.exit_code == 0 and "1 replicated" in res.output
    assert (repo / "data" / "studies" / SID / "replication.md").exists()


def test_statistics():
    assert permutation_p([1, 1.1, 0.9], [5, 5.1, 4.9]) < 0.15
    assert permutation_p([1], [2, 3]) is None
    assert r_squared([0, 0.1, 10, 10.1], ["a", "a", "b", "b"]) > 0.99
    assert r_squared([1, 2, 3], ["a", "a", "a"]) is None
