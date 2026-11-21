"""Tests for `whale differential` (src/whale/differential.py): the plan, its inputs
for dss.R and the summary. dss.R itself needs the DSS container and is not run here."""
from __future__ import annotations

import copy
import csv
import gzip
import shutil

import pytest
import yaml
from click.testing import CliRunner

from whale.cli import main
from whale.differential import (
    DifferentialError,
    _interactions,
    _terms,
    command,
    plan,
    write_inputs,
    write_summary,
)
from whale.validate import validate_study_file

MASK = {"mask_id": "GCF_X.snpmask.v1", "sha256": "0" * 64, "cpgs_masked": 1}
F = {"factor_original": "o", "factor_standardized": "genotype_or_population"}


def _study(repo, filled, sid, sheet, *, factors, flags=(), differential=None, mask=True):
    s = copy.deepcopy(filled)
    s["study_id"], s["factors"], s["quality_flags"] = sid, factors, list(flags)
    if differential:
        s["differential"] = differential
    (repo / "registry" / "studies" / f"{sid}.yaml").write_text(yaml.safe_dump(s))
    d = repo / "data" / "studies" / sid
    d.mkdir(parents=True)
    (d / "study_summary.yaml").write_text(yaml.safe_dump({
        "study_id": sid, "genome_assembly": "GCF_X", "snp_mask": MASK if mask else None, "nextflow_run": "run_x"}))
    cols = list(sheet[0])
    (d / "samplesheet.tsv").write_text("\t".join(cols) + "\n" + "".join(
        "\t".join(str(r[c]) for c in cols) + "\n" for r in sheet))
    (d / "genes").mkdir()
    with gzip.open(d / "genes" / "gene_body.meth.tsv.gz", "wt") as fh:
        fh.write("\t".join(["gene_id", *(r["sample_id"] for r in sheet)]) + "\n")
    return d


def _crossed(n=2):
    """ploidy x pH, n animals per cell; pH lives in a column named water_ph."""
    rows = []
    for p in ("diploid", "triploid"):
        for ph in ("8.2", "7.7"):
            for k in range(n):
                sid = f"L{p[0]}{ph}{k}"
                rows.append({"sample_id": sid, "sample_accession": f"S{sid}", "ploidy": p, "water_ph": ph})
    return rows


PLOIDY = {"factor_id": "ploidy", **F, "levels": ["diploid", "triploid"], "reference_level": "diploid"}
PH = {"factor_id": "ph", **F, "levels": ["8.2", "7.7"], "reference_level": "8.2"}


@pytest.fixture
def repo(tmp_path, repo_root):
    for rel in ("registry", "schemas"):
        shutil.copytree(repo_root / rel, tmp_path / rel)
    for path in (tmp_path / "registry" / "studies").glob("[!_]*.yaml"):
        path.unlink()
    (tmp_path / "data" / "studies").mkdir(parents=True)
    return tmp_path


def test_formula_terms():
    assert _terms("~ a * b + c") == ["a", "b", "c"]
    assert _terms("~ individual + time") == ["individual", "time"]
    assert _interactions("~ a * b + c") == ["a:b"]
    assert _interactions("~ a + b") == []


def test_default_plan_crossed_two_factors(repo, filled_study):
    _study(repo, filled_study, "R2022_X_WGBS", _crossed(), factors=[PLOIDY, PH])
    p = plan(repo, "R2022_X_WGBS")
    assert [f["column"] for f in p["factors"]] == ["ploidy", "water_ph"]
    assert [(j["job_id"], j["formula"], j["tests"]) for j in p["jobs"]] == [
        ("main", "~ ploidy + ph", ["ploidy", "ph"]),
        ("interaction", "~ ploidy * ph", ["ploidy:ph"])]
    assert len(p["jobs"][0]["units"]) == 8 and p["skipped"] == []


def test_default_plan_no_interaction_when_a_cell_is_thin(repo, filled_study):
    sheet = [r for r in _crossed() if not (r["ploidy"] == "triploid" and r["water_ph"] == "7.7" and r["sample_id"].endswith("1"))]
    _study(repo, filled_study, "R2022_X_WGBS", sheet, factors=[PLOIDY, PH])
    p = plan(repo, "R2022_X_WGBS")
    assert [j["job_id"] for j in p["jobs"]] == ["main"]
    assert any("interaction ploidy:ph" in s for s in p["skipped"])


def test_technical_replicates_are_one_unit(repo, filled_study):
    """Two libraries of one BioSample per group: one unit each, so nothing is testable."""
    sheet = [{"sample_id": f"L{k}", "sample_accession": acc, "origin": o}
             for k, (acc, o) in enumerate((("S1", "w"), ("S1", "w"), ("S2", "i"), ("S2", "i")))]
    origin = {"factor_id": "origin", **F, "levels": ["w", "i"], "reference_level": "w"}
    _study(repo, filled_study, "W2014_X_WGBS", sheet, factors=[origin])
    p = plan(repo, "W2014_X_WGBS")
    assert [u["libraries"] for u in p["units"]] == [["L0", "L1"], ["L2", "L3"]]
    assert p["jobs"] == [] and "fewer than two levels" in p["skipped"][0]


def test_shared_biosample_keeps_libraries_apart(repo, filled_study):
    sheet = [{"sample_id": f"L{k}", "sample_accession": "SHARED", "origin": o}
             for k, o in enumerate("wwii")]
    origin = {"factor_id": "origin", **F, "levels": ["w", "i"], "reference_level": "w"}
    _study(repo, filled_study, "W2020_X_WGBS", sheet, factors=[origin], flags=["shared_biosample"])
    p = plan(repo, "W2020_X_WGBS")
    assert len(p["units"]) == 4
    assert [(j["formula"], len(j["units"])) for j in p["jobs"]] == [("~ origin", 4)]


def test_factor_without_reference_is_not_tested_by_default(repo, filled_study):
    pop = {"factor_id": "population", **F, "levels": ["A", "B"], "reference_level": None}
    sheet = [dict(r, population="AB"[k % 2]) for k, r in enumerate(_crossed())]
    _study(repo, filled_study, "R2022_X_WGBS", sheet, factors=[PLOIDY, pop])
    p = plan(repo, "R2022_X_WGBS")
    assert [j["formula"] for j in p["jobs"]] == ["~ ploidy"]
    assert any("population: no reference level" in s for s in p["skipped"])


def test_no_factors_and_unmasked(repo, filled_study):
    _study(repo, filled_study, "P2022_X_WGBS", _crossed(), factors=[])
    assert plan(repo, "P2022_X_WGBS")["skipped"] == ["no factors registered"]
    _study(repo, filled_study, "U2022_X_WGBS", _crossed(), factors=[PLOIDY], mask=False)
    with pytest.raises(DifferentialError, match="not SNP-masked"):
        plan(repo, "U2022_X_WGBS")


def _repeated():
    """10 animals (5 R, 5 S, one per population), each biopsied at T0 and T1."""
    rows = []
    for pop in range(5):
        for cls in ("S", "R"):
            for t in ("T0", "T1"):
                sid = f"P{pop}{cls}{t}"
                rows.append({"sample_id": sid, "sample_accession": f"S{sid}", "population": f"GT{pop}",
                             "phenotype_class": cls, "time": t, "code": "07"})
    return rows


REPEATED = [{"factor_id": "phenotype_class", **F, "levels": ["S", "R"], "reference_level": "S"},
            {"factor_id": "time", **F, "levels": ["T0", "T1"], "reference_level": "T0"},
            {"factor_id": "population", **F, "levels": [f"GT{k}" for k in range(5)], "reference_level": None}]


def test_named_models_with_strata_and_individual(repo, filled_study):
    spec = {"individual": ["population", "phenotype_class", "code"],
            "models": [{"id": "phenotype", "formula": "~ population + phenotype_class", "within": ["time"]},
                       {"id": "time", "formula": "~ individual + time", "test": ["time"]}]}
    _study(repo, filled_study, "V2025_X_EMSEQ", _repeated(), factors=REPEATED, differential=spec)
    p = plan(repo, "V2025_X_EMSEQ")
    assert [(j["job_id"], j["stratum"], j["tests"], len(j["units"])) for j in p["jobs"]] == [
        ("phenotype.T0", "time=T0", ["phenotype_class"], 10),
        ("phenotype.T1", "time=T1", ["phenotype_class"], 10),
        ("time", "", ["time"], 20)]
    assert len({u["individual"] for u in p["units"]}) == 10

    dest = repo / "out"
    write_inputs(p, dest)
    design = list(csv.DictReader((dest / "design.tsv").read_text().splitlines(), delimiter="\t"))
    assert "individual" in design[0] and design[0]["libraries"] == "P0ST0"
    jobs = list(csv.DictReader((dest / "jobs.tsv").read_text().splitlines(), delimiter="\t"))
    assert jobs[0]["levels"] == "population=GT0,GT1,GT2,GT3,GT4;phenotype_class=S,R"
    assert jobs[2]["levels"] == "time=T0,T1"


def test_named_model_with_unknown_term(repo, filled_study):
    spec = {"models": [{"id": "m", "formula": "~ individual + time"}]}
    _study(repo, filled_study, "V2025_X_EMSEQ", _repeated(), factors=REPEATED, differential=spec)
    with pytest.raises(DifferentialError, match="differential.individual is not set"):
        plan(repo, "V2025_X_EMSEQ")


def test_validate_study_checks_model_terms(tmp_path, repo_root, filled_study):
    s = copy.deepcopy(filled_study)
    s["differential"] = {"models": [{"id": "m", "formula": "~ origin + nonsense"}]}
    path = tmp_path / f"{s['study_id']}.yaml"
    path.write_text(yaml.safe_dump(s))
    res = validate_study_file(path, repo_root)
    assert any("'nonsense' is not a factor_id" in e for e in res.errors)
    s["differential"] = {"models": [{"id": "m", "formula": "~ origin * heat", "within": []}]}
    path.write_text(yaml.safe_dump(s))
    assert validate_study_file(path, repo_root).errors == []


def test_write_summary(repo, filled_study):
    _study(repo, filled_study, "R2022_X_WGBS", _crossed(), factors=[PLOIDY, PH])
    p = plan(repo, "R2022_X_WGBS")
    dest = repo / "out"
    write_inputs(p, dest)
    cols = ["job_id", "analysis", "units", "outliers", "region", "genes_tested", "coef", "n_sig", "up", "down",
            "large", "lambda", "note"]
    rows = [["main", "all_units", 8, "", "gene_body", 100, "ploidytriploid", 3, 2, 1, 0, 1.02, ""],
            ["main", "without_pca_outliers", 7, "SLd8.20", "gene_body", 101, "ploidytriploid", 2, 2, 0, 0, 0.99, ""],
            ["interaction", "all_units", 8, "", "gene_body", 100, "ploidytriploid:ph7.7", 1, "", "", "", 1.1, ""]]
    (dest / "tests.tsv").write_text("\n".join("\t".join(map(str, r)) for r in [cols, *rows]) + "\n")
    (dest / "run_info.tsv").write_text("key\tvalue\ndss_version\t2.58.0\nmin_cpgs\t5\nmin_calls\t10\n")
    s = write_summary(p, dest, today="2026-11-20")
    main_job, inter = s["jobs"]
    t = main_job["all_units"]["regions"]["gene_body"]["tests"]["ploidytriploid"]
    assert t == {"fdr_lt_0.05": 3, "lambda": 1.02, "up": 2, "down": 1, "and_abs_diff_ge_10": 0}
    assert main_job["without_pca_outliers"]["outliers"] == ["SLd8.20"]
    assert inter["all_units"]["regions"]["gene_body"]["tests"]["ploidytriploid:ph7.7"] == {
        "fdr_lt_0.05": 1, "lambda": 1.1}
    assert yaml.safe_load((dest / "summary.yaml").read_text())["method"]["dss_version"] == "2.58.0"


def test_command_engines(repo):
    assert command(repo, "S", "local") == ["Rscript", "workflows/differential/dss.R", "S", "data/studies/S/differential"]
    assert command(repo, "S", "docker")[:3] == ["docker", "run", "--rm"]
    assert command(repo, "S", "apptainer")[0] == "apptainer"
    with pytest.raises(DifferentialError):
        command(repo, "S", "conda")


def test_cli_plan(repo, filled_study, monkeypatch):
    _study(repo, filled_study, "R2022_X_WGBS", _crossed(), factors=[PLOIDY, PH])
    monkeypatch.setenv("WHALE_REPO", str(repo))
    res = CliRunner().invoke(main, ["differential", "--all", "--plan"])
    assert res.exit_code == 0, res.output
    assert "main: ~ ploidy + ph, 8 units; tests ploidy, ph" in res.output
    assert not (repo / "data" / "studies" / "R2022_X_WGBS" / "differential").exists()
