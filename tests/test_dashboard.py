"""Tests for `whale build-dashboard` (src/whale/dashboard.py). No network access."""
from __future__ import annotations

import json
import shutil
import subprocess

import pytest
import yaml
from click.testing import CliRunner

from whale.cli import main
from whale.dashboard import build, load_findings, mirror_summary, write

URL = "https://github.com/RobertsLab/project-white-whale"


@pytest.fixture
def repo(tmp_path, repo_root, filled_study):
    """A copy of the registry with one study registered on main (WANG2020) and its status."""
    for rel in ("registry", "schemas"):
        shutil.copytree(repo_root / rel, tmp_path / rel)
    # Only the study registered below: studies registered on the real main must not change the counts.
    for path in (tmp_path / "registry" / "studies").glob("[!_]*.yaml"):
        path.unlink()
    (tmp_path / "data" / "reference").mkdir(parents=True)
    (tmp_path / "registry" / "studies" / "WANG2020_TIDAL_HEAT_WGBS.yaml").write_text(yaml.safe_dump(filled_study))
    return tmp_path


def _claim(acc, sid, stage, machine="raven", part=None):
    return {"study_id": sid, "part": part or "", "accession": acc, "stage": stage, "profile": machine,
            "machine": machine, "claimed_by": "t", "claimed_on": "2026-10-23", "updated": "2026-10-24",
            "branch": f"study/{acc}/{sid}" + (f".part-{part}" if part else "")}


def _by_acc(data):
    return {s["accession"]: s for s in data["studies"]}


def test_build_merges_backlog_registry_claims_and_mirror(repo):
    claims = [_claim("PRJNA173440", "WANG2014_FIRST_METHYLOMES_WGBS", "fastq_staged")]
    ledger = {"SRR1_1.fastq.gz": {"run_accession": "SRR1", "bytes": "2000000000", "source": "ncbi_sra",
                                  "mirrored_on": "2026-10-23"},
              "SRR1_2.fastq.gz": {"run_accession": "SRR1", "bytes": "1000000000", "source": "ncbi_sra",
                                  "mirrored_on": "2026-10-24"}}
    data = build(repo, claims=claims, mirrors={"PRJNA173440": ledger, "PRJNA682817": None},
                 repo_url=URL, now="2026-10-24T00:00:00Z")
    s = _by_acc(data)
    assert len(data["studies"]) == 21 and data["totals"]["backlog"] == 21

    w = s["PRJNA173440"]
    assert (w["state"], w["stage"], w["machine"]) == ("in progress", "fastq_staged", "raven")
    assert w["stage_index"] == [st["id"] for st in data["stages"]].index("fastq_staged")
    assert w["branch_url"] == f"{URL}/tree/study/PRJNA173440/WANG2014_FIRST_METHYLOMES_WGBS"
    assert w["mirror"] == {"files": 2, "runs": 1, "gb": 3.0, "sources": {"ncbi_sra": 2}, "last_mirrored": "2026-10-24"}

    r = s["PRJNA562805"]   # registered on main, not claimed
    assert (r["registered"], r["state"], r["stage"]) == (True, "registered", None)
    assert r["paired_rnaseq"] == "same_animals" and r["goals"] == [1, 2] and r["tissues"] == ["gill"]

    n = s["PRJNA807732"]   # backlog only, no name yet
    assert n["study_id"] is None and n["state"] == "not started"

    t = data["totals"]
    assert (t["registered"], t["in_progress"], t["done"]) == (1, 1, 0)
    assert t["mirrored_gb"] == 3.0 and t["mirrored_studies"] == 1
    assert t["by_stage"]["fastq_staged"] == 1
    assert data["claims_available"] and data["mirror_available"]
    assert "excluded" not in [st["id"] for st in data["stages"]]


def test_merged_study_reads_its_status_on_main(repo):
    status_dir = repo / "data" / "studies" / "WANG2020_TIDAL_HEAT_WGBS"
    status_dir.mkdir(parents=True)
    (status_dir / "status.yaml").write_text(yaml.safe_dump({
        "study_id": "WANG2020_TIDAL_HEAT_WGBS", "stage": "summarized", "profile": "raven", "machine": "raven",
        "history": [{"stage": "summarized", "date": "2026-11-02", "machine": "raven"}]}))
    data = build(repo, claims=[], mirrors={}, repo_url=URL)
    r = _by_acc(data)["PRJNA562805"]
    assert (r["state"], r["stage"], r["updated"]) == ("done", "summarized", "2026-11-02")
    assert data["totals"]["done"] == 1


def test_parts_report_the_furthest_stage(repo):
    claims = [_claim("PRJEB60400", "POMS2023_ADAPTATION_WGBS", "aligned", part="gill"),
              _claim("PRJEB60400", "POMS2023_ADAPTATION_WGBS", "fastq_staged", "klone", part="mantle")]
    p = _by_acc(build(repo, claims=claims, mirrors=None, repo_url=URL))["PRJEB60400"]
    assert p["stage"] == "aligned" and p["parts"] == ["gill", "mantle"]


def test_unavailable_sources_are_flagged(repo):
    data = build(repo, claims=None, mirrors=None, repo_url=URL)
    assert not data["claims_available"] and not data["mirror_available"]
    assert all(s["mirror"] is None for s in data["studies"])


def test_findings(repo):
    (repo / "site").mkdir()
    (repo / "site" / "findings.yaml").write_text(yaml.safe_dump({"findings": [
        {"date": "2026-11-01", "title": "A", "summary": "x ", "goal": 1, "studies": ["S1"]},
        {"date": "2026-12-01", "title": "B", "summary": "y"}]}))
    f = load_findings(repo)
    assert [x["title"] for x in f] == ["B", "A"] and f[1]["summary"] == "x" and f[0]["studies"] == []
    (repo / "site" / "findings.yaml").write_text(yaml.safe_dump({"findings": [{"date": "2026-11-01", "title": "A"}]}))
    with pytest.raises(ValueError, match="lacks"):
        load_findings(repo)


def test_committed_findings_file_is_valid(repo_root):
    load_findings(repo_root)


def test_mirror_summary_empty():
    assert mirror_summary(None) is None and mirror_summary({}) is None


def test_cli_writes_data_json(repo, tmp_path, monkeypatch):
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    monkeypatch.setenv("WHALE_REPO", str(repo))
    monkeypatch.setenv("GITHUB_REPOSITORY", "RobertsLab/project-white-whale")
    out = tmp_path / "out" / "data.json"
    res = CliRunner().invoke(main, ["build-dashboard", "--no-claims", "--no-mirror", "--out", str(out)])
    assert res.exit_code == 0, res.output
    data = json.loads(out.read_text())
    assert data["repo_url"] == URL and data["totals"]["backlog"] == 21


def test_write(tmp_path):
    p = write({"a": 1}, tmp_path / "x" / "data.json")
    assert json.loads(p.read_text()) == {"a": 1}


def test_page_reads_every_key_it_uses(repo_root):
    """The page and the builder agree on field names."""
    html = (repo_root / "site" / "index.html").read_text()
    data = build(repo_root, claims=[], mirrors={}, repo_url=URL)
    for key in ("generated_at", "commit", "repo_url", "claims_available", "mirror_available", "stages",
                "totals", "studies", "findings"):
        assert key in data and f"data.{key}" in html, key
    for key in data["totals"]:
        assert f"t.{key}" in html, key
    for key in ("study_id", "accession", "chemistry", "runs", "fastq_gb", "mirror", "stage", "stage_index",
                "state", "machine", "profile", "paired_rnaseq", "parts", "branch_url", "updated", "notes"):
        assert f"s.{key}" in html, key
