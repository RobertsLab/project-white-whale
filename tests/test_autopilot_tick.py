"""`whale autopilot tick --dry-run` end to end, on a fixture repository with study branches."""
import subprocess

import pytest
import yaml
from autopilot_fakes import nf_log
from click.testing import CliRunner
from conftest import git

from whale.autopilot.digest import lines, markdown
from whale.autopilot.tick import AutopilotError, tick
from whale.cli import main

BACKLOG = (
    "accession\tsuggested_study_id\tchemistry\truns\tfastq_gb\tpaired_rnaseq\tsuggested_profile\tsource\tnotes\n"
    "PRJNA999901\tTESTER2099_DONE_WGBS\twgbs\t4\t10\tnone\traven\ttest\t\n"
    "PRJNA999902\tTESTER2099_FAILED_WGBS\twgbs\t4\t10\tnone\tdesktop\ttest\t\n"
    "PRJNA999903\tTESTER2099_CLAIMED_WGBS\twgbs\t4\t10\tnone\traven\ttest\tAREE registration to import.\n"
    "PRJNA999904\tTESTER2099_NEXT_WGBS\twgbs\t4\t10\tnone\traven\ttest\t\n"
    "PRJNA999905\tTESTER2099_LATER_WGBS\twgbs\t4\t10\tnone\traven\ttest\t\n"
    "PRJNA999906\t(name on claim)\twgbs\t4\t10\tnone\traven\ttest\t\n"
)
MACHINES = {"version": "0.1.0", "machines": {
    "Hete": {"profile": "desktop", "max_concurrent_runs": 1, "env": {"WHALE_WORKDIR": "/w/work"}},
    "raven": {"profile": "raven", "max_concurrent_runs": 2, "gh": False},
}}
QC = ("sample_id\tmapping_pct\tnon_cpg_meth_pct\tqc_pass\n"
      + "".join(f"S{i}\t55.{i}\t0.3\ttrue\n" for i in range(4)))


def _status(sid, acc, stage, machine, profile):
    return yaml.safe_dump({"study_id": sid, "accession": acc, "part": None, "stage": stage, "profile": profile,
                           "machine": machine, "history": [{"stage": stage, "note": None}]})


def _branch(clone, acc, sid, stage, machine, profile, extra=None, merge=False):
    branch = f"study/{acc}/{sid}"
    git("checkout", "--quiet", "-b", branch, "main", cwd=clone)
    files = {f"data/studies/{sid}/status.yaml": _status(sid, acc, stage, machine, profile),
             f"registry/studies/{sid}.yaml": yaml.safe_dump({"study_id": sid, "chemistry": "wgbs",
                                                              "accessions": {"bioproject": acc}}),
             f"config/{sid}.config": f"params {{\n    study_id = '{sid}'\n    chemistry = 'wgbs'\n}}\n",
             **(extra or {})}
    for rel, text in files.items():
        (clone / rel).parent.mkdir(parents=True, exist_ok=True)
        (clone / rel).write_text(text)
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", f"{sid}: {stage}", cwd=clone)
    git("push", "--quiet", "origin", branch, cwd=clone)
    git("checkout", "--quiet", "main", cwd=clone)
    if merge:
        git("merge", "--quiet", "--no-ff", "-m", f"merge {sid}", branch, cwd=clone)
        git("push", "--quiet", "origin", "main", cwd=clone)
    return branch


@pytest.fixture
def fleet(claim_env):
    """main with a backlog; one merged study, one failed run on this machine, one claimed."""
    clone, tmp = claim_env["clone"], claim_env["tmp"]
    for p in (clone / "registry/studies").glob("*.yaml"):
        if not p.name.startswith("_"):
            p.unlink()
    (clone / "registry/backlog.tsv").write_text(BACKLOG)
    (clone / "config/autopilot.yaml").write_text(yaml.safe_dump(MACHINES))
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", "fixture backlog", cwd=clone)
    git("push", "--quiet", "origin", "main", cwd=clone)

    summary = yaml.safe_dump({"libraries": 4, "passed": 4, "genome_assembly": "GCF_963853765.1",
                              "snp_mask": {"mask_id": "m"}, "methylomes": {"on_data_host": True}})
    _branch(clone, "PRJNA999901", "TESTER2099_DONE_WGBS", "summarized", "raven", "raven", merge=True,
            extra={"data/studies/TESTER2099_DONE_WGBS/study_summary.yaml": summary,
                   "data/studies/TESTER2099_DONE_WGBS/qc_summary.tsv": QC})
    failed = _branch(clone, "PRJNA999902", "TESTER2099_FAILED_WGBS", "fastq_staged", "Hete", "desktop")
    _branch(clone, "PRJNA999903", "TESTER2099_CLAIMED_WGBS", "claimed", "raven", "raven")

    # The failed study's worktree on this machine, with its run log.
    wt = tmp / "white-whale-TESTER2099_FAILED_WGBS"
    git("worktree", "add", "--quiet", str(wt), failed, cwd=clone)
    (wt / ".nextflow.log").write_text(nf_log(failure={"process": "BISMARK_DEDUPLICATE", "tag": "S1", "exit": 137}))
    return claim_env


def by_key(res):
    return {d.key: d for d in res.decisions}


def test_tick_dry_run_on_fixture_repo(fleet):
    clone = fleet["clone"]
    head_before = git("rev-parse", "HEAD", cwd=clone)
    status_before = git("status", "--porcelain", cwd=clone)
    res = tick(clone, dry_run=True, here="Hete")
    d = by_key(res)
    assert d["TESTER2099_DONE_WGBS"].state == "merged" and d["TESTER2099_DONE_WGBS"].kind == "done"
    assert d["TESTER2099_FAILED_WGBS"].playbook == "oom_killed"
    claimed = d["TESTER2099_CLAIMED_WGBS"]
    assert claimed.kind == "blocked" and "AREE" in claimed.action
    # raven has 1 study in progress and room for 2: the next backlog study is claimed, one per tick.
    assert d["TESTER2099_NEXT_WGBS"].next_state == "claimed" and d["TESTER2099_NEXT_WGBS"].machine == "raven"
    assert d["TESTER2099_LATER_WGBS"].kind == "waiting"
    assert d["PRJNA999906"].kind == "blocked"
    # Nothing changed: same commit, same working tree, no new branches pushed.
    assert git("rev-parse", "HEAD", cwd=clone) == head_before
    assert git("status", "--porcelain", cwd=clone) == status_before
    assert len(lines(res)) == 6
    md = markdown(res, today="2026-11-21")
    assert "## Would advance" in md and "## Blocked: needs a person" in md and "oom_killed" in md


def test_tick_refuses_without_dry_run(fleet):
    with pytest.raises(AutopilotError):
        tick(fleet["clone"], dry_run=False)


def test_tick_respects_pause(fleet):
    (fleet["clone"] / "AUTOPILOT_PAUSED").write_text("klone maintenance\n")
    res = tick(fleet["clone"], dry_run=True, here="Hete")
    assert res.paused == "klone maintenance" and res.decisions == []
    assert "Paused" in markdown(res)


def test_cli(fleet, monkeypatch, tmp_path):
    monkeypatch.setenv("WHALE_REPO", str(fleet["clone"]))
    runner = CliRunner()
    r = runner.invoke(main, ["autopilot", "tick"])
    assert r.exit_code != 0 and "--dry-run" in r.output
    digest = tmp_path / "digest.md"
    r = runner.invoke(main, ["autopilot", "tick", "--dry-run", "--no-fetch", "--machine", "Hete",
                             "--digest", str(digest)])
    assert r.exit_code == 0, r.output
    assert "TESTER2099_FAILED_WGBS" in r.output and "# Autopilot digest" in digest.read_text()
    (fleet["clone"] / "AUTOPILOT_PAUSED").write_text("")
    r = runner.invoke(main, ["autopilot", "tick", "--dry-run", "--no-fetch"])
    assert r.exit_code == 0 and r.output.startswith("Paused")


def test_no_remote_still_works(tmp_path):
    """A checkout without a remote (main only) reads the backlog from main."""
    repo = tmp_path / "solo"
    subprocess.run(["git", "init", "--quiet", "-b", "main", str(repo)], check=True)
    (repo / "registry").mkdir()
    (repo / "registry/backlog.tsv").write_text(BACKLOG)
    git("-c", "user.name=t", "-c", "user.email=t@e", "add", "-A", cwd=repo)
    git("-c", "user.name=t", "-c", "user.email=t@e", "commit", "--quiet", "-m", "x", cwd=repo)
    res = tick(repo, dry_run=True, do_fetch=False, here="Hete", machines_path=tmp_path / "none.yaml")
    assert len(res.decisions) == 6 and any("no machines" in w for w in res.warnings)
