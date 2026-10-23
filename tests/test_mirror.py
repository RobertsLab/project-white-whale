"""Offline tests for `whale mirror-fastq` (src/whale/mirror.py).

A fake ENA serves a two-run project (read-run report and FASTQ bytes, with
Range support and injectable faults); the data host is a local directory.
No network access.
"""
from __future__ import annotations

import hashlib
import io
import shutil
import urllib.error
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

import whale.mirror as mirror_mod
from whale.cli import main
from whale.mirror import (
    LocalStore,
    MirrorError,
    SSHStore,
    download,
    mirror,
    plan_mirror,
    wanted_files,
)

ACC = "PRJNA000001"
BASE = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR000"


def _fq(run, mate):
    return (f"@{run}.{mate}\nACGT\n+\nIIII\n" * 500).encode()


class FakeENA:
    def __init__(self):
        self.files = {f"{run}_{m}.fastq.gz": _fq(run, m) for run in ("SRR1", "SRR2") for m in (1, 2)}
        self.cut_once: set[str] = set()     # fail these mid-stream once
        self.corrupt: set[str] = set()      # serve wrong bytes every time
        self.ignore_range = False
        self.requests: list[tuple[str, str]] = []

    def md5(self, name):
        return hashlib.md5(self.files[name]).hexdigest()

    def report(self, accession):
        assert accession == ACC
        rows = []
        for run, strategy in (("SRR1", "Bisulfite-Seq"), ("SRR2", "Bisulfite-Seq"), ("SRR9", "RNA-Seq")):
            names = [f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"]
            have = [n for n in names if n in self.files]
            rows.append({
                "run_accession": run, "library_strategy": strategy, "library_layout": "PAIRED",
                "tax_id": "29159",
                "fastq_ftp": ";".join(f"{BASE}/{run}/{n}" for n in have),
                "fastq_md5": ";".join(self.md5(n) for n in have),
                "fastq_bytes": ";".join(str(len(self.files[n])) for n in have),
            })
        return rows

    def open(self, req):
        name = req.full_url.rsplit("/", 1)[-1]
        rng = req.get_header("Range") or ""
        self.requests.append((name, rng))
        assert req.full_url.startswith("https://")
        if name not in self.files:
            raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, None)
        body = b"x" * len(self.files[name]) if name in self.corrupt else self.files[name]
        status = 200
        if rng and not self.ignore_range:
            body, status = body[int(rng.removeprefix("bytes=").rstrip("-")):], 206
        if name in self.cut_once:
            self.cut_once.discard(name)
            return _Resp(body[: len(body) // 3], status, fail=True)
        return _Resp(body, status)


class _Resp:
    def __init__(self, body, status, fail=False):
        self.buf, self.status, self.fail = io.BytesIO(body), status, fail

    def read(self, n):
        chunk = self.buf.read(min(n, 700))
        if not chunk and self.fail:
            raise ConnectionResetError("cut")
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr(mirror_mod, "_sleep", lambda s: None)


@pytest.fixture
def ena(monkeypatch):
    fake = FakeENA()
    monkeypatch.setattr(mirror_mod, "fetch_run_report", fake.report)
    monkeypatch.setattr(mirror_mod, "_urlopen", fake.open)
    return fake


@pytest.fixture
def repo(tmp_path, repo_root):
    r = tmp_path / "repo"
    shutil.copytree(repo_root / "registry", r / "registry")
    shutil.copytree(repo_root / "data", r / "data", ignore=shutil.ignore_patterns("raw", "studies"))
    return r


def _run(repo, tmp_path, **kw):
    return mirror(repo, [ACC], LocalStore(tmp_path / "host"), staging=tmp_path / "staging",
                  log=lambda *_: None, today="2026-10-23", **kw)


def test_wanted_files_keeps_methylation_runs(repo, ena):
    files, notes = wanted_files(repo, ACC)
    assert [f["file"] for f in files] == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz", "SRR2_1.fastq.gz", "SRR2_2.fastq.gz"]
    assert files[0]["url"] == f"ftp://{BASE}/SRR1/SRR1_1.fastq.gz"
    assert any("RNA-Seq" in n for n in notes)


def test_registered_study_uses_its_runs(repo, ena):
    (repo / "registry/studies/TEST_STUDY.yaml").write_text(yaml.safe_dump(
        {"study_id": "TEST_STUDY", "accessions": {"bioproject": ACC}}))
    d = repo / "data/studies/TEST_STUDY"
    d.mkdir(parents=True)
    (d / "fastq_manifest.tsv").write_text("run_accession\tmate\tfile\nSRR2\t1\tSRR2_1.fastq.gz\nSRR7\t1\tx\n")
    files, notes = wanted_files(repo, ACC)
    assert {f["run_accession"] for f in files} == {"SRR2"}
    assert any("SRR7" in n for n in notes)


def test_mirror_end_to_end_and_rerun(repo, ena, tmp_path):
    res = _run(repo, tmp_path)
    assert res["failed"] == [] and len(res["mirrored"]) == 4
    host = tmp_path / "host/raw" / ACC
    assert (host / "SRR2_2.fastq.gz").read_bytes() == ena.files["SRR2_2.fastq.gz"]
    ledger = LocalStore(tmp_path / "host").read_ledger(ACC)
    assert sorted(ledger) == sorted(ena.files)
    assert ledger["SRR1_1.fastq.gz"]["md5"] == ena.md5("SRR1_1.fastq.gz")
    assert ledger["SRR1_1.fastq.gz"]["mirrored_on"] == "2026-10-23"
    assert not list((tmp_path / "staging" / ACC).glob("*.fastq.gz*"))

    # Nothing left to do on a rerun.
    ena.requests.clear()
    assert _run(repo, tmp_path)["mirrored"] == [] and ena.requests == []

    # A file missing from the host, or reprocessed by ENA, is fetched again.
    (host / "SRR1_1.fastq.gz").unlink()
    ena.files["SRR2_1.fastq.gz"] = b"@new\nA\n+\nI\n"
    plan = plan_mirror(repo, [ACC], LocalStore(tmp_path / "host"))[0]
    assert {f["file"] for f in plan["todo"]} == {"SRR1_1.fastq.gz", "SRR2_1.fastq.gz"}
    assert any("MD5 changed" in n for n in plan["notes"])
    assert len(_run(repo, tmp_path)["mirrored"]) == 2
    assert LocalStore(tmp_path / "host").read_ledger(ACC)["SRR2_1.fastq.gz"]["md5"] == ena.md5("SRR2_1.fastq.gz")


def test_interrupted_download_resumes(repo, ena, tmp_path):
    ena.cut_once = {"SRR1_1.fastq.gz"}
    res = _run(repo, tmp_path, jobs=1)
    assert res["failed"] == []
    ranges = [r for n, r in ena.requests if n == "SRR1_1.fastq.gz"]
    assert ranges[0] == "" and ranges[1].startswith("bytes=") and ranges[1] != "bytes=0-"
    assert (tmp_path / "host/raw" / ACC / "SRR1_1.fastq.gz").read_bytes() == ena.files["SRR1_1.fastq.gz"]


def test_partial_file_from_earlier_run_resumes(ena, tmp_path):
    name = "SRR1_2.fastq.gz"
    dest = tmp_path / name
    dest.with_name(name + ".part").write_bytes(ena.files[name][:1000])
    md5 = download(f"https://{BASE}/SRR1/{name}", dest, expected_bytes=len(ena.files[name]))
    assert md5 == ena.md5(name) and ena.requests == [(name, "bytes=1000-")]

    # A server that ignores Range sends the whole file; it is not appended.
    ena.ignore_range, ena.requests = True, []
    dest.unlink()
    dest.with_name(name + ".part").write_bytes(ena.files[name][:1000])
    assert download(f"https://{BASE}/SRR1/{name}", dest) == ena.md5(name)
    assert dest.read_bytes() == ena.files[name]


def test_md5_mismatch_is_not_mirrored(repo, ena, tmp_path):
    ena.corrupt = {"SRR2_1.fastq.gz"}
    res = _run(repo, tmp_path)
    assert [f for f, _ in res["failed"]] == [f"{ACC}/SRR2_1.fastq.gz"]
    assert "does not match ENA" in res["failed"][0][1]
    assert not (tmp_path / "host/raw" / ACC / "SRR2_1.fastq.gz").exists()
    assert "SRR2_1.fastq.gz" not in LocalStore(tmp_path / "host").read_ledger(ACC)
    assert not (tmp_path / "staging" / ACC / "SRR2_1.fastq.gz").exists()
    assert len(res["mirrored"]) == 3


def test_max_gb_limits_a_run(repo, ena, tmp_path):
    one_file_gb = len(ena.files["SRR1_1.fastq.gz"]) / 1e9
    res = _run(repo, tmp_path, max_gb=one_file_gb * 1.5, jobs=1)
    assert len(res["mirrored"]) == 2


def test_second_run_exits_while_locked(repo, ena, tmp_path):
    with mirror_mod._Lock(tmp_path / "staging" / "mirror.lock") as got:
        assert got
        assert _run(repo, tmp_path)["locked"] is True
    assert _run(repo, tmp_path)["locked"] is False


def test_ssh_target_must_name_host_and_path():
    with pytest.raises(MirrorError):
        SSHStore("gannet.fish.washington.edu")
    s = SSHStore("sr320@gannet.fish.washington.edu:/volume1/v1_web/project-white-whale/")
    assert (s.host, s._dir(ACC)) == ("sr320@gannet.fish.washington.edu",
                                     f"/volume1/v1_web/project-white-whale/raw/{ACC}")


def test_cli_plan_then_mirror(repo, ena, tmp_path, monkeypatch):
    monkeypatch.setenv("WHALE_REPO", str(repo))
    host, staging = tmp_path / "host", tmp_path / "staging"
    runner = CliRunner()
    out = runner.invoke(main, ["mirror-fastq", "--accession", ACC, "--dest", str(host), "--plan"])
    assert out.exit_code == 0, out.output
    assert "4" in out.output and "RNA-Seq" in out.output and not Path(host / "raw").exists()

    out = runner.invoke(main, ["mirror-fastq", "--accession", ACC, "--dest", str(host),
                               "--staging", str(staging)])
    assert out.exit_code == 0, out.output
    assert "Mirrored 4 file(s)" in out.output and "complete" in out.output

    ena.corrupt = {"SRR1_1.fastq.gz"}
    (host / "raw" / ACC / "SRR1_1.fastq.gz").unlink()
    out = runner.invoke(main, ["mirror-fastq", "--accession", ACC, "--dest", str(host),
                               "--staging", str(staging)])
    assert out.exit_code == 1 and "1 failed" in out.output
