"""Offline tests for `whale mirror-fastq` (src/whale/mirror.py).

A fake ENA serves a two-run project (read-run report and FASTQ bytes, with
Range support and injectable faults); a fake NCBI serves the locator and
.sra files, and a fake fasterq-dump unpacks them. The data host is a local
directory. No network access.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import shutil
import sys
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
    record_in_ledger,
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
        self.short_once: set[str] = set()   # end these early once, without an error
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
                "tax_id": "29159", "read_count": "500", "base_count": "4000",
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
        if name in self.short_once:
            self.short_once.discard(name)
            return _Resp(body[: len(body) // 3], status)
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
    kw.setdefault("source", "ena")
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


def test_short_response_resumes(ena, tmp_path):
    """A response that ends early without an error is resumed, not restarted from zero."""
    name = "SRR1_1.fastq.gz"
    ena.short_once = {name}
    dest = tmp_path / name
    md5 = download(f"https://{BASE}/SRR1/{name}", dest, expected_bytes=len(ena.files[name]))
    assert md5 == ena.md5(name) and dest.read_bytes() == ena.files[name]
    third = len(ena.files[name]) // 3
    assert ena.requests == [(name, ""), (name, f"bytes={third}-")]


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
    assert [f for f, _ in res["failed"]] == [f"{ACC}/SRR2"]
    assert "does not match ENA" in res["failed"][0][1]
    assert not (tmp_path / "host/raw" / ACC / "SRR2_1.fastq.gz").exists()
    assert "SRR2_1.fastq.gz" not in LocalStore(tmp_path / "host").read_ledger(ACC)
    assert not (tmp_path / "staging" / ACC / "SRR2_1.fastq.gz").exists()
    assert len(res["mirrored"]) == 2


def test_max_gb_limits_a_run(repo, ena, tmp_path):
    one_file_gb = len(ena.files["SRR1_1.fastq.gz"]) / 1e9
    res = _run(repo, tmp_path, max_gb=one_file_gb * 1.5, jobs=1)
    assert len(res["mirrored"]) == 2


def test_second_run_exits_while_locked(repo, ena, tmp_path):
    lock = tmp_path / "staging" / "mirror.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text("12345\n")   # left by an earlier run; replaced, not appended to
    with mirror_mod._Lock(lock) as got:
        assert got
        holder = lock.read_text()
        assert holder == f"{os.getpid()}\n"
        assert _run(repo, tmp_path)["locked"] is True
        # The run that was turned away leaves the holder's PID in place.
        assert lock.read_text() == holder
    assert _run(repo, tmp_path)["locked"] is False


def _row(name):
    return {"run_accession": name.split("_")[0], "file": name, "md5": "0" * 32, "bytes": "1"}


def test_ledger_update_keeps_rows_written_by_others(tmp_path):
    """Writers with no shared memory (two machines, two staging dirs) never drop each other's rows."""
    import threading

    root = tmp_path / "host"
    names = [f"SRR{w}{i:02d}_1.fastq.gz" for w in (1, 2, 3) for i in range(15)]

    def writer(prefix):
        store = LocalStore(root)          # a store of its own, as another process would have
        for name in names:
            if name.startswith(prefix):
                record_in_ledger(store, ACC, _row(name), tmp_path / f"staging-{prefix}")

    threads = [threading.Thread(target=writer, args=(f"SRR{w}",)) for w in (1, 2, 3)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    assert sorted(LocalStore(root).read_ledger(ACC)) == sorted(names)
    assert not (root / "raw" / ACC / mirror_mod.LEDGER_LOCK).exists()


def test_held_ledger_lock_times_out_with_its_path(tmp_path):
    store = LocalStore(tmp_path / "host")
    lock = tmp_path / "host" / "raw" / ACC / mirror_mod.LEDGER_LOCK
    lock.mkdir(parents=True)
    with pytest.raises(MirrorError, match="remove that directory"):
        record_in_ledger(store, ACC, _row("SRR1_1.fastq.gz"), tmp_path / "staging", wait=0)
    assert lock.is_dir()     # someone else's lock is left alone


def test_ledger_lock_is_released_when_the_write_fails(tmp_path, monkeypatch):
    store = LocalStore(tmp_path / "host")

    def broken(*a, **k):
        raise OSError("disk full")
    monkeypatch.setattr(store, "write_ledger", broken)
    with pytest.raises(OSError):
        record_in_ledger(store, ACC, _row("SRR1_1.fastq.gz"), tmp_path / "staging")
    assert not (tmp_path / "host" / "raw" / ACC / mirror_mod.LEDGER_LOCK).exists()


def test_ssh_ledger_lock_runs_on_the_host(tmp_path, monkeypatch):
    """The remote lock script, run by a local shell standing in for ssh."""
    import subprocess

    real_run = subprocess.run

    def fake_run(args, **kw):
        assert args[:3] == ["ssh", "-o", "BatchMode=yes"]
        return real_run(["sh", "-c", args[-1]], **kw)
    monkeypatch.setattr(mirror_mod.subprocess, "run", fake_run)

    store = SSHStore(f"me@host:{tmp_path}/share")
    lock = tmp_path / "share" / "raw" / ACC / mirror_mod.LEDGER_LOCK
    store.lock_ledger(ACC, wait=5)
    assert lock.is_dir()
    with pytest.raises(MirrorError, match="was held for 1s"):
        store.lock_ledger(ACC, wait=1)
    store.unlock_ledger(ACC)
    assert not lock.exists()


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
                               "--staging", str(staging), "--source", "ena"])
    assert out.exit_code == 0, out.output
    assert "Mirrored 4 file(s)" in out.output and "complete" in out.output

    ena.corrupt = {"SRR1_1.fastq.gz"}
    (host / "raw" / ACC / "SRR1_1.fastq.gz").unlink()
    out = runner.invoke(main, ["mirror-fastq", "--accession", ACC, "--dest", str(host),
                               "--staging", str(staging), "--source", "ena"])
    assert out.exit_code == 1 and "1 failed" in out.output


# --- NCBI source ---------------------------------------------------------------

S3 = "https://sra-pub-run-odp.s3.amazonaws.com/sra"

FAKE_FASTERQ = """#!{python}
import json, sys
from pathlib import Path
args = sys.argv[1:]
out = Path(args[args.index("--outdir") + 1])
sra = Path(args[-1])
out.mkdir(parents=True, exist_ok=True)
for suffix, text in json.loads(sra.read_text()).items():
    (out / f"{{sra.stem}}{{suffix}}.fastq").write_text(text)
"""


class FakeNCBI:
    """Each .sra is JSON {output suffix: FASTQ text}, which the fake fasterq-dump writes out."""

    def __init__(self, ena):
        self.ena = ena
        self.sra = {run: json.dumps({f"_{m}": _fq(run, m).decode() for m in (1, 2)}).encode()
                    for run in ("SRR1", "SRR2")}
        self.missing: set[str] = set()
        self.bad_md5: set[str] = set()

    def locate(self, url):
        run = url.split("acc=")[1].split("&")[0]
        if run in self.missing:
            return {"result": [{"bundle": run, "status": 404, "msg": "no data"}]}
        md5 = "0" * 32 if run in self.bad_md5 else hashlib.md5(self.sra[run]).hexdigest()
        return {"result": [{"files": [{"type": "sra", "size": len(self.sra[run]), "md5": md5,
                                       "locations": [{"service": "s3", "link": f"{S3}/{run}/{run}"}]}]}]}

    def open(self, req):
        if req.full_url.startswith(S3):
            run = req.full_url.rsplit("/", 1)[-1]
            self.ena.requests.append((f"{run}.sra", req.get_header("Range") or ""))
            return _Resp(self.sra[run], 200)
        return self.ena.open(req)


@pytest.fixture
def ncbi(ena, monkeypatch, tmp_path):
    fake = FakeNCBI(ena)
    monkeypatch.setattr(mirror_mod, "_get_json", fake.locate)
    monkeypatch.setattr(mirror_mod, "_urlopen", fake.open)
    kit = tmp_path / "sratoolkit"
    (kit / "bin").mkdir(parents=True)
    tool = kit / "bin" / "fasterq-dump"
    tool.write_text(FAKE_FASTERQ.format(python=sys.executable))
    tool.chmod(0o755)
    monkeypatch.setenv("WHALE_SRATOOLKIT", str(kit))
    monkeypatch.delenv("WHALE_PIGZ", raising=False)
    return fake


def test_ncbi_converts_and_records(repo, ncbi, ena, tmp_path):
    res = _run(repo, tmp_path, source="ncbi")
    assert res["failed"] == [] and res["sources"] == {"ncbi_sra": 4, "ena": 0}
    assert {n for n, _ in ena.requests} == {"SRR1.sra", "SRR2.sra"}
    host = tmp_path / "host/raw" / ACC
    assert gzip.decompress((host / "SRR2_1.fastq.gz").read_bytes()) == _fq("SRR2", 1)
    row = LocalStore(tmp_path / "host").read_ledger(ACC)["SRR2_1.fastq.gz"]
    assert row["source"] == "ncbi_sra" and row["read_count"] == "500"
    assert row["md5"] == hashlib.md5((host / "SRR2_1.fastq.gz").read_bytes()).hexdigest()
    assert row["ena_md5"] == ena.md5("SRR2_1.fastq.gz") != row["md5"]
    assert row["sra_md5"] == hashlib.md5(ncbi.sra["SRR2"]).hexdigest()
    assert row["url"] == f"{S3}/SRR2/SRR2"
    assert not (tmp_path / "staging" / ACC / "SRR2").exists()

    # Converted files are skipped on a rerun although their MD5 is not ENA's.
    ena.requests.clear()
    assert _run(repo, tmp_path, source="ncbi")["mirrored"] == [] and ena.requests == []


def test_ncbi_falls_back_to_ena_per_run(repo, ncbi, ena, tmp_path):
    ncbi.missing = {"SRR1"}
    ncbi.bad_md5 = {"SRR2"}
    res = _run(repo, tmp_path, source="ncbi")
    assert res["failed"] == [] and res["sources"] == {"ncbi_sra": 0, "ena": 4}
    ledger = LocalStore(tmp_path / "host").read_ledger(ACC)
    assert {r["source"] for r in ledger.values()} == {"fastq_ftp"}
    assert ledger["SRR1_1.fastq.gz"]["md5"] == ena.md5("SRR1_1.fastq.gz")


def test_ncbi_read_count_mismatch_falls_back(repo, ncbi, ena, tmp_path):
    short = _fq("SRR1", 2)[: -len(b"@SRR1.2\nACGT\n+\nIIII\n")]
    ncbi.sra["SRR1"] = json.dumps({"_1": _fq("SRR1", 1).decode(), "_2": short.decode()}).encode()
    res = _run(repo, tmp_path, source="ncbi")
    assert res["sources"] == {"ncbi_sra": 2, "ena": 2}
    assert LocalStore(tmp_path / "host").read_ledger(ACC)["SRR1_2.fastq.gz"]["source"] == "fastq_ftp"


def test_ncbi_unpaired_leftover_counts_then_is_dropped(repo, ncbi, ena, tmp_path, monkeypatch):
    lone = "@SRR1.x\nACGT\n+\nIIII\n"
    ncbi.sra["SRR1"] = json.dumps({"_1": _fq("SRR1", 1).decode(), "_2": _fq("SRR1", 2).decode(),
                                   "": lone}).encode()
    report = ena.report

    def with_lone(acc):
        rows = report(acc)
        rows[0].update(read_count="501", base_count="4004")
        return rows
    monkeypatch.setattr(mirror_mod, "fetch_run_report", with_lone)
    res = _run(repo, tmp_path, source="ncbi")
    assert res["sources"] == {"ncbi_sra": 4, "ena": 0}
    assert sorted(p.name for p in (tmp_path / "host/raw" / ACC).glob("SRR1*")) == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]


def test_ncbi_needs_fasterq_dump(repo, ena, tmp_path, monkeypatch):
    monkeypatch.setenv("WHALE_SRATOOLKIT", str(tmp_path / "nowhere"))
    with pytest.raises(MirrorError, match="fasterq-dump not found"):
        _run(repo, tmp_path, source="ncbi")


def test_cli_default_source_is_ncbi(repo, ncbi, tmp_path, monkeypatch):
    monkeypatch.setenv("WHALE_REPO", str(repo))
    out = CliRunner().invoke(main, ["mirror-fastq", "--accession", ACC, "--dest", str(tmp_path / "host"),
                                    "--staging", str(tmp_path / "staging")])
    assert out.exit_code == 0, out.output
    assert "4 converted from NCBI, 0 from ENA" in out.output
