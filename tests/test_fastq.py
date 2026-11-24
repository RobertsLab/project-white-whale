"""Offline tests for `whale fetch-fastq` (src/whale/fastq.py).

A local HTTP server with byte-range support stands in for ENA. No network
access beyond localhost.
"""
from __future__ import annotations

import csv
import hashlib
import http.server
import threading
import urllib.error

import pytest
from click.testing import CliRunner

import whale.fastq as fq
from whale.cli import main
from whale.fastq import (
    FastqError,
    ManifestFile,
    fetch_all,
    fetch_one,
    plan,
    read_manifest,
    record_staged,
)

PAYLOAD = {f"SRR1_{m}.fastq.gz": bytes(range(256)) * (40 + m) for m in (1, 2)}


class _Handler(http.server.BaseHTTPRequestHandler):
    ignore_range = False
    requests: list[str] = []
    routes: dict[str, bytes] = {}   # full path -> body, checked before PAYLOAD by file name

    def do_GET(self):  # noqa: N802
        name = self.path.rsplit("/", 1)[-1]
        type(self).requests.append(f"{self.path if self.path in self.routes else name} "
                                   f"{self.headers.get('Range') or ''}".strip())
        data = self.routes.get(self.path, PAYLOAD.get(name) if not self.path.startswith("/raw/") else None)
        if data is None:
            self.send_error(404)
            return
        rng = self.headers.get("Range")
        if rng and not self.ignore_range:
            start = int(rng.split("=")[1].split("-")[0])
            body = data[start:]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
        else:
            body = data
            self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


@pytest.fixture
def server(monkeypatch):
    # Clusters that route web traffic through a proxy (klone sets http_proxy,
    # with 127.0.0.1 missing from no_proxy) would send these requests to the
    # proxy, which cannot reach the test server.
    for var in ("no_proxy", "NO_PROXY"):
        monkeypatch.setenv(var, "127.0.0.1,localhost")
    _Handler.requests = []
    _Handler.ignore_range = False
    _Handler.routes = {}
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}/vol1/fastq/SRR1", _Handler
    httpd.shutdown()


def _files(base, md5_override=None):
    return [ManifestFile("SRR1", name, f"{base}/{name}",
                         (md5_override or {}).get(name, hashlib.md5(data).hexdigest()), len(data))
            for name, data in PAYLOAD.items()]


def test_download_verify_and_skip(server, tmp_path):
    base, handler = server
    files = _files(base)
    res = fetch_all(files, tmp_path, jobs=2)
    assert res == {"downloaded": 2, "skipped": 0, "failed": {}}
    for name, data in PAYLOAD.items():
        assert (tmp_path / name).read_bytes() == data
        assert (tmp_path / f"{name}.md5").read_text().startswith(hashlib.md5(data).hexdigest())
        assert not (tmp_path / f"{name}.part").exists()
    handler.requests.clear()
    assert fetch_all(files, tmp_path)["skipped"] == 2
    assert handler.requests == []
    assert plan(files, tmp_path)["remaining_bytes"] == 0


def test_resume_from_partial(server, tmp_path):
    base, handler = server
    f = _files(base)[0]
    (tmp_path / f"{f.file}.part").write_bytes(PAYLOAD[f.file][:1000])
    (tmp_path / f"{f.file}.part.src").write_text(f.https_url + "\n")
    assert plan([f], tmp_path)["remaining_bytes"] == f.bytes - 1000
    assert fetch_one(f, tmp_path) == "downloaded"
    assert handler.requests == [f"{f.file} bytes=1000-"]
    assert (tmp_path / f.file).read_bytes() == PAYLOAD[f.file]


def test_server_ignoring_range_restarts_cleanly(server, tmp_path):
    base, handler = server
    handler.ignore_range = True
    f = _files(base)[0]
    (tmp_path / f"{f.file}.part").write_bytes(PAYLOAD[f.file][:1000])
    (tmp_path / f"{f.file}.part.src").write_text(f.https_url + "\n")
    assert fetch_one(f, tmp_path) == "downloaded"
    assert (tmp_path / f.file).read_bytes() == PAYLOAD[f.file]


def test_md5_mismatch_deletes(server, tmp_path):
    base, _ = server
    f = _files(base, {"SRR1_1.fastq.gz": "0" * 32})[0]
    with pytest.raises(FastqError, match="Deleted"):
        fetch_one(f, tmp_path)
    assert not any(tmp_path.iterdir())


def test_unverified_file_in_place_is_checked_not_refetched(server, tmp_path):
    base, handler = server
    f = _files(base)[0]
    (tmp_path / f.file).write_bytes(PAYLOAD[f.file])
    assert fetch_one(f, tmp_path) == "skipped"
    assert handler.requests == [] and (tmp_path / f"{f.file}.md5").exists()
    # A wrong file in place is replaced.
    (tmp_path / f"{f.file}.md5").unlink()
    (tmp_path / f.file).write_bytes(b"x" * f.bytes)
    assert fetch_one(f, tmp_path) == "downloaded"
    assert (tmp_path / f.file).read_bytes() == PAYLOAD[f.file]


def test_retries_then_succeeds_and_gives_up(tmp_path):
    data = b"ACGT" * 100
    f = ManifestFile("SRR9", "SRR9_1.fastq.gz", "ftp://host/SRR9_1.fastq.gz", hashlib.md5(data).hexdigest(), len(data))
    calls = []

    def flaky(url, part, offset):
        calls.append((url, offset))
        if len(calls) == 1:
            with open(part, "wb") as out:
                out.write(data[:150])
            raise urllib.error.URLError("reset")
        with open(part, "ab") as out:
            out.write(data[offset:])

    assert fetch_one(f, tmp_path, fetch=flaky, sleep=lambda s: None) == "downloaded"
    assert calls == [("https://host/SRR9_1.fastq.gz", 0), ("https://host/SRR9_1.fastq.gz", 150)]

    def dead(url, part, offset):
        raise ConnectionError("down")
    g = ManifestFile("SRR9", "SRR9_2.fastq.gz", "ftp://host/SRR9_2.fastq.gz", "a" * 32, 10)
    with pytest.raises(FastqError, match="after 4 attempts"):
        fetch_one(g, tmp_path, fetch=dead, sleep=lambda s: None)


def test_read_manifest(tmp_path):
    path = tmp_path / "m.tsv"
    path.write_text("run_accession\tmate\tfile\turl\tmd5\tbytes\tsource\n"
                    "SRR1\t1\tSRR1_1.fastq.gz\tftp://h/SRR1_1.fastq.gz\tabc\t10\tfastq_ftp\n"
                    "SRR2\t1\tSRR2_1.fastq.gz\tftp://h/SRR2_1.fastq.gz\tdef\t20\tfastq_ftp\n")
    files = read_manifest(path)
    assert [f.bytes for f in files] == [10, 20] and files[0].https_url == "https://h/SRR1_1.fastq.gz"
    assert [f.run_accession for f in read_manifest(path, ["SRR2"])] == ["SRR2"]
    with pytest.raises(FastqError, match="not in the manifest"):
        read_manifest(path, ["SRR3"])
    path.write_text(path.read_text().replace("\tabc\t", "\t\t"))
    with pytest.raises(FastqError, match="no MD5"):
        read_manifest(path)


def test_record_staged(server, tmp_path):
    base, _ = server
    files = _files(base)
    fetch_one(files[0], tmp_path / "raw")
    record = tmp_path / "study" / "fastq_staged.tsv"
    record_staged(files, tmp_path / "raw", record, machine="raven", today="2026-10-23")
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    assert [(r["file"], r["machine"]) for r in rows] == [("SRR1_1.fastq.gz", "raven")]
    fetch_one(files[1], tmp_path / "raw")
    record_staged(files, tmp_path / "raw", record, machine="mac", today="2026-10-24")
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    assert [(r["file"], r["machine"]) for r in rows] == [("SRR1_1.fastq.gz", "raven"), ("SRR1_2.fastq.gz", "mac")]


def test_cli(server, tmp_path, repo_root, monkeypatch):
    base, _ = server
    monkeypatch.setenv("WHALE_REPO", str(tmp_path))
    (tmp_path / "registry" / "controlled_vocabularies").mkdir(parents=True)
    study = tmp_path / "data" / "studies" / "X2020_Y"
    study.mkdir(parents=True)
    with open(study / "fastq_manifest.tsv", "w") as fh:
        fh.write("run_accession\tmate\tfile\turl\tmd5\tbytes\tsource\n")
        for f in _files(base):
            fh.write(f"SRR1\t1\t{f.file}\t{f.url}\t{f.md5}\t{f.bytes}\tfastq_ftp\n")
    monkeypatch.setenv("WHALE_RAW", str(tmp_path / "raw"))

    res = CliRunner().invoke(main, ["fetch-fastq", "--study", "X2020_Y", "--dry-run"])
    assert res.exit_code == 0, res.output
    assert "to download      : 2" in res.output and not (tmp_path / "raw").exists()

    res = CliRunner().invoke(main, ["fetch-fastq", "--study", "X2020_Y"])
    assert res.exit_code == 0, res.output
    assert "2 downloaded" in res.output
    assert (tmp_path / "raw" / "X2020_Y" / "SRR1_2.fastq.gz").exists()
    assert (study / "fastq_staged.tsv").exists()

    monkeypatch.setattr(fq.shutil, "disk_usage", lambda p: type("U", (), {"free": 0})())
    (tmp_path / "raw" / "X2020_Y" / "SRR1_2.fastq.gz.md5").unlink()
    (tmp_path / "raw" / "X2020_Y" / "SRR1_2.fastq.gz").unlink()
    res = CliRunner().invoke(main, ["fetch-fastq", "--study", "X2020_Y"])
    assert res.exit_code != 0 and "not enough free space" in res.output


def test_cli_manifest_outside_repo_keeps_record_beside_it(server, tmp_path, monkeypatch):
    base, _ = server
    repo = tmp_path / "repo"
    (repo / "registry" / "controlled_vocabularies").mkdir(parents=True)
    monkeypatch.setenv("WHALE_REPO", str(repo))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    manifest = elsewhere / "fastq_manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("run_accession\tmate\tfile\turl\tmd5\tbytes\tsource\n")
        for f in _files(base):
            fh.write(f"SRR1\t1\t{f.file}\t{f.url}\t{f.md5}\t{f.bytes}\tfastq_ftp\n")
    res = CliRunner().invoke(main, ["fetch-fastq", "--study", "X2020_Y", "--manifest", str(manifest),
                                    "--dest", str(tmp_path / "raw")])
    assert res.exit_code == 0, res.output
    assert (elsewhere / "fastq_staged.tsv").exists()
    assert not (repo / "data").exists()


def test_second_downloader_into_same_dir_is_refused(server, tmp_path):
    import fcntl
    base, _ = server
    tmp_path.mkdir(exist_ok=True)
    with open(tmp_path / ".whale-fetch.lock", "w") as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(FastqError, match="already downloading"):
            fetch_all(_files(base), tmp_path)
    assert fetch_all(_files(base), tmp_path)["downloaded"] == 2


# --- --launch -----------------------------------------------------------------

LAUNCH_ACC, LAUNCH_SID = "PRJNA999903", "TESTER2099_LAUNCH_WGBS"


def _claimed_worktree(claim_env, server, monkeypatch, profile):
    """A claimed study worktree with a manifest, and a fake `nextflow` on PATH that logs its call."""
    import os
    import stat

    from whale.claims import start_study

    base, _ = server
    wt = claim_env["tmp"] / "wt"
    start_study(LAUNCH_ACC, LAUNCH_SID, repo=claim_env["clone"], profile=profile, worktree_dir=wt,
                registered_by="Tester", machine="hete", open_pr=False)
    with open(wt / "data" / "studies" / LAUNCH_SID / "fastq_manifest.tsv", "w") as fh:
        fh.write("run_accession\tmate\tfile\turl\tmd5\tbytes\tsource\n")
        for f in _files(base):
            fh.write(f"SRR1\t1\t{f.file}\t{f.url}\t{f.md5}\t{f.bytes}\tfastq_ftp\n")
    bin_dir = claim_env["tmp"] / "bin"
    bin_dir.mkdir()
    log = claim_env["tmp"] / "nextflow.log"
    fake = bin_dir / "nextflow"
    fake.write_text(f'#!/bin/sh\necho "$PWD|$WHALE_RAW|$*" > {log}\nexit 3\n')
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("WHALE_REPO", str(wt))
    monkeypatch.setenv("WHALE_RAW", str(claim_env["tmp"] / "raw"))
    return wt, log


def _set_config_chemistry(wt, value):
    cfg = wt / "config" / f"{LAUNCH_SID}.config"
    text = cfg.read_text()
    assert "chemistry    = null" in text
    cfg.write_text(text.replace("chemistry    = null", f"chemistry    = {value}"))


def test_cli_launch_refuses_a_config_without_the_yaml_chemistry(claim_env, server, monkeypatch):
    import yaml

    wt, log = _claimed_worktree(claim_env, server, monkeypatch, "desktop")
    # As start-study writes it: chemistry = null.
    res = CliRunner().invoke(main, ["fetch-fastq", "--source", "ena", "--launch"])
    assert res.exit_code != 0 and "sets no chemistry" in res.output, res.output

    yml = wt / "registry" / "studies" / f"{LAUNCH_SID}.yaml"
    study = yaml.safe_load(yml.read_text())
    study["chemistry"] = "wgbs"
    yml.write_text(yaml.safe_dump(study, sort_keys=False))
    _set_config_chemistry(wt, "'rrbs'")
    res = CliRunner().invoke(main, ["fetch-fastq", "--source", "ena", "--launch"])
    assert res.exit_code != 0 and "the study YAML says 'wgbs'" in res.output, res.output
    assert not (claim_env["tmp"] / "raw").exists() and not log.exists()


def test_cli_launch_stages_then_runs_nextflow(claim_env, server, monkeypatch):
    import yaml

    wt, log = _claimed_worktree(claim_env, server, monkeypatch, "desktop")
    _set_config_chemistry(wt, "'wgbs'")
    res = CliRunner().invoke(main, ["fetch-fastq", "--source", "ena", "--launch", "--dry-run"])
    assert res.exit_code == 0, res.output
    assert "Would then launch: nextflow run" in res.output and not log.exists()

    res = CliRunner().invoke(main, ["fetch-fastq", "--source", "ena", "--launch"])
    # The run's exit status is fetch-fastq's.
    assert res.exit_code == 3, res.output
    status = yaml.safe_load((wt / "data" / "studies" / LAUNCH_SID / "status.yaml").read_text())
    assert status["stage"] == "fastq_staged"
    cwd, raw, args = log.read_text().strip().split("|")
    assert cwd == str(wt.resolve()) and raw == str((claim_env["tmp"] / "raw").resolve())
    assert args == (f"run workflows/methylation -c config/{LAUNCH_SID}.config -profile desktop,docker")


def test_cli_launch_refusals_happen_before_download(claim_env, server, monkeypatch):
    wt, log = _claimed_worktree(claim_env, server, monkeypatch, "klone")
    for extra, msg in ([[], "sbatch head job"], [["--run", "SRR1"], "drop --run"],
                       [["--dest", str(claim_env["tmp"] / "elsewhere")], "folder named"]):
        res = CliRunner().invoke(main, ["fetch-fastq", "--source", "ena", "--launch", *extra])
        assert res.exit_code != 0 and msg in res.output, res.output
    assert not (claim_env["tmp"] / "raw").exists() and not log.exists()
