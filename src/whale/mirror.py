"""Mirror the backlog's FASTQ to the lab data host before anyone claims a study.

Downloading is the slowest step of a study and needs no compute, so it runs
ahead of the claims, a little at a time, from raven. The file list (names,
MD5s, read and base counts) always comes from ENA's read-run report. Each run
is then fetched from one of two sources:

* NCBI (default): the run's .sra from NCBI's AWS open-data bucket, found with
  NCBI's locator service and checked against its MD5, converted here with
  fasterq-dump, gzipped under ENA's file names, and checked against ENA's
  read and base counts. Measured from raven, ~25 MB/s against ENA's ~1 MB/s.
* ENA (fallback, or ``--source ena``): each FASTQ over HTTPS, checked against
  ENA's MD5. Used for a run NCBI has no full-quality copy of, one whose NCBI
  checks fail, or one only partly missing.

Downloads go to a staging directory and resume a partial file. Each finished
file is copied to the data host (data/data_host.yaml) with rsync, its MD5
optionally re-checked there, and recorded. On the data host:

* ``raw/<ACCESSION>/<file>``               every FASTQ, named as ENA names it
* ``raw/<ACCESSION>/mirror_manifest.tsv``  one row per verified file: run, mate,
                                           source URL, MD5 and bytes of the file
                                           as stored, ENA's MD5, the .sra MD5
                                           and read count when converted, date

A converted FASTQ holds the same reads as ENA's but is not byte-identical, so
its MD5 differs from ENA's: check staged copies against the ledger's ``md5``.
ENA's MD5 is kept as ``ena_md5`` so a run ENA reprocesses (new MD5) is fetched
again. For a registered study the runs are those in its committed
fastq_manifest.tsv; otherwise every run with a methylation library strategy.
A file already in the ledger against ENA's current MD5, and present on the
host at its recorded size, is skipped, so the command can be rerun (from
cron) until the backlog is complete.

Nothing here touches a study branch or claim; staging FASTQ for a run still
happens per study (docs/adding_a_study.md).
"""
from __future__ import annotations

import csv
import fcntl
import hashlib
import io
import json
import os
import shlex
import shutil
import subprocess
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from http.client import HTTPException
from pathlib import Path

import yaml

from whale.ena import (
    EXPECTED_TAXID,
    METHYLATION_STRATEGIES,
    ENAError,
    fetch_run_report,
    run_fastq_files,
)
from whale.registry import load_backlog, load_studies

LEDGER_NAME = "mirror_manifest.tsv"
LEDGER_COLUMNS = ["run_accession", "mate", "file", "url", "md5", "bytes", "source", "ena_md5", "sra_md5",
                  "read_count", "mirrored_on"]
NCBI_LOCATE = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
TIMEOUT = 120
CHUNK = 1 << 20
RETRIES = 5
BACKOFF = [10, 30, 60, 120, 300]

_sleep = time.sleep


class MirrorError(RuntimeError):
    """A file could not be mirrored and verified."""


# --- what to mirror ------------------------------------------------------------

def https_url(url: str) -> str:
    """ENA's ftp:// FASTQ URLs are served over HTTPS too, which allows resuming with Range."""
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    if url.startswith(("https://", "http://")):
        return url
    return "https://" + url


def registered_runs(repo: Path, accession: str) -> set[str] | None:
    """Run accessions in the committed fastq_manifest.tsv of the study registered for this accession."""
    for sid, study in load_studies(repo).items():
        if (study.get("accessions") or {}).get("bioproject") != accession:
            continue
        manifest = Path(repo) / "data" / "studies" / sid / "fastq_manifest.tsv"
        if manifest.exists():
            with manifest.open(newline="") as fh:
                return {r["run_accession"] for r in csv.DictReader(fh, delimiter="\t")}
    return None


def wanted_files(repo: Path, accession: str, library_strategies: list[str] | None = None,
                 fetch_runs: Callable[[str], list[dict]] | None = None) -> tuple[list[dict], list[str]]:
    """The FASTQ files to mirror for one accession, with ENA's current MD5s, and notes on what was left out."""
    runs = (fetch_runs or fetch_run_report)(accession)
    notes = []
    keep_runs = registered_runs(repo, accession)
    if keep_runs is not None:
        kept = [r for r in runs if r["run_accession"] in keep_runs]
        gone = keep_runs - {r["run_accession"] for r in kept}
        if gone:
            notes.append(f"{len(gone)} registered run(s) no longer in ENA's report: {sorted(gone)[:5]}")
    else:
        keep = list(library_strategies or METHYLATION_STRATEGIES)
        kept = [r for r in runs if r.get("library_strategy") in keep]
        other: dict[str, int] = {}
        for r in runs:
            if r.get("library_strategy") not in keep:
                other[r.get("library_strategy") or "?"] = other.get(r.get("library_strategy") or "?", 0) + 1
        if other:
            notes.append(f"not mirrored (library strategy): {dict(sorted(other.items()))}")
        wrong_taxon = [r["run_accession"] for r in kept if r.get("tax_id") and r["tax_id"] != EXPECTED_TAXID]
        if wrong_taxon:
            notes.append(f"not mirrored (not taxid {EXPECTED_TAXID}): {wrong_taxon[:5]}")
            kept = [r for r in kept if r["run_accession"] not in wrong_taxon]

    files, no_fastq = [], []
    for r in sorted(kept, key=lambda x: x["run_accession"]):
        source, found, _dropped = run_fastq_files(r)
        if source is None:
            no_fastq.append(r["run_accession"])
            continue
        for idx, (url, md5, size) in enumerate(found):
            files.append({"run_accession": r["run_accession"], "mate": str(idx + 1),
                          "file": url.rsplit("/", 1)[-1], "url": f"ftp://{url}", "md5": md5,
                          "bytes": size, "source": source, "read_count": r.get("read_count") or "",
                          "base_count": r.get("base_count") or ""})
    if no_fastq:
        notes.append(f"{len(no_fastq)} run(s) have no usable FASTQ: {no_fastq[:5]}")
    names = [f["file"] for f in files]
    clashes = sorted({n for n in names if names.count(n) > 1})
    if clashes:
        raise MirrorError(f"{accession}: several runs have files named {clashes[:5]}; "
                          "they cannot share raw/<ACCESSION>/")
    return files, notes


# --- where it goes -------------------------------------------------------------

def _read_ledger_text(text: str) -> dict[str, dict]:
    if not text.strip():
        return {}
    return {r["file"]: r for r in csv.DictReader(io.StringIO(text), delimiter="\t")}


def _ledger_text(rows: dict[str, dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=LEDGER_COLUMNS, delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows[k] for k in sorted(rows))
    return buf.getvalue()


class LocalStore:
    """A directory on this machine (a mounted share, or a test)."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self.label = str(self.root)

    def _dir(self, acc: str) -> Path:
        return self.root / "raw" / acc

    def listing(self, acc: str) -> dict[str, int]:
        d = self._dir(acc)
        return {p.name: p.stat().st_size for p in d.iterdir() if p.is_file()} if d.is_dir() else {}

    def read_ledger(self, acc: str) -> dict[str, dict]:
        p = self._dir(acc) / LEDGER_NAME
        return _read_ledger_text(p.read_text()) if p.exists() else {}

    def write_ledger(self, acc: str, rows: dict[str, dict], staging: Path) -> None:
        d = self._dir(acc)
        d.mkdir(parents=True, exist_ok=True)
        tmp = d / (LEDGER_NAME + ".tmp")
        tmp.write_text(_ledger_text(rows))
        tmp.replace(d / LEDGER_NAME)

    def put(self, src: Path, acc: str, name: str) -> None:
        d = self._dir(acc)
        d.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(d / name))

    def md5(self, acc: str, name: str) -> str:
        return md5_file(self._dir(acc) / name)


class SSHStore:
    """``user@host:/path`` reached over SSH (key auth) and rsync, as gannet is from raven."""

    def __init__(self, target: str):
        host, sep, path = target.partition(":")
        if not sep or not host or not path:
            raise MirrorError(f"upload target must look like user@host:/path, not {target!r}")
        self.host, self.path, self.label = host, path.rstrip("/"), target

    def _dir(self, acc: str) -> str:
        return f"{self.path}/raw/{acc}"

    def _ssh(self, command: str) -> str:
        proc = subprocess.run(["ssh", "-o", "BatchMode=yes", self.host, command],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise MirrorError(f"ssh {self.host} failed: {proc.stderr.strip() or proc.returncode}")
        return proc.stdout

    def _rsync(self, src: Path, dest: str) -> None:
        proc = subprocess.run(
            ["rsync", "-tp", "--partial", "--chmod=D755,F644", "-e", "ssh -o BatchMode=yes",
             str(src), f"{self.host}:{dest}"], capture_output=True, text=True)
        if proc.returncode != 0:
            raise MirrorError(f"rsync to {self.host}:{dest} failed: {proc.stderr.strip()}")

    def listing(self, acc: str) -> dict[str, int]:
        d = shlex.quote(self._dir(acc))
        out = self._ssh(f"if [ -d {d} ]; then find {d} -maxdepth 1 -type f -printf '%f\\t%s\\n'; fi")
        return {name: int(size) for name, size in (ln.split("\t") for ln in out.splitlines() if ln)}

    def read_ledger(self, acc: str) -> dict[str, dict]:
        p = shlex.quote(f"{self._dir(acc)}/{LEDGER_NAME}")
        return _read_ledger_text(self._ssh(f"if [ -f {p} ]; then cat {p}; fi"))

    def write_ledger(self, acc: str, rows: dict[str, dict], staging: Path) -> None:
        local = Path(staging) / acc / LEDGER_NAME
        local.parent.mkdir(parents=True, exist_ok=True)
        local.write_text(_ledger_text(rows))
        self._ssh(f"mkdir -p {shlex.quote(self._dir(acc))}")
        self._rsync(local, f"{self._dir(acc)}/{LEDGER_NAME}")

    def put(self, src: Path, acc: str, name: str) -> None:
        self._ssh(f"mkdir -p {shlex.quote(self._dir(acc))}")
        self._rsync(src, f"{self._dir(acc)}/{name}")
        src.unlink()

    def md5(self, acc: str, name: str) -> str:
        return self._ssh(f"md5sum {shlex.quote(f'{self._dir(acc)}/{name}')}").split()[0]


def load_data_host(repo: Path) -> dict:
    return yaml.safe_load((Path(repo) / "data" / "data_host.yaml").read_text())


def default_staging(repo: Path) -> Path:
    if os.environ.get("WHALE_MIRROR_STAGING"):
        return Path(os.environ["WHALE_MIRROR_STAGING"])
    raw = os.environ.get("WHALE_RAW")
    return (Path(raw) if raw else Path(repo) / "data" / "raw") / ".mirror-staging"


# --- downloading ---------------------------------------------------------------

def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


class Throttle:
    """A total download rate shared by every worker thread (bytes/s; None = unlimited)."""

    def __init__(self, rate: float | None):
        self.rate = rate
        self._lock = threading.Lock()
        self._next = time.monotonic()

    def take(self, n: int) -> None:
        if not self.rate:
            return
        with self._lock:
            now = time.monotonic()
            self._next = max(self._next, now) + n / self.rate
            wait = self._next - now
        if wait > 0:
            _sleep(wait)


def _urlopen(req):
    return urllib.request.urlopen(req, timeout=TIMEOUT)  # noqa: S310 - ENA URLs from ENA's report


def download(url: str, dest: Path, *, expected_bytes: int | None = None,
             throttle: Throttle | None = None, opener=None) -> str:
    """Fetch url to dest via dest.part, resuming a partial file; return the MD5 of what was written.

    Network errors are retried with back-off, each time resuming where the
    partial file ends, and so is a response that ends early without an error
    (when the size is known). A server that ignores the Range header restarts
    the file.
    """
    part = dest.with_name(dest.name + ".part")
    if part.exists() and expected_bytes is not None and part.stat().st_size > expected_bytes:
        part.unlink()
    for attempt in range(RETRIES + 1):
        h = hashlib.md5()
        have = 0
        if part.exists():
            with open(part, "rb") as fh:
                for chunk in iter(lambda: fh.read(CHUNK), b""):
                    h.update(chunk)
                    have += len(chunk)
        if expected_bytes is not None and have == expected_bytes:
            break
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with (opener or _urlopen)(req) as resp:
                status = getattr(resp, "status", 200)
                if have and status != 206:
                    have, h = 0, hashlib.md5()
                with open(part, "ab" if have else "wb") as out:
                    for chunk in iter(lambda: resp.read(CHUNK), b""):
                        out.write(chunk)
                        h.update(chunk)
                        have += len(chunk)
                        if throttle:
                            throttle.take(len(chunk))
            if expected_bytes is not None and have < expected_bytes:
                # Retried below like a reset, resuming from here rather than from zero.
                raise ConnectionError(f"response ended at {have:,} of {expected_bytes:,} bytes")
            break
        except (urllib.error.URLError, HTTPException, OSError) as exc:
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (403, 404, 410):
                raise MirrorError(f"{url}: HTTP {exc.code}") from None
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 416:
                part.unlink(missing_ok=True)
            if attempt == RETRIES:
                raise MirrorError(f"{url}: gave up after {RETRIES + 1} attempts ({exc})") from None
            _sleep(BACKOFF[min(attempt, len(BACKOFF) - 1)])
    part.replace(dest)
    return h.hexdigest()


# --- the mirror ----------------------------------------------------------------

class _Lock:
    """One mirror at a time per staging directory, so overlapping cron runs exit quietly."""

    def __init__(self, path: Path):
        self.path = path

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Append mode, so a run that finds the lock taken leaves the holder's PID in the file.
        self.fh = open(self.path, "a")  # noqa: SIM115 - held for the run
        try:
            fcntl.flock(self.fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.fh.close()
            return False
        self.fh.truncate(0)
        self.fh.write(f"{os.getpid()}\n")
        self.fh.flush()
        return True

    def __exit__(self, *exc):
        if not self.fh.closed:
            fcntl.flock(self.fh, fcntl.LOCK_UN)
            self.fh.close()


def _ena_md5(done: dict) -> str:
    """The ENA MD5 a ledger row was mirrored against (rows converted from NCBI hold their own MD5 in md5)."""
    return done.get("ena_md5") or done.get("md5", "")


def plan_mirror(repo: Path, accessions: list[str], store, *, library_strategies=None,
                fetch_runs=None) -> list[dict]:
    """Per accession: the files wanted, which are already mirrored and verified, and what is left."""
    plans = []
    for acc in accessions:
        try:
            files, notes = wanted_files(repo, acc, library_strategies, fetch_runs)
        except (ENAError, MirrorError, urllib.error.URLError, OSError) as exc:
            plans.append({"accession": acc, "error": str(exc), "files": [], "todo": [], "notes": [],
                          "ledger": {}})
            continue
        ledger = store.read_ledger(acc)
        present = store.listing(acc)
        todo = []
        for f in files:
            done = ledger.get(f["file"])
            if done and _ena_md5(done) == f["md5"] and f["file"] in present \
                    and str(present[f["file"]]) == str(done.get("bytes") or present[f["file"]]):
                continue
            if done and _ena_md5(done) != f["md5"]:
                notes.append(f"{f['file']}: ENA's MD5 changed since it was mirrored; fetching again")
            todo.append(f)
        plans.append({"accession": acc, "files": files, "todo": todo, "notes": notes, "ledger": ledger})
    return plans


# --- NCBI: .sra from NCBI's AWS open-data bucket, converted on this machine -----------

class NcbiUnavailable(MirrorError):
    """The run cannot be taken from NCBI (not on S3, check failed); ENA is used instead."""


def _get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:  # noqa: S310 - fixed NCBI host
        return json.loads(fh.read().decode())


def locate_sra(run: str) -> dict:
    """The run's full-quality .sra on S3, with NCBI's MD5 and size, from NCBI's locator service."""
    try:
        doc = _get_json(f"{NCBI_LOCATE}?acc={run}&location=s3.us-east-1&filetype=run")
    except (urllib.error.URLError, HTTPException, OSError, ValueError) as exc:
        raise NcbiUnavailable(f"NCBI locator failed ({exc})") from None
    for res in doc.get("result", []):
        for f in res.get("files", []):
            if f.get("type") != "sra" or not f.get("md5"):
                continue
            for loc in f.get("locations", []):
                if loc.get("service") == "s3" and loc.get("link"):
                    return {"url": loc["link"], "md5": f["md5"], "bytes": int(f["size"])}
    raise NcbiUnavailable("no full-quality .sra on NCBI's S3 bucket")


def find_tools() -> dict:
    """fasterq-dump ($WHALE_SRATOOLKIT/bin, else PATH) and a gzip compressor (pigz if present)."""
    kit = os.environ.get("WHALE_SRATOOLKIT")
    fasterq = str(Path(kit) / "bin" / "fasterq-dump") if kit else shutil.which("fasterq-dump")
    if not fasterq or not Path(fasterq).exists():
        raise MirrorError("fasterq-dump not found: set WHALE_SRATOOLKIT to an sratoolkit directory "
                          "(on raven /home/shared/sratoolkit.3.0.2-ubuntu64), or use --source ena")
    return {"fasterq_dump": fasterq, "pigz": os.environ.get("WHALE_PIGZ") or shutil.which("pigz")}


def _count_fastq(src: Path, dest: Path | None, compress_cmd: list[str] | None) -> tuple[int, int]:
    """(reads, bases) of a FASTQ, gzipping it to dest on the way when dest is given."""
    proc = out = None
    if dest is not None:
        out = open(dest, "wb")  # noqa: SIM115 - closed below
        proc = subprocess.Popen(compress_cmd, stdin=subprocess.PIPE, stdout=out)
    lines = bases = 0
    try:
        with open(src, "rb") as fh:
            while True:
                batch = fh.readlines(1 << 22)
                if not batch:
                    break
                for k in range((1 - lines) % 4, len(batch), 4):
                    bases += len(batch[k].rstrip(b"\r\n"))
                lines += len(batch)
                if proc:
                    proc.stdin.write(b"".join(batch))
    finally:
        if proc:
            proc.stdin.close()
            rc = proc.wait()
            out.close()
            if rc != 0:
                raise NcbiUnavailable(f"{' '.join(compress_cmd)} exited {rc}")
    if lines % 4:
        raise NcbiUnavailable(f"{src.name}: {lines} lines is not whole FASTQ records")
    return lines // 4, bases


def convert_sra(sra: Path, outdir: Path, files: list[dict], *, tools: dict, threads: int,
                tmp: Path) -> list[dict]:
    """fasterq-dump --split-3 the .sra, gzip each mate under ENA's file name, and check the counts.

    Reads and bases must equal ENA's read_count / base_count for the run (an
    unpaired leftover file counts towards them and is then discarded, as the
    ENA path drops it). Returns [{file, path, reads}] in mate order.
    """
    run = sra.stem
    tmp.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([tools["fasterq_dump"], "--split-3", "--threads", str(threads), "--temp", str(tmp),
                           "--outdir", str(outdir), "--force", str(sra)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise NcbiUnavailable(f"fasterq-dump exited {proc.returncode}: {proc.stderr.strip()[-300:]}")
    paired = len(files) == 2
    mates = [outdir / f"{run}_{i + 1}.fastq" for i in range(2)] if paired else [outdir / f"{run}.fastq"]
    if not all(m.exists() for m in mates):
        raise NcbiUnavailable(f"fasterq-dump did not write {[m.name for m in mates]}")
    compress = ([tools["pigz"], "-p", str(max(1, threads // len(mates))), "-c"] if tools.get("pigz")
                else ["gzip", "-c"])
    with ThreadPoolExecutor(max_workers=len(mates)) as pool:
        counts = list(pool.map(lambda mf: _count_fastq(mf[0], outdir / mf[1]["file"], compress),
                               zip(mates, files, strict=True)))
    unpaired = outdir / f"{run}.fastq"
    extra = _count_fastq(unpaired, None, None) if paired and unpaired.exists() else (0, 0)
    for m in [*mates, unpaired]:
        m.unlink(missing_ok=True)

    spots = counts[0][0] + extra[0]
    bases = sum(b for _, b in counts) + extra[1]
    want_reads, want_bases = files[0].get("read_count"), files[0].get("base_count")
    if paired and counts[0][0] != counts[1][0]:
        raise NcbiUnavailable(f"mates differ: {counts[0][0]} vs {counts[1][0]} reads")
    if want_reads and spots != int(want_reads):
        raise NcbiUnavailable(f"{spots} spots, ENA reports {want_reads}")
    if want_bases and bases != int(want_bases):
        raise NcbiUnavailable(f"{bases} bases, ENA reports {want_bases}")
    return [{"file": f["file"], "path": outdir / f["file"], "reads": c[0]}
            for f, c in zip(files, counts, strict=True)]


# --- the mirror ----------------------------------------------------------------

def mirror(repo: Path, accessions: list[str], store, *, staging: Path, jobs: int = 4,
           rate_mb: float | None = None, max_gb: float | None = None, verify_remote: bool = True,
           source: str = "ncbi", threads: int = 4, tmp: Path | None = None,
           library_strategies=None, fetch_runs=None, opener=None,
           log=print, today: str | None = None) -> dict:
    """Download, verify and upload what the plan says is missing; record each file in the ledger.

    Work is done per run. With source "ncbi", a run whose files are all missing
    is taken as one .sra from NCBI's S3 bucket (MD5-checked), converted with
    fasterq-dump and checked against ENA's read and base counts; if any of that
    fails, or only some of a run's files are missing, its files come from ENA.

    Returns {"locked": bool, "plans", "mirrored": [...], "failed": [(what, reason)], "bytes": int,
    "sources": {"ncbi_sra": n, "ena": n}}.
    """
    staging = Path(staging)
    tmp = Path(tmp) if tmp else staging / "tmp"
    today = today or date.today().isoformat()
    result = {"locked": False, "plans": [], "mirrored": [], "failed": [], "bytes": 0,
              "sources": {"ncbi_sra": 0, "ena": 0}}
    tools = find_tools() if source == "ncbi" else None
    with _Lock(staging / "mirror.lock") as got:
        if not got:
            result["locked"] = True
            return result
        plans = plan_mirror(repo, accessions, store, library_strategies=library_strategies,
                            fetch_runs=fetch_runs)
        result["plans"] = plans
        for p in plans:
            if p.get("error"):
                result["failed"].append((p["accession"], p["error"]))

        tasks, budget = [], (max_gb * 1e9 if max_gb else None)
        for p in plans:
            by_run: dict[str, list[dict]] = {}
            for f in p["files"]:
                by_run.setdefault(f["run_accession"], []).append(f)
            todo_runs: dict[str, list[dict]] = {}
            for f in p["todo"]:
                todo_runs.setdefault(f["run_accession"], []).append(f)
            for run, fs in todo_runs.items():
                if budget is not None and budget <= 0:
                    break
                tasks.append((p["accession"], run, fs, len(fs) == len(by_run[run])))
                if budget is not None:
                    budget -= sum(int(f["bytes"] or 0) for f in fs)

        ledgers = {p["accession"]: dict(p["ledger"]) for p in plans}
        ledger_lock = threading.Lock()
        throttle = Throttle(rate_mb * 1e6 if rate_mb else None)

        def publish(acc: str, row: dict, local: Path) -> int:
            size = local.stat().st_size
            store.put(local, acc, row["file"])
            if verify_remote:
                remote = store.md5(acc, row["file"])
                if remote != row["md5"]:
                    raise MirrorError(f"MD5 on {store.label} is {remote}, not {row['md5']}")
            with ledger_lock:
                ledgers[acc][row["file"]] = {**row, "bytes": str(size), "mirrored_on": today}
                store.write_ledger(acc, ledgers[acc], staging)
            return size

        def from_ena(acc: str, f: dict) -> int:
            if not f["md5"]:
                raise MirrorError(f"{f['file']}: ENA lists no MD5 for it, so it cannot be verified")
            local = staging / acc / f["file"]
            local.parent.mkdir(parents=True, exist_ok=True)
            expected = int(f["bytes"]) if f["bytes"] else None
            got_md5 = download(https_url(f["url"]), local, expected_bytes=expected, throttle=throttle,
                               opener=opener)
            if got_md5 != f["md5"]:
                local.unlink()
                got_md5 = download(https_url(f["url"]), local, expected_bytes=expected,
                                   throttle=throttle, opener=opener)
                if got_md5 != f["md5"]:
                    local.unlink()
                    raise MirrorError(f"{f['file']}: MD5 {got_md5} does not match ENA's {f['md5']} (twice)")
            return publish(acc, {**f, "ena_md5": f["md5"]}, local)

        def from_ncbi(acc: str, run: str, fs: list[dict]) -> list[tuple[str, int]]:
            info = locate_sra(run)
            work = staging / acc / run
            work.mkdir(parents=True, exist_ok=True)
            sra = work / f"{run}.sra"
            try:
                for attempt in (1, 2):
                    try:
                        got_md5 = download(info["url"], sra, expected_bytes=info["bytes"], throttle=throttle,
                                           opener=opener)
                    except MirrorError as exc:
                        raise NcbiUnavailable(str(exc)) from None
                    if got_md5 == info["md5"]:
                        break
                    sra.unlink()
                    if attempt == 2:
                        raise NcbiUnavailable(f".sra MD5 {got_md5} does not match NCBI's {info['md5']} (twice)")
                outputs = convert_sra(sra, work, fs, tools=tools, threads=threads, tmp=tmp / run)
                sra.unlink()
                done = []
                for f, out in zip(fs, outputs, strict=True):
                    row = {**f, "url": info["url"], "md5": md5_file(out["path"]), "source": "ncbi_sra",
                           "ena_md5": f["md5"], "sra_md5": info["md5"], "read_count": str(out["reads"])}
                    done.append((f["file"], publish(acc, row, out["path"])))
                return done
            finally:
                shutil.rmtree(work, ignore_errors=True)
                shutil.rmtree(tmp / run, ignore_errors=True)

        def one_run(acc: str, run: str, fs: list[dict], whole: bool) -> list[tuple[str, int, str]]:
            if source == "ncbi" and whole:
                try:
                    return [(name, size, "ncbi_sra") for name, size in from_ncbi(acc, run, fs)]
                except NcbiUnavailable as exc:
                    log(f"note {acc}/{run}: NCBI not used ({exc}); fetching from ENA")
            return [(f["file"], from_ena(acc, f), "ena") for f in fs]

        with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
            futures = {pool.submit(one_run, *t): t for t in tasks}
            for fut in as_completed(futures):
                acc, run = futures[fut][:2]
                try:
                    done = fut.result()
                except (MirrorError, OSError, subprocess.SubprocessError) as exc:
                    result["failed"].append((f"{acc}/{run}", str(exc)))
                    log(f"FAIL {acc}/{run}: {exc}")
                    continue
                for name, size, src in done:
                    result["mirrored"].append(f"{acc}/{name}")
                    result["bytes"] += size
                    result["sources"][src] += 1
                    log(f"ok   {acc}/{name} ({size / 1e9:.2f} GB, {'NCBI .sra' if src == 'ncbi_sra' else 'ENA'})")
    return result


def backlog_accessions(repo: Path) -> list[str]:
    return [b["accession"] for b in load_backlog(repo)]
