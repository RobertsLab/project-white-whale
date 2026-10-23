"""Download a study's FASTQ from its manifest, resuming and MD5-verifying every file.

Reads data/studies/<STUDY_ID>/fastq_manifest.tsv (written by
`whale fetch-samplesheet`) and stages each file under
``$WHALE_RAW/<STUDY_ID>/`` (default ``data/raw/<STUDY_ID>/``), the directory
the workflow reads.

* Files already mirrored on gannet (`whale mirror-fastq`, docs/fastq_mirror.md)
  are copied from there by default, much faster than ENA, and checked against
  the mirror ledger's MD5: a copy converted from NCBI's .sra holds the same
  reads as ENA's file but not the same bytes. A mirrored file is used only if
  the ledger mirrored it against the manifest's current ENA MD5. A file not
  on the mirror, or whose mirror copy fails, comes from ENA (ftp:// URLs are
  fetched over HTTPS from the same host).
* A file is written to ``<name>.part`` and resumed from its size with an HTTP
  range request, so an interrupted download continues rather than restarts.
* A file is complete only when its MD5 matches the manifest. It is then
  renamed into place and a ``<name>.md5`` sidecar is written, so a re-run
  skips it without re-hashing gigabytes. A mismatch deletes the file.
* Free space is checked against what is still to download before starting.

Each verified file is recorded in data/studies/<STUDY_ID>/fastq_staged.tsv
(file, bytes, md5, source, ENA's MD5, machine, date), which is committed on
the study branch.
"""
from __future__ import annotations

import csv
import fcntl
import hashlib
import shutil
import socket
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date
from pathlib import Path

TIMEOUT = 120
CHUNK = 1 << 20
RETRIES = 4


class FastqError(RuntimeError):
    pass


@dataclass
class ManifestFile:
    run_accession: str
    file: str
    url: str
    md5: str
    bytes: int
    source: str = "ena"                       # "ena" or "mirror"
    ena_md5: str = ""                         # the manifest's MD5 (same as md5 for ENA files)
    fallback: ManifestFile | None = None      # the ENA file to use if a mirror copy fails

    def accepted(self) -> dict[str, int]:
        """MD5 -> size of every copy that counts as this file (mirror and ENA)."""
        out = {self.md5: self.bytes}
        if self.fallback is not None:
            out[self.fallback.md5] = self.fallback.bytes
        return out

    @property
    def https_url(self) -> str:
        return "https://" + self.url[len("ftp://"):] if self.url.startswith("ftp://") else self.url


def read_manifest(path: Path, runs: list[str] | None = None) -> list[ManifestFile]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise FastqError(f"{path} lists no files")
    out = []
    for r in rows:
        if not r.get("md5") or not r.get("bytes"):
            raise FastqError(f"{path}: {r.get('file')} has no MD5 or size; regenerate with whale fetch-samplesheet")
        out.append(ManifestFile(r["run_accession"], r["file"], r["url"], r["md5"], int(r["bytes"]),
                                ena_md5=r["md5"]))
    if runs:
        unknown = sorted(set(runs) - {f.run_accession for f in out})
        if unknown:
            raise FastqError(f"run(s) not in the manifest: {unknown}")
        out = [f for f in out if f.run_accession in runs]
    names = [f.file for f in out]
    if len(set(names)) != len(names):
        raise FastqError(f"{path}: file names repeat, so they would overwrite each other")
    return out


def md5_file(path: Path) -> str:
    h = hashlib.md5()  # noqa: S324 - ENA publishes MD5; integrity, not security
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _sidecar(dest: Path, f: ManifestFile) -> tuple[str, str] | None:
    """(md5, source) recorded when the file was verified, or None."""
    path = dest / f"{f.file}.md5"
    if not path.exists():
        return None
    parts = path.read_text().split()
    return (parts[0], parts[2] if len(parts) > 2 else "ena") if parts else None


def is_staged(dest: Path, f: ManifestFile) -> bool:
    """Verified earlier, from the mirror or from ENA: the sidecar's MD5 is one of
    the file's accepted copies and the file has that copy's size."""
    side, path = _sidecar(dest, f), dest / f.file
    accepted = f.accepted()
    return (side is not None and side[0] in accepted and path.exists()
            and path.stat().st_size == accepted[side[0]])


def load_ledger(url: str, get=None) -> dict[str, dict] | None:
    """The mirror ledger at `url` (file -> row), or None if there is none (HTTP 404)."""
    from whale.mirror import _read_ledger_text

    get = get or _http_get_text
    try:
        return _read_ledger_text(get(url))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise


def _http_get_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:  # noqa: S310 - data host from data_host.yaml
        return resp.read().decode()


def apply_mirror(files: list[ManifestFile], ledger: dict[str, dict] | None, base_url: str) -> tuple[list[ManifestFile], dict]:
    """Point each file at its mirror copy when the ledger has it against the manifest's ENA MD5.

    Returns (files, counts): counts of files taken from the mirror, files the
    mirror lacks, and files whose mirror copy is stale (mirrored against an
    older ENA MD5, so ENA has reprocessed the run since).
    """
    from whale.mirror import _ena_md5

    counts = {"mirror": 0, "not_mirrored": 0, "stale": 0}
    out = []
    for f in files:
        row = (ledger or {}).get(f.file)
        if row is None:
            counts["not_mirrored"] += 1
            out.append(f)
        elif _ena_md5(row) != f.md5:
            counts["stale"] += 1
            out.append(f)
        else:
            counts["mirror"] += 1
            out.append(ManifestFile(f.run_accession, f.file, f"{base_url.rstrip('/')}/{f.file}", row["md5"],
                                    int(row["bytes"]), source="mirror", ena_md5=f.md5, fallback=f))
    return out, counts


def plan(files: list[ManifestFile], dest: Path) -> dict:
    todo = [f for f in files if not is_staged(dest, f)]
    remaining = 0
    for f in todo:
        part = dest / f"{f.file}.part"
        remaining += f.bytes - (part.stat().st_size if part.exists() and part.stat().st_size <= f.bytes else 0)
    probe = dest if dest.exists() else next(p for p in dest.parents if p.exists())
    return {"files": len(files), "staged": len(files) - len(todo), "todo": todo,
            "total_bytes": sum(f.bytes for f in files), "remaining_bytes": remaining,
            "free_bytes": shutil.disk_usage(probe).free}


def _http_fetch(url: str, part: Path, offset: int) -> None:
    """Append url's bytes from `offset` to `part`."""
    req = urllib.request.Request(url, headers={"Range": f"bytes={offset}-"} if offset else {})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:  # noqa: S310 - manifest URLs from ENA
        if offset and resp.status != 206:
            # Server ignored the range: start over rather than append a second copy.
            part.unlink(missing_ok=True)
            offset = 0
        with open(part, "ab" if offset else "wb") as out:
            shutil.copyfileobj(resp, out, length=CHUNK)


def _write_sidecar(dest: Path, f: ManifestFile) -> None:
    (dest / f"{f.file}.md5").write_text(f"{f.md5}  {f.file}  {f.source}\n")


def _fetch_from(f: ManifestFile, dest: Path, *, fetch, retries: int, sleep, log) -> str:
    """Download, resume and verify one copy of a file (mirror or ENA)."""
    final, part, part_src = dest / f.file, dest / f"{f.file}.part", dest / f"{f.file}.part.src"
    # A partial file is only resumed from the URL it was started from: a mirror
    # copy and ENA's file are different bytes.
    if part.exists() and (not part_src.exists() or part_src.read_text().strip() != f.https_url):
        part.unlink()
    part_src.write_text(f.https_url + "\n")
    for attempt in range(1, retries + 1):
        offset = part.stat().st_size if part.exists() else 0
        if offset > f.bytes:
            part.unlink()
            offset = 0
        try:
            if offset < f.bytes:
                log(f"{'resume' if offset else 'get   '} {f.file} from {f.source}, {offset:,} of {f.bytes:,} bytes")
                fetch(f.https_url, part, offset)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            if attempt == retries:
                raise FastqError(f"{f.file}: download from {f.source} failed after {retries} attempts: {exc}") from None
            log(f"retry  {f.file} ({exc}); attempt {attempt + 1} of {retries}")
            sleep(min(60, 5 * 2 ** (attempt - 1)))
            continue
        size = part.stat().st_size if part.exists() else 0
        if size < f.bytes:
            if attempt == retries:
                raise FastqError(f"{f.file}: stopped at {size:,} of {f.bytes:,} bytes after {retries} attempts")
            continue
        got = md5_file(part)
        if size != f.bytes or got != f.md5:
            part.unlink()
            part_src.unlink(missing_ok=True)
            raise FastqError(f"{f.file}: {size:,} bytes, MD5 {got}; {f.source} says {f.bytes:,} bytes, "
                             f"MD5 {f.md5}. Deleted; re-run to fetch it again.")
        part.replace(final)
        part_src.unlink(missing_ok=True)
        _write_sidecar(dest, f)
        return "downloaded"
    raise FastqError(f"{f.file}: not downloaded")  # unreachable


def fetch_one(f: ManifestFile, dest: Path, *, fetch=_http_fetch, retries: int = RETRIES,
              sleep=time.sleep, log=lambda m: None) -> str:
    """Download, resume and verify one file, falling back to ENA if a mirror copy fails.

    Returns 'skipped' or 'downloaded'.
    """
    if is_staged(dest, f):
        return "skipped"
    dest.mkdir(parents=True, exist_ok=True)
    final = dest / f.file
    if final.exists():
        # Present but never verified here (copied in by hand, or a crash between rename and sidecar).
        size, accepted = final.stat().st_size, f.accepted()
        if size in accepted.values():
            got = md5_file(final)
            if accepted.get(got) == size:
                copy = f if got == f.md5 else f.fallback
                _write_sidecar(dest, copy)
                return "skipped"
        final.unlink()
    try:
        return _fetch_from(f, dest, fetch=fetch, retries=retries, sleep=sleep, log=log)
    except FastqError as exc:
        if f.fallback is None:
            raise
        log(f"mirror failed for {f.file} ({exc}); falling back to ENA")
        return _fetch_from(f.fallback, dest, fetch=fetch, retries=retries, sleep=sleep, log=log)


RECORD_COLUMNS = ["run_accession", "file", "bytes", "md5", "source", "ena_md5", "machine", "date_verified"]


def record_staged(files: list[ManifestFile], dest: Path, record: Path, machine: str | None = None,
                  today: str | None = None) -> None:
    """Write (or update) the committed record of verified files, one row per file."""
    machine = machine or socket.gethostname().split(".")[0]
    today = today or date.today().isoformat()
    rows: dict[str, dict] = {}
    if record.exists():
        with open(record, newline="") as fh:
            rows = {r["file"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for f in files:
        if not is_staged(dest, f):
            continue
        md5, source = _sidecar(dest, f)
        if rows.get(f.file, {}).get("md5") != md5:
            rows[f.file] = {"run_accession": f.run_accession, "file": f.file, "bytes": str(f.accepted()[md5]),
                            "md5": md5, "source": source, "ena_md5": f.ena_md5 or f.md5,
                            "machine": machine, "date_verified": today}
    record.parent.mkdir(parents=True, exist_ok=True)
    with open(record, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=RECORD_COLUMNS, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(sorted(rows.values(), key=lambda r: r["file"]))


def fetch_all(files: list[ManifestFile], dest: Path, *, jobs: int = 2, fetch=_http_fetch,
              sleep=time.sleep, log=lambda m: None) -> dict:
    """Fetch every file not yet staged, `jobs` at a time. Returns counts and failures.

    Holds an exclusive lock on <dest>/.whale-fetch.lock for the duration: two
    downloaders appending to the same .part files would corrupt them.
    """
    dest.mkdir(parents=True, exist_ok=True)
    result = {"downloaded": 0, "skipped": 0, "failed": {}}
    lock = open(dest / ".whale-fetch.lock", "w")  # noqa: SIM115 - held until the function returns
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise FastqError(f"another whale fetch-fastq is already downloading into {dest}") from None
    with lock, ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {pool.submit(fetch_one, f, dest, fetch=fetch, sleep=sleep, log=log): f for f in files}
        for fut in as_completed(futures):
            f = futures[fut]
            try:
                outcome = fut.result()
                result[outcome] += 1
                if outcome == "downloaded":
                    side = _sidecar(dest, f)
                    log(f"ok     {f.file} (MD5 verified, from {side[1] if side else f.source})")
            except FastqError as exc:
                result["failed"][f.file] = str(exc)
                log(f"FAILED {exc}")
    return result
