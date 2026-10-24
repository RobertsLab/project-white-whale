"""Record a finished pipeline run on its study branch: `whale record-run`.

The methylation workflow calls this when a run succeeds (workflow.onComplete
in workflows/methylation/main.nf), from the checkout it was launched in, so no
stage, summary or finding is ever typed in by hand. It:

1. checks the run is complete: run_manifest.json, qc_summary.tsv and one CpG
   file for every library in the QC table;
2. copies qc_summary.tsv and run_manifest.json into the study's data folder
   (data/studies/<STUDY_ID>/, or parts/<PART>/ for a part);
3. writes methylome_manifest.tsv (every per-sample file's size, SHA-256 and
   URL on the data host) and uploads those files to the data host
   (data/data_host.yaml `upload_target`, methylomes/<STUDY_ID>/<sample>/);
4. writes study_summary.yaml: QC medians, CpG methylation of the libraries
   that passed, and the same split by each design factor and condition;
5. moves status.yaml on through aligned, called, qc_passed and (once the
   methylomes are on the data host) summarized, then commits and pushes;
6. at summarized, marks the study PR ready for review.

The dashboard turns study_summary.yaml into the study's finding. Running it
again is safe: it rewrites the same files and only moves stages forward, so
a failed upload is finished by running it again once the host is reachable.
"""
from __future__ import annotations

import csv
import json
import shlex
import shutil
import statistics
import subprocess
from datetime import date
from pathlib import Path

import yaml

from whale.checksums import sha256_file
from whale.claims import (
    ClaimError,
    _run,
    advance_status,
    current_branch,
    dump_status,
    machine_name,
    parse_branch,
)
from whale.paths import SUMMARY_NAME, status_path, studies_dir
from whale.vocab import stage_order

CPG_SUFFIX = ".CpG_report.merged_CpG_evidence.cov.gz"
METHYLOME_MANIFEST = "methylome_manifest.tsv"
MANIFEST_COLUMNS = ["sample_id", "file", "bytes", "sha256", "url"]
RECORDED_STAGES = ["aligned", "called", "qc_passed", "summarized"]


class RecordError(RuntimeError):
    pass


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _float(v) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _stats(values: list[float]) -> dict | None:
    values = [v for v in values if v is not None]
    if not values:
        return None
    out = {"n": len(values), "median": round(statistics.median(values), 2),
           "mean": round(statistics.fmean(values), 2),
           "min": round(min(values), 2), "max": round(max(values), 2)}
    out["sd"] = round(statistics.stdev(values), 2) if len(values) > 1 else None
    return out


def _factor_column(factor: dict, columns: list[str], sheet: list[dict]) -> str | None:
    """The sample sheet column holding a design factor: the one named after it, else
    the one whose values are exactly the factor's levels (water_ph for factor ph)."""
    if factor["factor_id"] in columns:
        return factor["factor_id"]
    levels = {str(v) for v in factor.get("levels") or []}
    for col in columns:
        if levels and {r.get(col, "") for r in sheet} == levels:
            return col
    return None


def summarize(study: dict, sheet: list[dict], qc: list[dict], manifest: dict, *,
              part: str | None, run_name: str | None, today: str) -> dict:
    passed = [r for r in qc if r.get("qc_pass") == "true"]
    meth = {r["sample_id"]: _float(r.get("cpg_meth_pct")) for r in passed}
    params = manifest.get("params", {})
    out = {
        "study_id": study.get("study_id") or params.get("study_id"),
        "part": part,
        "generated_on": today,
        "nextflow_run": run_name or params.get("run_name"),
        "commit": params.get("commit"),
        "genome_assembly": params.get("genome_assembly"),
        "chemistry": params.get("chemistry"),
        "libraries": len(qc),
        "passed": len(passed),
        "failed": [{"sample_id": r["sample_id"], "reasons": r.get("qc_reasons", "")}
                   for r in qc if r.get("qc_pass") != "true"],
        "qc": {k: _stats([_float(r.get(k)) for r in qc])
               for k in ("mapping_pct", "duplication_pct", "conversion_pct", "cpgs_5x", "mean_depth")},
        "cpg_methylation_pct": _stats(list(meth.values())),
        "by_factor": [],
        "by_condition": [],
    }
    columns = list(sheet[0].keys()) if sheet else []
    by_id = {r["sample_id"]: r for r in sheet}
    for factor in study.get("factors") or []:
        col = _factor_column(factor, columns, sheet)
        if not col:
            continue
        levels = [str(v) for v in factor.get("levels") or []]
        rows = []
        for level in levels:
            s = _stats([v for sid, v in meth.items() if str(by_id.get(sid, {}).get(col, "")) == level])
            rows.append({"level": level, **(s or {"n": 0})})
        out["by_factor"].append({"factor_id": factor["factor_id"], "column": col,
                                 "reference_level": str(factor.get("reference_level")), "levels": rows})
    conditions = sorted({r.get("condition") for r in sheet if r.get("condition")})
    for c in conditions:
        s = _stats([v for sid, v in meth.items() if by_id.get(sid, {}).get("condition") == c])
        out["by_condition"].append({"condition": c, **(s or {"n": 0})})
    return out


def sample_files(results: Path, sample_id: str) -> list[Path]:
    """The files kept on the data host for one library: its CpG calls and its QC reports."""
    files = [results / "cpg" / f"{sample_id}{CPG_SUFFIX}"]
    qc = results / "qc"
    if qc.is_dir():
        files += sorted(p for p in qc.iterdir()
                        if p.is_file() and (p.name.startswith(f"{sample_id}_") or p.name.startswith(f"{sample_id}.")))
    return files


class DataHost:
    """Where methylomes go: ``user@host:/path`` (SSH + rsync, key auth) or a local folder."""

    def __init__(self, target: str):
        host, sep, path = target.partition(":")
        if sep and host and path.startswith("/"):
            self.host, self.path = host, path.rstrip("/")
        else:
            self.host, self.path = None, target.rstrip("/")
        self.label = target

    def _ssh(self, command: str) -> str:
        proc = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", self.host, command],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RecordError(f"ssh {self.host} failed: {proc.stderr.strip() or proc.returncode}")
        return proc.stdout

    def put(self, files: list[Path], subdir: str) -> None:
        dest = f"{self.path}/{subdir}"
        if self.host is None:
            Path(dest).mkdir(parents=True, exist_ok=True)
            for f in files:
                shutil.copy2(f, Path(dest) / f.name)
        else:
            self._ssh(f"mkdir -p {shlex.quote(dest)}")
            proc = subprocess.run(
                ["rsync", "-tp", "--partial", "--chmod=D755,F644", "-e", "ssh -o BatchMode=yes",
                 *map(str, files), f"{self.host}:{dest}/"], capture_output=True, text=True)
            if proc.returncode != 0:
                raise RecordError(f"rsync to {self.host}:{dest} failed: {proc.stderr.strip()}")
        sizes = self.sizes(subdir)
        wrong = [f.name for f in files if sizes.get(f.name) != f.stat().st_size]
        if wrong:
            raise RecordError(f"{len(wrong)} file(s) did not arrive whole in {dest}: {', '.join(wrong[:5])}")

    def sizes(self, subdir: str) -> dict[str, int]:
        dest = f"{self.path}/{subdir}"
        if self.host is None:
            d = Path(dest)
            return {p.name: p.stat().st_size for p in d.iterdir() if p.is_file()} if d.is_dir() else {}
        out = self._ssh(f"d={shlex.quote(dest)}; if [ -d \"$d\" ]; then find \"$d\" -maxdepth 1 -type f -printf '%f\\t%s\\n'; fi")
        return {n: int(s) for n, s in (ln.split("\t") for ln in out.splitlines() if ln)}


def _data_host(repo: Path) -> dict:
    path = repo / "data" / "data_host.yaml"
    return yaml.safe_load(path.read_text()) if path.exists() else {}


def _samplesheet(repo: Path, study_id: str, part: str | None) -> list[dict]:
    base = repo / "data" / "studies" / study_id
    for p in ([base / "parts" / part / "samplesheet.tsv"] if part else []) + [base / "samplesheet.tsv"]:
        if p.exists():
            return read_tsv(p)
    return []


def record_run(repo, results, *, run_name: str | None = None, upload: bool = True,
               upload_target: str | None = None, commit: bool = True, push: bool = True,
               remote: str = "origin", today: str | None = None, mark_ready: bool = True) -> dict:
    repo, results = Path(repo), Path(results)
    today = today or date.today().isoformat()
    branch = current_branch(repo)
    parsed = parse_branch(branch)
    if parsed is None:
        raise ClaimError(f"{branch!r} is not a study branch; run this inside the study's worktree.")
    _, study_id, part = parsed
    spath = status_path(repo, study_id, part)
    if not spath.exists():
        raise ClaimError(f"{spath.relative_to(repo)} is missing on this branch.")
    status = yaml.safe_load(spath.read_text())

    for name in ("run_manifest.json", "qc_summary.tsv"):
        if not (results / name).exists():
            raise RecordError(f"{results / name} is missing; did the run finish?")
    qc = read_tsv(results / "qc_summary.tsv")
    if not qc:
        raise RecordError(f"{results / 'qc_summary.tsv'} lists no libraries.")
    missing = [r["sample_id"] for r in qc if not (results / "cpg" / f"{r['sample_id']}{CPG_SUFFIX}").exists()]
    if missing:
        raise RecordError(f"no CpG calls in {results / 'cpg'} for {len(missing)} librar(ies): {', '.join(missing[:5])}")
    manifest = json.loads((results / "run_manifest.json").read_text())
    run_name = run_name or manifest.get("params", {}).get("run_name")

    dest = spath.parent
    shutil.copyfile(results / "qc_summary.tsv", dest / "qc_summary.tsv")
    shutil.copyfile(results / "run_manifest.json", dest / "run_manifest.json")

    host = _data_host(repo)
    url_base = (host.get("url_base") or "").rstrip("/")
    prefix = f"methylomes/{study_id}" + (f"/parts/{part}" if part else "")
    per_sample = {r["sample_id"]: sample_files(results, r["sample_id"]) for r in qc}
    rows = []
    for sid, files in per_sample.items():
        for f in files:
            rows.append({"sample_id": sid, "file": f.name, "bytes": f.stat().st_size, "sha256": sha256_file(f),
                         "url": f"{url_base}/{prefix}/{sid}/{f.name}" if url_base else ""})
    with open(dest / METHYLOME_MANIFEST, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    uploaded, upload_error = False, None
    target = upload_target or host.get("upload_target")
    if upload and target:
        store = DataHost(target)
        try:
            for sid, files in per_sample.items():
                store.put(files, f"{prefix}/{sid}")
            uploaded = True
        except (RecordError, OSError) as exc:
            upload_error = str(exc)
    elif upload:
        upload_error = "data/data_host.yaml has no upload_target"

    study_yaml = studies_dir(repo) / f"{study_id}.yaml"
    study = yaml.safe_load(study_yaml.read_text()) if study_yaml.exists() else {}
    summary = summarize(study, _samplesheet(repo, study_id, part), qc, manifest,
                        part=part, run_name=run_name, today=today)
    summary["methylomes"] = {"files": len(rows), "bytes": sum(r["bytes"] for r in rows),
                             "on_data_host": uploaded, "url": f"{url_base}/{prefix}/" if url_base else None}
    (dest / SUMMARY_NAME).write_text(
        "# Written by `whale record-run` from this study's qc_summary.tsv; read by the dashboard.\n"
        + yaml.safe_dump(summary, sort_keys=False))

    n, n_pass = len(qc), summary["passed"]
    mapping = (summary["qc"]["mapping_pct"] or {}).get("median")
    notes = {
        "aligned": f"{n} libraries aligned; median mapping {mapping}%",
        "called": f"CpG calls for {n} libraries (run {run_name})",
        "qc_passed": f"{n_pass} of {n} libraries pass the QC gate"
                     + (f"; failed: {', '.join(f['sample_id'] for f in summary['failed'])}" if summary["failed"] else ""),
        "summarized": f"summaries committed; {len(rows)} files on {host.get('host') or target}",
    }
    order = stage_order(repo)
    machine = machine_name()
    added = []
    for stage in RECORDED_STAGES:
        if order[stage] <= order.get(status["stage"], -1):
            continue
        if stage == "qc_passed" and n_pass == 0:
            break
        if stage == "summarized" and not uploaded:
            break
        status = advance_status(status, stage, repo=repo, machine=machine, today=today,
                                note=notes[stage], nextflow_run=run_name)
        added.append(stage)
    spath.write_text(dump_status(status))

    result = {"study_id": study_id, "part": part, "branch": branch, "stage": status["stage"], "added": added,
              "libraries": n, "passed": n_pass, "uploaded": uploaded, "upload_error": upload_error,
              "committed": False, "pushed": False, "pr_ready": False}
    if commit:
        files = [spath, dest / "qc_summary.tsv", dest / "run_manifest.json", dest / METHYLOME_MANIFEST,
                 dest / SUMMARY_NAME]
        _run(["git", "add", *[str(f.relative_to(repo)) for f in files]], repo)
        if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=repo).returncode != 0:
            label = f"{study_id}{f' part {part}' if part else ''}"
            msg = f"{label}: {status['stage']} (record-run {run_name})" if added else f"{label}: summaries from {run_name}"
            _run(["git", "commit", "--quiet", "-m", msg], repo)
            result["committed"] = True
        if push:
            _run(["git", "push", "--quiet", remote, branch], repo)
            result["pushed"] = True
    if mark_ready and push and status["stage"] == "summarized" and shutil.which("gh"):
        proc = subprocess.run(["gh", "pr", "ready", branch], cwd=repo, capture_output=True, text=True)
        result["pr_ready"] = proc.returncode == 0
    return result
