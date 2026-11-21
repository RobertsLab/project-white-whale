"""A study's paired RNA-seq: its sample sheet, and recording a finished run.

Goal 2 (PLAN.md §7.2) relates methylation to expression for the studies whose
registration names RNA-seq under ``paired_rnaseq.bioproject``. Their RNA-seq
lives beside the methylation data, in data/studies/<STUDY_ID>/rnaseq/:

* ``samplesheet.tsv``, ``fastq_manifest.tsv``, ``ena_provenance.json``
  written by `whale fetch-rnaseq-samplesheet` from ENA (whale.ena, keeping
  RNA-Seq runs), in the column contract workflows/rnaseq reads;
* after a run of workflows/rnaseq, `whale record-rnaseq` (called by the
  workflow when it succeeds) adds what is small enough to commit: the gene x
  sample matrices (``gene_tpm.tsv.gz``, ``gene_counts.tsv.gz``,
  ``gene_length.tsv.gz``), ``qc_summary.tsv``, ``run_manifest.json`` and
  ``rnaseq_summary.yaml``; and ``rnaseq_manifest.tsv``, the size, SHA-256 and
  data-host URL of everything else (each library's quant.sf and Salmon/fastp/
  FastQC reports, the transcript-level matrices, the MultiQC report), which it
  uploads to the data host under rnaseq/<STUDY_ID>/<sample_id>/ (study-level
  files under rnaseq/<STUDY_ID>/_study/).

Recording follows `whale record-run` (whale.results): it is safe to repeat,
writes gzip files byte-identically for the same input, and a failed upload is
finished by running it again. It does not move status.yaml: the pipeline
stages there are the methylation workflow's. It commits on the study's own
branch (study/<ACCESSION>/<STUDY_ID>) or on ``rnaseq/<STUDY_ID>`` for a study
whose methylation has already merged, never on another branch.
"""
from __future__ import annotations

import csv
import gzip
import json
import shutil
import subprocess
from datetime import date
from pathlib import Path

import yaml

from whale.checksums import sha256_file
from whale.claims import ClaimError, _run, current_branch, parse_branch
from whale.ena import build_samplesheet
from whale.paths import studies_dir, study_data_dir
from whale.results import DataHost, RecordError, _data_host, _float, _stats, read_tsv

RNASEQ_DIR = "rnaseq"
RNASEQ_STRATEGIES = ["RNA-Seq"]
GENE_MATRICES = ["gene_tpm.tsv.gz", "gene_counts.tsv.gz", "gene_length.tsv.gz"]
STUDY_FILES = ["matrices/transcript_tpm.tsv.gz", "matrices/transcript_counts.tsv.gz",
               "matrices/matrices_report.json", "multiqc_report.html"]
SALMON_FILES = ["quant.sf", "aux_info/meta_info.json", "lib_format_counts.json", "cmd_info.json",
                "logs/salmon_quant.log"]
MANIFEST_NAME = "rnaseq_manifest.tsv"
MANIFEST_COLUMNS = ["sample_id", "file", "bytes", "sha256", "url"]
SUMMARY_NAME = "rnaseq_summary.yaml"
STUDY_SAMPLE = "_study"
DETECTED_TPM = 1.0


def rnaseq_dir(repo: Path, study_id: str) -> Path:
    return study_data_dir(repo, study_id) / RNASEQ_DIR


def load_study(repo: Path, study_id: str) -> dict:
    path = studies_dir(repo) / f"{study_id}.yaml"
    if not path.exists():
        raise RecordError(f"{path.relative_to(repo)} does not exist; register the study first.")
    return yaml.safe_load(path.read_text()) or {}


def paired_bioproject(study: dict) -> tuple[str, str]:
    """(RNA-seq BioProject, relationship) from a study's paired_rnaseq block."""
    block = study.get("paired_rnaseq") or {}
    bioproject, relationship = block.get("bioproject"), block.get("relationship") or "unknown"
    if not bioproject or relationship == "none":
        raise RecordError(f"{study.get('study_id')}: paired_rnaseq names no BioProject "
                          f"(relationship {relationship!r}); record it in the study YAML first.")
    return bioproject, relationship


def fetch_rnaseq_samplesheet(repo: Path, study_id: str, *, condition_attributes: list[str],
                             bioproject: str | None = None, library_strategies: list[str] | None = None,
                             **kwargs) -> dict:
    """Write data/studies/<STUDY_ID>/rnaseq/{samplesheet.tsv, fastq_manifest.tsv, ena_provenance.json}.

    The BioProject defaults to the study's ``paired_rnaseq.bioproject``. Other
    keyword arguments go to whale.ena.build_samplesheet (extra_attributes,
    label_field, label_pattern, design_table, design_join, design_source).
    """
    study = load_study(repo, study_id)
    relationship = (study.get("paired_rnaseq") or {}).get("relationship") or "unknown"
    if bioproject is None:
        bioproject, relationship = paired_bioproject(study)
    report = build_samplesheet(bioproject, study_id=study_id, out_dir=rnaseq_dir(repo, study_id),
                               condition_attributes=condition_attributes,
                               library_strategies=library_strategies or RNASEQ_STRATEGIES, repo=repo, **kwargs)
    return {**report, "relationship": relationship}


# --- recording a run ------------------------------------------------------------

def sample_files(results: Path, sample_id: str) -> list[Path]:
    """One library's files for the data host: Salmon's quantification and logs, fastp and FastQC reports."""
    salmon = results / "salmon" / sample_id
    files = [salmon / f for f in SALMON_FILES if (salmon / f).exists()]
    qc = results / "qc"
    if qc.is_dir():
        files += sorted(p for p in qc.iterdir()
                        if p.is_file() and (p.name.startswith(f"{sample_id}_") or p.name.startswith(f"{sample_id}.")))
    return files


def _matrix(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    with gzip.open(path, "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0][1:], [r[0] for r in rows[1:]], [[float(v) for v in r[1:]] for r in rows[1:]]


def summarize(qc: list[dict], manifest: dict, results: Path, *, run_name: str | None, today: str) -> dict:
    params = manifest.get("params", {})
    passed = [r for r in qc if r.get("qc_pass") == "true"]
    samples, genes, tpm = _matrix(results / "matrices" / "gene_tpm.tsv.gz")
    keep = [i for i, s in enumerate(samples) if s in {r["sample_id"] for r in passed}]
    detected_all = sum(all(row[i] >= DETECTED_TPM for i in keep) for row in tpm) if keep else 0
    strands: dict[str, int] = {}
    for r in qc:
        key = r.get("strandedness") or "unknown"
        strands[key] = strands.get(key, 0) + 1
    index = manifest.get("salmon_index", {})
    return {
        "study_id": params.get("study_id"),
        "generated_on": today,
        "nextflow_run": run_name or params.get("run_name"),
        "commit": params.get("commit"),
        "genome_assembly": params.get("genome_assembly"),
        "quantifier": {"salmon": index.get("salmon"), "decoys": index.get("decoys"), "kmer": index.get("kmer")},
        "libraries": len(qc),
        "passed": len(passed),
        "failed": [{"sample_id": r["sample_id"], "reasons": r.get("qc_reasons", "")}
                   for r in qc if r.get("qc_pass") != "true"],
        "layouts": sorted({r.get("layout") for r in qc if r.get("layout")}),
        "strandedness": dict(sorted(strands.items())),
        "mixed_strandedness": len(set(strands) - {"unknown"}) > 1,
        "qc": {k: _stats([_float(r.get(k)) for r in qc])
               for k in ("reads_in", "mapping_pct", "decoy_pct", "fragments_mapped", "genes_tpm_ge1")},
        "genes": {"universe": len(genes), "tpm_ge1_in_every_passing_library": detected_all,
                  "matrices": GENE_MATRICES},
    }


def check_branch(repo: Path, study_id: str) -> str:
    branch = current_branch(repo)
    parsed = parse_branch(branch)
    if (parsed and parsed[1] == study_id) or branch == f"{RNASEQ_DIR}/{study_id}":
        return branch
    raise ClaimError(f"{branch!r} is neither {study_id}'s study branch nor {RNASEQ_DIR}/{study_id}; "
                     f"check one out (git switch -c {RNASEQ_DIR}/{study_id} origin/main), or pass --no-commit.")


def record_rnaseq(repo, results, *, run_name: str | None = None, upload: bool = True,
                  upload_target: str | None = None, commit: bool = True, push: bool = True,
                  remote: str = "origin", today: str | None = None) -> dict:
    repo, results = Path(repo), Path(results)
    today = today or date.today().isoformat()
    for name in ("run_manifest.json", "qc_summary.tsv", *(f"matrices/{m}" for m in GENE_MATRICES)):
        if not (results / name).exists():
            raise RecordError(f"{results / name} is missing; did the run finish?")
    manifest = json.loads((results / "run_manifest.json").read_text())
    params = manifest.get("params", {})
    if params.get("workflow") != "rnaseq":
        raise RecordError(f"{results / 'run_manifest.json'} is not from workflows/rnaseq.")
    study_id = params.get("study_id")
    if not study_id:
        raise RecordError(f"{results / 'run_manifest.json'} names no study_id.")
    if commit:
        branch = check_branch(repo, study_id)
    else:
        try:
            branch = current_branch(repo)
        except ClaimError:      # not a git checkout, or no commit yet: nothing to commit to anyway
            branch = None
    qc = read_tsv(results / "qc_summary.tsv")
    if not qc:
        raise RecordError(f"{results / 'qc_summary.tsv'} lists no libraries.")
    missing = [r["sample_id"] for r in qc if not (results / "salmon" / r["sample_id"] / "quant.sf").exists()]
    if missing:
        raise RecordError(f"no quant.sf under {results / 'salmon'} for {len(missing)} librar(ies): {', '.join(missing[:5])}")
    run_name = run_name or params.get("run_name")

    dest = rnaseq_dir(repo, study_id)
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("qc_summary.tsv", "run_manifest.json"):
        shutil.copyfile(results / name, dest / name)
    for name in GENE_MATRICES:
        shutil.copyfile(results / "matrices" / name, dest / name)

    host = _data_host(repo)
    url_base = (host.get("url_base") or "").rstrip("/")
    prefix = f"{RNASEQ_DIR}/{study_id}"
    groups = {r["sample_id"]: sample_files(results, r["sample_id"]) for r in qc}
    groups[STUDY_SAMPLE] = [results / f for f in STUDY_FILES if (results / f).exists()]
    rows = []
    for sid, files in groups.items():
        names = [f.name for f in files]
        dupes = sorted({n for n in names if names.count(n) > 1})
        if dupes:
            raise RecordError(f"{sid}: several files named {dupes} would overwrite each other on the data host")
        for f in files:
            rows.append({"sample_id": sid, "file": f.name, "bytes": f.stat().st_size, "sha256": sha256_file(f),
                         "url": f"{url_base}/{prefix}/{sid}/{f.name}" if url_base else ""})
    with open(dest / MANIFEST_NAME, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    uploaded, upload_error = False, None
    target = upload_target or host.get("upload_target")
    if upload and target:
        store = DataHost(target)
        try:
            for sid, files in groups.items():
                if files:
                    store.put(files, f"{prefix}/{sid}")
            uploaded = True
        except (RecordError, OSError) as exc:
            upload_error = str(exc)
    elif upload:
        upload_error = "data/data_host.yaml has no upload_target"

    summary = summarize(qc, manifest, results, run_name=run_name, today=today)
    summary["files"] = {"files": len(rows), "bytes": sum(r["bytes"] for r in rows), "on_data_host": uploaded,
                        "url": f"{url_base}/{prefix}/" if url_base else None}
    (dest / SUMMARY_NAME).write_text(
        "# Written by `whale record-rnaseq` from this study's RNA-seq run (workflows/rnaseq).\n"
        + yaml.safe_dump(summary, sort_keys=False))

    result = {"study_id": study_id, "branch": branch, "libraries": len(qc), "passed": summary["passed"],
              "files": len(rows), "uploaded": uploaded, "upload_error": upload_error,
              "mixed_strandedness": summary["mixed_strandedness"], "committed": False, "pushed": False}
    if commit:
        files = [dest / n for n in ("qc_summary.tsv", "run_manifest.json", MANIFEST_NAME, SUMMARY_NAME,
                                    *GENE_MATRICES)]
        _run(["git", "add", *[str(f.relative_to(repo)) for f in files]], repo)
        if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=repo).returncode != 0:
            _run(["git", "commit", "--quiet", "-m",
                  f"{study_id}: RNA-seq quantified (record-rnaseq {run_name}); "
                  f"{summary['passed']} of {len(qc)} libraries pass"], repo)
            result["committed"] = True
        if push:
            _run(["git", "push", "--quiet", "-u", remote, branch], repo)
            result["pushed"] = True
    return result
