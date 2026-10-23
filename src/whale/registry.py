"""The study registry: every study YAML, the generated index, and the backlog."""
from __future__ import annotations

import csv
from pathlib import Path

import yaml

from whale.paths import studies_dir

INDEX_COLUMNS = [
    "study_id", "bioproject", "chemistry", "n_runs", "n_libraries", "sample_unit",
    "tissues", "life_stages", "paired_rnaseq", "rnaseq_bioproject", "goals",
    "fastq_size_gb", "quality_flags", "registered_by", "date_registered",
]


def study_files(repo: Path) -> list[Path]:
    return sorted(p for p in studies_dir(repo).glob("*.yaml") if not p.name.startswith("_"))


def load_studies(repo: Path) -> dict[str, dict]:
    out = {}
    for path in study_files(repo):
        study = yaml.safe_load(path.read_text()) or {}
        out[study.get("study_id") or path.stem] = study
    return out


def index_row(study: dict) -> dict:
    s = study.get("samples") or {}
    r = study.get("paired_rnaseq") or {}
    p = study.get("provenance") or {}
    return {
        "study_id": study.get("study_id"),
        "bioproject": (study.get("accessions") or {}).get("bioproject"),
        "chemistry": study.get("chemistry"),
        "n_runs": s.get("n_runs"),
        "n_libraries": s.get("n_libraries"),
        "sample_unit": s.get("sample_unit"),
        "tissues": ";".join(s.get("tissues") or []),
        "life_stages": ";".join(s.get("life_stages") or []),
        "paired_rnaseq": r.get("relationship"),
        "rnaseq_bioproject": r.get("bioproject"),
        "goals": ";".join(str(g) for g in study.get("goals") or []),
        "fastq_size_gb": (study.get("data_availability") or {}).get("fastq_size_gb"),
        "quality_flags": ";".join(study.get("quality_flags") or []),
        "registered_by": p.get("registered_by"),
        "date_registered": p.get("date_registered"),
    }


def build_index(repo: Path, out: Path | None = None) -> Path:
    out = out or (repo / "registry" / "study_registry.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=INDEX_COLUMNS)
        w.writeheader()
        for study in load_studies(repo).values():
            w.writerow({k: ("" if v is None else v) for k, v in index_row(study).items()})
    return out


def load_backlog(repo: Path) -> list[dict]:
    path = repo / "registry" / "backlog.tsv"
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
