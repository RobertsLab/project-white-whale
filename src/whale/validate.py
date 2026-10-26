"""Validate study registrations and status files.

Errors are hard failures: schema violations, unknown vocabulary terms,
inconsistent fields. Warnings are fields still null or combinations a curator
should look at twice. A freshly claimed study has many warnings and no errors.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import jsonschema
import yaml

from whale.paths import schemas_dir
from whale.vocab import load_vocab


@dataclass
class ValidationResult:
    path: str
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _schema_errors(doc, schema_name: str, repo: Path) -> list[str]:
    schema = json.loads((schemas_dir(repo) / schema_name).read_text())
    validator = jsonschema.Draft7Validator(schema)
    out = []
    for err in sorted(validator.iter_errors(doc), key=lambda e: list(e.path)):
        loc = "/".join(str(p) for p in err.path) or "<root>"
        out.append(f"{loc}: {err.message}")
    return out


def _check_term(errors: list[str], repo: Path, vocab: str, value, where: str) -> None:
    if value is not None and value not in load_vocab(repo, vocab):
        errors.append(f"{where}: {value!r} is not a term in {vocab}.yaml")


def validate_study(study: dict, repo: Path, *, filename_stem: str | None = None) -> tuple[list[str], list[str]]:
    errors = _schema_errors(study, "study.schema.json", repo)
    if errors:
        return errors, []
    warnings: list[str] = []

    if filename_stem and filename_stem != "_TEMPLATE" and study["study_id"] != filename_stem:
        errors.append(f"study_id {study['study_id']!r} does not match the file name {filename_stem!r}")

    _check_term(errors, repo, "assay_chemistry", study["chemistry"], "chemistry")
    samples = study["samples"]
    _check_term(errors, repo, "sample_unit", samples["sample_unit"], "samples/sample_unit")
    for t in samples["tissues"]:
        _check_term(errors, repo, "tissue", t, "samples/tissues")
    for s in samples["life_stages"]:
        _check_term(errors, repo, "life_stage", s, "samples/life_stages")
    _check_term(errors, repo, "rnaseq_relationship", study["paired_rnaseq"]["relationship"],
                "paired_rnaseq/relationship")
    for f in study["factors"]:
        _check_term(errors, repo, "stressor", f["factor_standardized"],
                    f"factors/{f['factor_id']}/factor_standardized")
        ref = f.get("reference_level")
        if ref is not None and ref not in f["levels"]:
            errors.append(f"factors/{f['factor_id']}: reference_level {ref!r} is not one of its levels")
    for q in study["quality_flags"]:
        _check_term(errors, repo, "quality_flags", q, "quality_flags")

    ids = [f["factor_id"] for f in study["factors"]]
    if len(ids) != len(set(ids)):
        errors.append("factors: factor_id values must be unique")

    if samples["sample_unit"] == "pool" and samples.get("pool_size") is None:
        warnings.append("samples/sample_unit is pool but pool_size is not recorded")
    if samples["sample_unit"] == "pool" and "pooled_libraries" not in study["quality_flags"]:
        errors.append("samples/sample_unit is pool, so quality_flags must include pooled_libraries")

    chem = study["chemistry"]
    if chem in load_vocab(repo, "assay_chemistry"):
        tier = load_vocab(repo, "assay_chemistry")[chem]["tier"]
        if tier == "enrichment" and "enrichment_assay" not in study["quality_flags"]:
            errors.append(f"chemistry {chem} is an enrichment assay, so quality_flags must include enrichment_assay")
        if chem == "rrbs" and "rrbs_subset" not in study["quality_flags"]:
            errors.append("chemistry rrbs needs quality flag rrbs_subset")

    rel = study["paired_rnaseq"]["relationship"]
    goals = study.get("goals") or []
    if 2 in goals and rel in ("none", "unknown"):
        errors.append(f"goals includes 2 (methylation and expression) but paired_rnaseq/relationship is {rel}")
    if rel in ("same_animals", "same_design", "related") and not study["paired_rnaseq"].get("bioproject"):
        warnings.append("paired_rnaseq has a relationship but no bioproject")

    n_runs, n_libs = samples["n_runs"], samples["n_libraries"]
    if n_runs is not None and n_libs is not None and n_libs > n_runs:
        errors.append(f"samples: n_libraries ({n_libs}) cannot exceed n_runs ({n_runs})")

    for path, value in [
        ("accessions/bioproject", study["accessions"]["bioproject"]),
        ("citation", study["citation"]),
        ("chemistry", chem),
        ("library_layout", study["library_layout"]),
        ("samples/n_runs", n_runs),
        ("samples/sample_unit", samples["sample_unit"]),
        ("data_availability/raw_data_available", study["data_availability"]["raw_data_available"]),
        ("data_availability/fastq_size_gb", study["data_availability"]["fastq_size_gb"]),
        ("limitations", study["limitations"]),
    ]:
        if value is None:
            warnings.append(f"{path} is not filled in")
    if not samples["tissues"]:
        warnings.append("samples/tissues is empty")
    if not study["factors"]:
        warnings.append("factors is empty: record the design, even if it is a single group")
    if not goals:
        warnings.append("goals is empty")
    return errors, warnings


def validate_study_file(path: Path, repo: Path) -> ValidationResult:
    path = Path(path)
    res = ValidationResult(str(path))
    try:
        study = yaml.safe_load(path.read_text())
    except yaml.YAMLError as exc:
        res.errors.append(f"could not parse YAML: {exc}")
        return res
    if not isinstance(study, dict):
        res.errors.append("file is not a YAML mapping")
        return res
    res.errors, res.warnings = validate_study(study, repo, filename_stem=path.stem)
    if not res.errors:
        res.errors += config_chemistry_errors(study, repo)
    return res


def config_chemistry_errors(study: dict, repo: Path) -> list[str]:
    """The study's run configs (whole and parts) must set the YAML's chemistry.

    `whale start-study` writes them with `chemistry = null`, and the workflow
    refuses that only once launched, after FASTQ staging.
    """
    from whale.claims import config_chemistry

    want = study["chemistry"]
    if want is None:
        return []
    errors = []
    sid = study["study_id"]
    for cfg in sorted([Path(repo) / "config" / f"{sid}.config",
                       *(Path(repo) / "config").glob(f"{sid}.part-*.config")]):
        if not cfg.exists():
            continue
        got = config_chemistry(cfg.read_text())
        if got != want:
            errors.append(f"config/{cfg.name} sets chemistry = {got!r} but this YAML says {want!r}; "
                          f"set chemistry = '{want}' there")
    return errors


def validate_status(status: dict, repo: Path) -> list[str]:
    errors = _schema_errors(status, "status.schema.json", repo)
    if errors:
        return errors
    _check_term(errors, repo, "pipeline_stage", status["stage"], "stage")
    _check_term(errors, repo, "compute_profile", status["profile"], "profile")
    for i, h in enumerate(status["history"]):
        _check_term(errors, repo, "pipeline_stage", h["stage"], f"history/{i}/stage")
    if status["history"] and status["history"][-1]["stage"] != status["stage"]:
        errors.append("stage does not match the last history entry")
    return errors
