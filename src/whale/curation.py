"""Curated per-library sample metadata: data/studies/<STUDY_ID>/sample_curation.tsv.

The sample sheet is the deposit as ENA and BioSample have it; its columns differ
by study and rarely hold the attributes the atlas compares (tissue, life stage,
sex, ploidy, pooling) in the shared vocabularies. A curated table holds them, one
row per library of the sample sheet, read from the deposit and the paper, with
the source of every value. When a study has one, the atlas takes these
attributes from it (src/whale/atlas.py, `library_metadata`); see
docs/sample_curation.md.

`whale validate-curation` checks every curated table: the same libraries as the
sample sheet, vocabulary terms, pool sizes, and a source for each row.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

from whale.vocab import load_vocab

FILENAME = "sample_curation.tsv"
COLUMNS = ["sample_id", "animal_id", "tissue", "life_stage", "age", "sex", "ploidy", "sample_unit",
           "pool_size", "population", "treatment", "confidence", "sources", "notes"]
UNKNOWN = "unknown"
SEXES = {"female", "male", "hermaphrodite", UNKNOWN}
PLOIDIES = {"diploid", "triploid", "tetraploid", UNKNOWN}
CONFIDENCE = {"high", "medium", "low"}
# column -> vocabulary file; `unknown` is accepted in each as well
VOCAB_COLUMNS = {"tissue": "tissue", "life_stage": "life_stage", "sample_unit": "sample_unit"}


def curation_path(study_dir: Path) -> Path:
    return Path(study_dir) / FILENAME


def read_rows(path: Path) -> tuple[list[str], list[dict]]:
    """(header, rows) of a curated table; ([], []) when there is none."""
    path = Path(path)
    if not path.exists():
        return [], []
    reader = csv.DictReader(io.StringIO(path.read_text()), delimiter="\t")
    return list(reader.fieldnames or []), list(reader)


def read_curation(study_dir: Path) -> dict[str, dict]:
    """sample_id -> curated row; empty when the study has no curated table."""
    _, rows = read_rows(curation_path(study_dir))
    return {r["sample_id"]: r for r in rows}


def validate_curation(path: Path, sheet_ids: list[str], repo: Path) -> list[str]:
    """Errors in one curated table, checked against its study's sample sheet."""
    header, rows = read_rows(path)
    if not header:
        return [f"{path}: missing or empty"]
    errors = []
    if header != COLUMNS:
        errors.append(f"columns must be exactly {', '.join(COLUMNS)}; got {', '.join(header)}")
        return errors
    ids = [r["sample_id"] for r in rows]
    if ids != list(sheet_ids):
        extra = [i for i in ids if i not in set(sheet_ids)]
        missing = [i for i in sheet_ids if i not in set(ids)]
        if extra or missing:
            errors.append("libraries differ from samplesheet.tsv"
                          + (f"; not in the sample sheet: {', '.join(extra[:5])}" if extra else "")
                          + (f"; not curated: {', '.join(missing[:5])}" if missing else ""))
        elif len(set(ids)) != len(ids):
            errors.append("a sample_id appears more than once")
        else:
            errors.append("rows must be in the sample sheet's order")
    vocabs = {col: set(load_vocab(repo, name)) | {UNKNOWN} for col, name in VOCAB_COLUMNS.items()}
    for n, r in enumerate(rows, start=2):
        where = f"line {n} ({r['sample_id']})"
        for col, terms in vocabs.items():
            if r[col] not in terms:
                errors.append(f"{where}: {col} {r[col]!r} is not a {VOCAB_COLUMNS[col]}.yaml term or {UNKNOWN!r}")
        if r["sex"] not in SEXES:
            errors.append(f"{where}: sex {r['sex']!r} must be one of {', '.join(sorted(SEXES))}")
        if r["ploidy"] not in PLOIDIES:
            errors.append(f"{where}: ploidy {r['ploidy']!r} must be one of {', '.join(sorted(PLOIDIES))}")
        if r["confidence"] not in CONFIDENCE:
            errors.append(f"{where}: confidence {r['confidence']!r} must be high, medium or low")
        pool = r["pool_size"].strip()
        if pool and not (pool.isdigit() and int(pool) >= 2):
            errors.append(f"{where}: pool_size {pool!r} must be a whole number of at least 2")
        if r["sample_unit"] == "pool" and not pool:
            errors.append(f"{where}: sample_unit is pool but pool_size is empty")
        if pool and r["sample_unit"] != "pool":
            errors.append(f"{where}: pool_size is set but sample_unit is {r['sample_unit']!r}, not pool")
        if not r["sources"].strip():
            errors.append(f"{where}: no sources")
    return errors
