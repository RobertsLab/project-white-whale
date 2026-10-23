"""Keep a study PR inside the files its study (or study part) owns.

Study branches are worked on different machines at the same time. They merge
in any order without conflicts only if none edits a shared file, so a study
PR may touch only:

    whole study                               part <PART>
    registry/studies/<STUDY_ID>.yaml          data/studies/<STUDY_ID>/parts/<PART>/**
    config/<STUDY_ID>.config                  config/<STUDY_ID>.part-<PART>.config
    data/studies/<STUDY_ID>/** (not parts/)
    tests/test_<study_id>.py

A code, vocabulary or docs change the study needs goes in its own PR, which
merges first; the study branch then rebases onto main.

Adapted from AREE's validation/study_scope.py (github.com/sr320/AREE).
"""
from __future__ import annotations

from whale.claims import accessions_in, config_name, parse_branch


def owned(study_id: str, part: str | None) -> tuple[set[str], str, str | None]:
    """(exact files, directory prefix, excluded sub-prefix) a study branch may change."""
    if part:
        return {config_name(study_id, part)}, f"data/studies/{study_id}/parts/{part}/", None
    files = {
        f"registry/studies/{study_id}.yaml",
        config_name(study_id),
        f"tests/test_{study_id.lower()}.py",
    }
    return files, f"data/studies/{study_id}/", f"data/studies/{study_id}/parts/"


def check_scope(branch: str, changed: list[str], study: dict | None,
                registered: dict[str, set[str]]) -> tuple[list[str], list[str]]:
    """(errors, warnings) for a study PR.

    `changed` are paths changed relative to the base, `study` is the parsed
    study YAML (from the branch for a whole-study claim, from the base for a
    part; None if missing), and `registered` maps study_id -> accessions on the base.
    """
    parsed = parse_branch(branch)
    if parsed is None:
        return [f"{branch!r} is not a study branch (expected study/<ACCESSION>/<STUDY_ID>[.part-<PART>])."], []
    accession, study_id, part = parsed
    errors, warnings = [], []

    files, prefix, excluded = owned(study_id, part)
    outside = sorted(
        p for p in changed
        if p not in files and not (p.startswith(prefix) and not (excluded and p.startswith(excluded)))
    )
    if outside:
        errors.append(f"changes files outside {study_id}{f' part {part}' if part else ''}'s own paths; "
                      "move these to a separate PR:\n    " + "\n    ".join(outside))

    if study is None:
        errors.append(f"registry/studies/{study_id}.yaml is missing"
                      f"{' on the base branch (register the study before claiming parts)' if part else ''}.")
        return errors, warnings
    if study.get("study_id") != study_id:
        errors.append(f"branch names study {study_id} but the YAML has study_id {study.get('study_id')!r}.")
    if accession not in accessions_in(study):
        errors.append(f"branch accession {accession} does not appear in the YAML's `accessions` block.")

    for sid, accs in sorted(registered.items()):
        if sid != study_id and accession in accs:
            warnings.append(f"{accession} is also an accession of {sid} on main. That is fine only "
                            "if one deposit genuinely backs two studies; say so in the YAML.")
    return errors, warnings
