"""Locate the repository a command is acting on.

Study claims put each study in its own git worktree, and the package is
usually installed editable from the main checkout. So the repository is found
from the current directory (walking up to the first directory that has
`registry/controlled_vocabularies/`), not from where the package lives.
`WHALE_REPO` overrides both.
"""
from __future__ import annotations

import os
from pathlib import Path

MARKER = Path("registry/controlled_vocabularies")


def find_repo(start: Path | None = None) -> Path:
    env = os.environ.get("WHALE_REPO")
    if env:
        return Path(env).resolve()
    here = (start or Path.cwd()).resolve()
    for d in (here, *here.parents):
        if (d / MARKER).is_dir():
            return d
    package_repo = Path(__file__).resolve().parents[2]
    if (package_repo / MARKER).is_dir():
        return package_repo
    raise FileNotFoundError(
        f"No White Whale repository at or above {here}; run from a checkout or set WHALE_REPO."
    )


def vocab_dir(repo: Path) -> Path:
    return repo / MARKER


def studies_dir(repo: Path) -> Path:
    return repo / "registry" / "studies"


def schemas_dir(repo: Path) -> Path:
    return repo / "schemas"


def study_data_dir(repo: Path, study_id: str) -> Path:
    return repo / "data" / "studies" / study_id


def status_path(repo: Path, study_id: str, part: str | None = None) -> Path:
    base = study_data_dir(repo, study_id)
    return (base / "parts" / part / "status.yaml") if part else (base / "status.yaml")
