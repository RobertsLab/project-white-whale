"""Controlled vocabularies in registry/controlled_vocabularies/<name>.yaml."""
from __future__ import annotations

from functools import cache
from pathlib import Path

import yaml

from whale.paths import vocab_dir


@cache
def _load(path: Path) -> dict:
    doc = yaml.safe_load(path.read_text()) or {}
    terms = doc.get("terms")
    if not isinstance(terms, dict):
        raise ValueError(f"{path} has no `terms` mapping")
    return terms


def load_vocab(repo: Path, name: str) -> dict:
    """term id -> term attributes."""
    return _load(vocab_dir(repo) / f"{name}.yaml")


def stage_order(repo: Path) -> dict[str, int]:
    return {k: v["order"] for k, v in load_vocab(repo, "pipeline_stage").items()}
