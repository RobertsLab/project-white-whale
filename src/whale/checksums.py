"""File checksums: MD5 to match what ENA, NCBI and the mirror ledger publish
(integrity, not security), SHA-256 for what this project records itself."""
from __future__ import annotations

import hashlib
from pathlib import Path

CHUNK = 1 << 20


def file_digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_file(path: Path) -> str:
    return file_digest(path, "md5")


def sha256_file(path: Path) -> str:
    return file_digest(path, "sha256")
