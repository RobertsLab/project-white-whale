"""Read-only access to White Whale's committed results, cached for a long-lived server.

Everything here reads files; nothing writes, and nothing calls git or the network.
Sources, all committed on main (PLAN.md §3, §6.2, §6.3, §7.1):

* registry/studies/*.yaml                         design, chemistry, tissues, quality flags
* data/studies/<ID>/study_summary.yaml            the recorded run (`whale record-run`): QC, mask
* data/studies/<ID>/status.yaml                   pipeline stage
* data/studies/<ID>/samplesheet.tsv, qc_summary.tsv
* data/studies/<ID>/genes/*.tsv.gz                pooled into the atlas by whale.atlas.pool
* data/studies/<ID>/sample_curation.tsv           curated sample attributes, when present
* data/studies/<ID>/differential/                 shared DSS model (`whale differential`), when present
* data/studies/<ID>/replication.json              replication verdicts (`whale replicate`), when present
* data/reference/gene_info/<ASSEMBLY>.tsv.gz      gene symbols and product descriptions
* site/atlas/data/atlas.json                      `whale build-atlas` output, when built locally

The pooled gene matrices are rebuilt in memory with the same code as `whale
build-atlas` (whale.atlas.pool and units), so the numbers match the atlas pages:
libraries of one BioSample are summed into one unit, and a unit's gene value
needs at least GENE_MIN_CPGS covered CpGs.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import threading
from functools import cached_property
from pathlib import Path

import numpy as np
import yaml

from whale import atlas
from whale.dashboard import GENE_MIN_CPGS, summaries_on_main
from whale.reference import gene_info_path, read_gene_info
from whale.registry import load_studies

REGIONS = ("gene_body", "promoter")
CURATED_FIELDS = ("tissue", "life_stage", "sex", "ploidy", "sample_unit", "pool_size")


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))


def _yaml(path: Path) -> dict:
    return (yaml.safe_load(path.read_text()) or {}) if path.exists() else {}


def parse_factors(text: str) -> dict[str, str]:
    """'ploidy=diploid;ph=8.2' -> {'ploidy': 'diploid', 'ph': '8.2'}"""
    out = {}
    for part in (text or "").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


class Atlas:
    """The pooled gene tables, per unit, with % methylation per region."""

    def __init__(self, repo: Path, curation: dict[str, dict[str, dict]]):
        pooled = atlas.pool(repo)
        u = atlas.units(pooled)
        self.genome_assembly = pooled["genome_assembly"]
        self.genes: list[str] = pooled["genes"]
        self.gene_index = {g: i for i, g in enumerate(self.genes)}
        self.studies: list[dict] = pooled["studies"]
        self.excluded: list[dict] = pooled["excluded"]
        self.libraries: list[dict] = pooled["libraries"]
        self.units: list[dict] = u["units"]
        self.pct = {r: atlas._pct(u["matrices"][(r, "meth")], u["matrices"][(r, "calls")],
                                  u["matrices"][(r, "cpgs")]) for r in REGIONS}
        lib_meta = {m["library_id"]: m for m in self.libraries}
        for unit in self.units:
            unit["factor_levels"] = parse_factors(unit.get("factors", ""))
            first = lib_meta.get(unit["libraries"][0], {})
            unit["chemistry"] = first.get("chemistry", "")
            unit["metadata_source"] = first.get("metadata_source") or "deposit"
            # Curated attributes win (whale build-atlas does the same once sample_curation lands).
            cur = curation.get(unit["study_id"], {}).get(unit["libraries"][0])
            if cur:
                for f in CURATED_FIELDS:
                    v = (cur.get(f) or "").strip()
                    unit[f] = "" if v == "unknown" else v
                unit["animal_id"] = cur.get("animal_id", "")
                unit["metadata_source"] = "curated"
        self.study_ids = [s["study_id"] for s in self.studies]
        self.columns = {sid: [j for j, m in enumerate(self.units) if m["study_id"] == sid]
                        for sid in self.study_ids}


class Store:
    """One repository's committed results. Reloads when the files it read change."""

    def __init__(self, repo: Path):
        self.repo = Path(repo).resolve()
        self._lock = threading.Lock()
        self._atlas: Atlas | None = None
        self._atlas_sig = None
        self._diff_cache: dict[tuple, tuple] = {}

    # ------------------------------------------------------------------ paths

    def study_dir(self, study_id: str) -> Path:
        return self.repo / "data" / "studies" / study_id

    def _signature(self) -> tuple:
        base = self.repo / "data" / "studies"
        files = sorted([*base.glob("*/study_summary.yaml"), *base.glob("*/parts/*/study_summary.yaml"),
                        *base.glob("*/genes/*.tsv.gz"), *base.glob("*/sample_curation.tsv"),
                        *base.glob("*/samplesheet.tsv"), *base.glob("*/qc_summary.tsv"),
                        *(self.repo / "registry" / "studies").glob("*.yaml")])
        return tuple((str(p), p.stat().st_mtime_ns) for p in files)

    # ------------------------------------------------------------------ registry and summaries

    def registry(self) -> dict[str, dict]:
        return load_studies(self.repo)

    def summaries(self) -> dict[str, list[dict]]:
        """study_id -> its study_summary.yaml entries (one per part)."""
        out: dict[str, list[dict]] = {}
        for e in summaries_on_main(self.repo):
            sid = e["summary"].get("study_id")
            if sid:
                out.setdefault(sid, []).append({**e["summary"], "_path": e["path"]})
        return out

    def status(self, study_id: str) -> dict:
        return _yaml(self.study_dir(study_id) / "status.yaml")

    def samplesheet(self, study_id: str) -> list[dict]:
        return read_tsv(self.study_dir(study_id) / "samplesheet.tsv")

    def qc(self, study_id: str) -> list[dict]:
        return read_tsv(self.study_dir(study_id) / "qc_summary.tsv")

    def curation(self, study_id: str) -> dict[str, dict]:
        return {r["sample_id"]: r for r in read_tsv(self.study_dir(study_id) / "sample_curation.tsv")}

    def replication(self, study_id: str) -> dict | None:
        path = self.study_dir(study_id) / "replication.json"
        return json.loads(path.read_text()) if path.exists() else None

    def atlas_json(self) -> dict | None:
        path = self.repo / "site" / "atlas" / "data" / "atlas.json"
        return json.loads(path.read_text()) if path.exists() else None

    @cached_property
    def gene_info(self) -> dict[str, dict]:
        assemblies = sorted((self.repo / "data" / "reference" / "gene_info").glob("*.tsv.gz"))
        if not assemblies:
            return {}
        # The atlas assembly, when the studies name one; else the only table there is.
        names = {s.get("genome_assembly") for es in self.summaries().values() for s in es}
        for name in names:
            path = gene_info_path(self.repo, name or "")
            if name and path.exists():
                return read_gene_info(path)
        return read_gene_info(assemblies[0])

    # ------------------------------------------------------------------ atlas

    def atlas(self) -> Atlas:
        with self._lock:
            sig = self._signature()
            if self._atlas is None or sig != self._atlas_sig:
                curation = {d.name: self.curation(d.name)
                            for d in sorted((self.repo / "data" / "studies").iterdir()) if d.is_dir()}
                self._atlas = Atlas(self.repo, curation)
                self._atlas_sig = sig
            return self._atlas

    # ------------------------------------------------------------------ differential

    def differential_dir(self, study_id: str) -> Path:
        return self.study_dir(study_id) / "differential"

    def differential_tests(self, study_id: str) -> list[dict]:
        return read_tsv(self.differential_dir(study_id) / "tests.tsv")

    def differential_summary(self, study_id: str) -> dict:
        return _yaml(self.differential_dir(study_id) / "summary.yaml")

    def differential_tables(self, study_id: str) -> list[tuple[str, str, Path]]:
        """(job_id, region, path) for each per-gene table."""
        out = []
        for p in sorted(self.differential_dir(study_id).glob("*.tsv.gz")):
            stem = p.name[: -len(".tsv.gz")]
            job, _, region = stem.rpartition(".")
            if region in REGIONS and job:
                out.append((job, region, p))
        return out

    def differential_table(self, path: Path) -> tuple[list[str], list[list[str]]]:
        """(header, rows) of one per-gene table; the last few are kept in memory."""
        key = (str(path), path.stat().st_mtime_ns)
        with self._lock:
            if key in self._diff_cache:
                return self._diff_cache[key]
        with gzip.open(path, "rt") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        value = (rows[0], rows[1:]) if rows else ([], [])
        with self._lock:
            if len(self._diff_cache) >= 6:
                self._diff_cache.pop(next(iter(self._diff_cache)))
            self._diff_cache[key] = value
        return value

    def differential_rows_for(self, path: Path, gene_ids: set[str]) -> tuple[list[str], list[list[str]]]:
        """The header and the rows for these genes, streamed (no caching)."""
        prefixes = tuple(g + "\t" for g in gene_ids)
        with gzip.open(path, "rt") as fh:
            header = next(csv.reader([fh.readline()], delimiter="\t"), [])
            rows = [next(csv.reader([line], delimiter="\t")) for line in fh if line.startswith(prefixes)]
        return header, rows


def nan_to_none(x, digits: int = 2):
    if x is None:
        return None
    x = float(x)
    return None if np.isnan(x) else round(x, digits)


__all__ = ["Store", "Atlas", "REGIONS", "GENE_MIN_CPGS", "nan_to_none", "parse_factors", "read_tsv"]
