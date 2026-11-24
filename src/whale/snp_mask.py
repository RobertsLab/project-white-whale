"""Register a C/T SNP mask built by workflows/snp_mask: `whale record-mask`.

The mask's files (cpg_mask.bed.gz, snps.tsv.gz, mask_manifest.json) go to the
data host under reference/<ASSEMBLY>/snp_mask/<MASK_ID>/, and an entry with
their sizes, SHA-256s and URLs, the studies and libraries the mask was called
from, its thresholds and its counts is written to data/reference/snp_masks.yaml.
That file is committed (in a normal PR); analyses name the mask they used by
its id. See docs/snp_mask.md.
"""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import yaml

from whale.checksums import sha256_file
from whale.results import DataHost, RecordError

MASK_FILES = ["cpg_mask.bed.gz", "snps.tsv.gz", "mask_manifest.json"]
REGISTRY = Path("data/reference/snp_masks.yaml")
HEADER = ("# C/T SNP masks (PLAN.md §6.1, docs/snp_mask.md), written by `whale record-mask`.\n"
          "# Each lists the reference CpGs on a SNP called from pooled WGBS libraries.\n")


def mask_entry(mask_dir: Path, url_base: str, today: str) -> dict:
    manifest = json.loads((mask_dir / "mask_manifest.json").read_text())
    p, c = manifest["params"], manifest["counts"]
    prefix = f"reference/{p['genome_assembly']}/snp_mask/{p['mask_id']}"
    for name, digest in manifest.get("outputs", {}).items():
        if sha256_file(mask_dir / name) != digest:
            raise RecordError(f"{name} does not match mask_manifest.json; was it changed after the run?")
    return {
        "mask_id": p["mask_id"],
        "genome_assembly": p["genome_assembly"],
        "recorded_on": today,
        "studies": p.get("studies", []),
        "libraries": len(p.get("samples", [])),
        "thresholds": {"min_alt_reads": p["min_alt_reads"], "min_alt_frac": p["min_alt_frac"],
                       "min_mapq": p.get("min_mapq"), "min_baseq": p.get("min_baseq")},
        "tool_versions": p.get("tool_versions", {}),
        "commit": p.get("commit"),
        "counts": {k: c[k] for k in ("snps", "snps_per_kb", "cpgs", "cpgs_masked", "cpgs_masked_frac")},
        "files": [{"file": f, "bytes": (mask_dir / f).stat().st_size, "sha256": sha256_file(mask_dir / f),
                   "url": f"{url_base.rstrip('/')}/{prefix}/{f}" if url_base else ""} for f in MASK_FILES],
    }


def record_mask(repo, mask_dir, *, upload: bool = True, upload_target: str | None = None,
                today: str | None = None) -> dict:
    repo, mask_dir = Path(repo), Path(mask_dir)
    missing = [f for f in MASK_FILES if not (mask_dir / f).exists()]
    if missing:
        raise RecordError(f"{mask_dir} lacks {', '.join(missing)}; did the snp_mask run finish?")
    host_path = repo / "data" / "data_host.yaml"
    host = yaml.safe_load(host_path.read_text()) if host_path.exists() else {}
    entry = mask_entry(mask_dir, host.get("url_base", ""), today or date.today().isoformat())
    target = upload_target or host.get("upload_target")
    if upload:
        if not target:
            raise RecordError("data/data_host.yaml has no upload_target; pass --upload-target or --no-upload")
        DataHost(target).put([mask_dir / f for f in MASK_FILES],
                             f"reference/{entry['genome_assembly']}/snp_mask/{entry['mask_id']}")
    path = repo / REGISTRY
    doc = (yaml.safe_load(path.read_text()) if path.exists() else None) or {"masks": []}
    doc["masks"] = [m for m in doc.get("masks") or [] if m["mask_id"] != entry["mask_id"]] + [entry]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + yaml.safe_dump(doc, sort_keys=False))
    return {"entry": entry, "uploaded": upload, "registry": path}
