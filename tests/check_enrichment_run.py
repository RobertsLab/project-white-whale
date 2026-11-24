"""Check a real (not stub) enrichment-workflow run on the fixtures against how they were simulated.

    python3 tests/check_enrichment_run.py results/test/TEST2026_FIXTURE_MEDIP

The fixture reads (tests/fixtures/make_enrichment_fixtures.py) are plain reads
captured in proportion to methylated CpGs, and only gene-body CpGs are
methylated (mosaic, as in the oyster). A run that aligns, deduplicates and
counts them correctly maps nearly every read, removes the ~5% duplicates,
finds every fixture gene body enriched (obs_exp > 1) and gene bodies far
above intergenic DNA. Exits non-zero, listing every failure.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

SAMPLES = ["E1", "E2"]
GENES = ["G1", "G2", "G3", "G4"]
MIN_MAPPING = 95.0
DUP_RANGE = (1.0, 15.0)
MIN_FRAGMENTS = 2000
MIN_GENE_OVER_INTERGENIC = 3.0
TILES = 30                  # chr1 20 kb + chr2 10 kb in 1 kb tiles


def _tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check(run: Path) -> list[str]:
    problems = []
    qc = {r["sample_id"]: r for r in _tsv(run / "qc_summary.tsv")}
    for sid in SAMPLES:
        r = qc.get(sid)
        if r is None:
            problems.append(f"{sid}: not in qc_summary.tsv")
            continue
        f = {k: float(r[k]) if r[k] else None for k in ("mapping_pct", "duplication_pct", "fragments", "relh",
                                                         "saturation_r", "gene_body_obs_exp", "intergenic_obs_exp")}
        if f["mapping_pct"] is None or f["mapping_pct"] < MIN_MAPPING:
            problems.append(f"{sid}: mapping {f['mapping_pct']}% < {MIN_MAPPING}%")
        if f["duplication_pct"] is None or not DUP_RANGE[0] <= f["duplication_pct"] <= DUP_RANGE[1]:
            problems.append(f"{sid}: duplication {f['duplication_pct']}% outside {DUP_RANGE}")
        if f["fragments"] is None or f["fragments"] < MIN_FRAGMENTS:
            problems.append(f"{sid}: {f['fragments']} fragments < {MIN_FRAGMENTS}")
        if f["relh"] is None or f["relh"] <= 1.0:
            problems.append(f"{sid}: relH {f['relh']} <= 1 (captured fragments not CpG-rich)")
        if f["saturation_r"] is None or f["saturation_r"] < 0.5:
            problems.append(f"{sid}: saturation r {f['saturation_r']} < 0.5")
        gb, ig = f["gene_body_obs_exp"], f["intergenic_obs_exp"]
        if gb is None or ig is None or gb < MIN_GENE_OVER_INTERGENIC * ig:
            problems.append(f"{sid}: gene bodies {gb} not >= {MIN_GENE_OVER_INTERGENIC}x intergenic {ig}")
        if r["qc_pass"] != "true":
            problems.append(f"{sid}: failed the QC gate ({r['qc_reasons']})")

        feats = run / "features" / f"{sid}.enrichment_features.tsv.gz"
        if not feats.exists():
            problems.append(f"{sid}: {feats} is missing")
        else:
            body = {x["gene_id"]: x for x in _tsv(feats) if x["feature"] == "gene_body"}
            if list(body) != GENES:
                problems.append(f"{sid}: gene bodies {list(body)} != {GENES}")
            for g in GENES:
                v = body.get(g, {}).get("obs_exp")
                if not v or float(v) <= 1.0:
                    problems.append(f"{sid}: gene body {g} obs_exp {v!r} <= 1 (it is methylated)")
        tiles = run / "tiles" / f"{sid}.tiles_1kb.tsv.gz"
        if not tiles.exists():
            problems.append(f"{sid}: {tiles} is missing")
        elif len(_tsv(tiles)) != TILES:
            problems.append(f"{sid}: {len(_tsv(tiles))} tiles, expected {TILES}")

    manifest = json.loads((run / "run_manifest.json").read_text())
    p = manifest.get("params", {})
    if p.get("tier") != "enrichment" or not (p.get("enrichment_method") or {}).get("name"):
        problems.append("run_manifest.json does not record tier: enrichment and the enrichment method")
    if "BOWTIE2_ALIGN" not in manifest.get("tool_versions", {}):
        problems.append("run_manifest.json has no bowtie2 version")
    return problems


def main() -> int:
    problems = check(Path(sys.argv[1]))
    for p in problems:
        print(f"FAIL {p}")
    if not problems:
        print("ok: enrichment run matches the simulated MeDIP fixtures")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
