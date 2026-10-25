"""Pool the merged studies into the atlas tables (PLAN.md §6.2) and compute the
first two cross-study views (§7.1): the methylation landscape and cross-study
concordance.

Inputs, all committed on main:

* data/studies/*/study_summary.yaml         which studies have a recorded run, their mask and assembly
* data/studies/*/genes/{class}.{metric}.tsv.gz  genes x libraries (meth, calls, cpgs)
* data/studies/*/samplesheet.tsv, qc_summary.tsv, registry/studies/*.yaml   sample metadata and QC

Only studies merged to main are pooled, so the atlas changes when a study PR
merges. A study is left out (and listed with the reason) when its gene tables
are not SNP-masked or are on another assembly.

Outputs (`write`): the pooled tables, genes x libraries for each class and
metric, `sample_metadata.tsv` and `sample_qc.tsv`, plus `atlas.json` for the
pages in site/atlas/. The tables are derived, so they are built in CI and
published with the site rather than committed.

Libraries of one biological sample (technical replicates: same BioSample,
unless the study shares one BioSample across libraries) are merged into one
unit, summing calls, before anything is compared, so replicates do not count
as independent animals.
"""
from __future__ import annotations

import csv
import gzip
import io
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from whale.dashboard import GENE_MIN_CPGS, HIST_BIN, _git, summaries_on_main
from whale.registry import load_studies

CLASSES = [("gene_body", "Gene body"), ("promoter", "Promoter")]
METRICS = ("meth", "calls", "cpgs")
LOW, HIGH = 10, 50          # % methylation: a gene under LOW is "unmethylated", at or over HIGH "methylated"
DENSITY_BINS = 25           # study-vs-study density grid, 4% cells

# sample_metadata.tsv columns, in order
METADATA_COLUMNS = ["library_id", "unit_id", "study_id", "part", "bioproject", "chemistry", "tissue",
                    "life_stage", "sample_unit", "pool_size", "ploidy", "sex", "condition", "factors",
                    "genome_assembly", "snp_mask", "nextflow_run"]


class AtlasError(RuntimeError):
    pass


# --------------------------------------------------------------------------- reading

def _tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))


def _matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """genes x libraries; an empty cell is 0."""
    with gzip.open(path, "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    libs, genes = rows[0][1:], [r[0] for r in rows[1:]]
    values = np.array([[int(v) if v else 0 for v in r[1:]] for r in rows[1:]], dtype=np.int64)
    return libs, genes, values.reshape(len(genes), len(libs))


def _single(values) -> str:
    values = [str(v) for v in values or [] if v not in (None, "")]
    return values[0] if len(values) == 1 else ""


def _sex(row: dict) -> str:
    v = (row.get("sex") or row.get("gender") or "").strip().lower()
    return v if v in ("female", "male", "hermaphrodite") else "unknown"


def library_metadata(entry: dict, study: dict, sheet: list[dict], libraries: list[str]) -> list[dict]:
    """One row per library in the gene tables (those that passed the QC gate)."""
    s = entry["summary"]
    samples = study.get("samples") or {}
    shared = "shared_biosample" in (study.get("quality_flags") or [])
    # factor id -> its sample-sheet column, as record-run resolved it
    columns = {f["factor_id"]: f.get("column") or f["factor_id"] for f in s.get("by_factor") or []}
    factors = [(f["factor_id"], columns.get(f["factor_id"], f["factor_id"])) for f in study.get("factors") or []]
    by_id = {r["sample_id"]: r for r in sheet}
    out = []
    for lib in libraries:
        r = by_id.get(lib, {})
        unit = lib if shared or not r.get("sample_accession") else r["sample_accession"]
        out.append({
            "library_id": lib,
            "unit_id": unit,
            "study_id": s["study_id"],
            "part": s.get("part") or "",
            "bioproject": (study.get("accessions") or {}).get("bioproject") or "",
            "chemistry": s.get("chemistry") or study.get("chemistry") or "",
            # Registry terms when the study has one; the deposit's own word otherwise.
            "tissue": _single(samples.get("tissues")) or r.get("tissue", ""),
            "life_stage": _single(samples.get("life_stages")) or r.get("life_stage", ""),
            "sample_unit": samples.get("sample_unit") or "",
            "pool_size": samples.get("pool_size") or "",
            "ploidy": r.get("ploidy") or _single(samples.get("ploidy")),
            "sex": _sex(r),
            "condition": r.get("condition", ""),
            "factors": ";".join(f"{f}={r[col]}" for f, col in factors if r.get(col, "") != ""),
            "genome_assembly": s.get("genome_assembly") or "",
            "snp_mask": (s.get("snp_mask") or {}).get("mask_id") or "",
            "nextflow_run": s.get("nextflow_run") or "",
        })
    return out


def pool(repo: Path) -> dict:
    """The merged studies' gene tables side by side, with sample metadata and QC.

    Returns {"genes", "libraries" (metadata rows), "qc" (rows), "matrices"
    {(class, metric): genes x libraries}, "studies" (included), "excluded"}."""
    repo = Path(repo)
    registered = load_studies(repo)
    entries = summaries_on_main(repo)
    assemblies = [e["summary"].get("genome_assembly") for e in entries if e["summary"].get("genome_assembly")]
    assembly = max(set(assemblies), key=assemblies.count) if assemblies else None

    genes, libraries, qc, studies, excluded = None, [], [], [], []
    blocks: dict[tuple[str, str], list[np.ndarray]] = {(c, m): [] for c, _ in CLASSES for m in METRICS}
    for e in entries:
        s = e["summary"]
        sid = s.get("study_id")
        base = repo / e["path"].rsplit("/", 1)[0]
        paths = {(c, m): base / "genes" / f"{c}.{m}.tsv.gz" for c, _ in CLASSES for m in METRICS}
        reason = None
        if not all(p.exists() for p in paths.values()):
            reason = "no gene tables (recorded before gene-level summaries)"
        elif not (s.get("snp_mask") or {}).get("mask_id"):
            reason = "gene tables not SNP-masked"
        elif s.get("genome_assembly") != assembly:
            reason = f"on {s.get('genome_assembly')}, not {assembly}"
        elif sid not in registered:
            reason = "not registered on main"
        if reason:
            excluded.append({"study_id": sid, "part": s.get("part"), "reason": reason})
            continue

        libs = None
        for key, path in paths.items():
            lib_ids, gene_ids, values = _matrix(path)
            if genes is None:
                genes = gene_ids
            if gene_ids != genes:
                raise AtlasError(f"{path}: genes differ from the other studies' (same assembly {assembly})")
            if libs is not None and lib_ids != libs:
                raise AtlasError(f"{path}: libraries differ from the study's other gene tables")
            libs = lib_ids
            blocks[key].append(values)
        sheet = _tsv(base / "samplesheet.tsv") or _tsv(repo / "data" / "studies" / sid / "samplesheet.tsv")
        meta = library_metadata(e, registered[sid], sheet, libs)
        libraries += meta
        qc += [{"study_id": sid, "part": s.get("part") or "", **r} for r in _tsv(base / "qc_summary.tsv")]
        studies.append({"study_id": sid, "part": s.get("part"), "libraries": len(libs),
                        "units": len({m["unit_id"] for m in meta}),
                        "tissues": sorted({m["tissue"] for m in meta}),
                        "life_stages": sorted({m["life_stage"] for m in meta}),
                        "chemistry": meta[0]["chemistry"] if meta else "",
                        "sample_unit": meta[0]["sample_unit"] if meta else "",
                        "snp_mask": meta[0]["snp_mask"] if meta else "",
                        "paired_rnaseq": (registered[sid].get("paired_rnaseq") or {}).get("relationship")})

    ids = [m["library_id"] for m in libraries]
    if len(set(ids)) != len(ids):
        raise AtlasError("a library id appears in more than one study")
    n = len(genes or [])
    matrices = {k: (np.hstack(v) if v else np.zeros((n, 0), dtype=np.int64)) for k, v in blocks.items()}
    return {"genome_assembly": assembly, "genes": genes or [], "libraries": libraries, "qc": qc,
            "matrices": matrices, "studies": studies, "excluded": excluded}


# --------------------------------------------------------------------------- units

def units(pooled: dict) -> dict:
    """Libraries merged per biological unit: methylated and total calls summed; CpGs
    with a call is the largest of the libraries' (the union is not in the tables)."""
    libs = pooled["libraries"]
    order: dict[tuple[str, str], int] = {}
    for m in libs:
        order.setdefault((m["study_id"], m["unit_id"]), len(order))
    idx = np.array([order[(m["study_id"], m["unit_id"])] for m in libs], dtype=np.int64)
    n_units, n_genes = len(order), len(pooled["genes"])
    out = {}
    for c, _ in CLASSES:
        for metric in METRICS:
            src = pooled["matrices"][(c, metric)]
            dst = np.zeros((n_genes, n_units), dtype=np.int64)
            for j in range(src.shape[1]):
                if metric == "cpgs":
                    np.maximum(dst[:, idx[j]], src[:, j], out=dst[:, idx[j]])
                else:
                    dst[:, idx[j]] += src[:, j]
            out[(c, metric)] = dst
    qc_by_lib = {r["sample_id"]: r for r in pooled["qc"]}
    meta = []
    for (sid, uid), k in order.items():
        members = [m for m, i in zip(libs, idx, strict=True) if i == k]
        first = members[0]
        cpg = [float(qc_by_lib[m["library_id"]]["cpg_meth_pct"]) for m in members
               if qc_by_lib.get(m["library_id"], {}).get("cpg_meth_pct")]
        meta.append({"unit_id": uid, "study_id": sid, "libraries": [m["library_id"] for m in members],
                     **{k2: first[k2] for k2 in ("tissue", "life_stage", "sample_unit", "ploidy", "sex",
                                                 "condition", "factors")},
                     "cpg_meth_pct": round(sum(cpg) / len(cpg), 2) if cpg else None})
    return {"units": meta, "matrices": out}


def _pct(meth: np.ndarray, calls: np.ndarray, cpgs: np.ndarray, min_cpgs: int = GENE_MIN_CPGS) -> np.ndarray:
    """% methylated calls; NaN where fewer than min_cpgs CpGs have a call."""
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = 100.0 * meth / calls
    pct[(cpgs < min_cpgs) | (calls == 0)] = np.nan
    return pct


# --------------------------------------------------------------------------- statistics

def _rank(x: np.ndarray) -> np.ndarray:
    """Ranks with ties averaged (many genes sit at exactly 0%)."""
    order = np.argsort(x, kind="mergesort")
    _, start, counts = np.unique(x[order], return_index=True, return_counts=True)
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.repeat(start + (counts - 1) / 2.0, counts)
    return ranks


def spearman(a: np.ndarray, b: np.ndarray) -> tuple[float | None, int]:
    """Spearman's rho over the genes measured in both, and how many there are."""
    ok = ~np.isnan(a) & ~np.isnan(b)
    n = int(ok.sum())
    if n < 3:
        return None, n
    ra, rb = _rank(a[ok]), _rank(b[ok])
    if ra.std() == 0 or rb.std() == 0:
        return None, n
    return float(np.corrcoef(ra, rb)[0, 1]), n


def _r(x, d=3):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), d)


def _shape(v: np.ndarray) -> dict:
    """Distribution of one sample's (or study's) gene values."""
    v = v[~np.isnan(v)]
    if not len(v):
        return {"genes": 0}
    bins = np.minimum((v // HIST_BIN).astype(int), 100 // HIST_BIN - 1)
    return {"genes": int(len(v)), "median": _r(np.median(v), 2), "mean": _r(v.mean(), 2),
            "low_frac": _r((v < LOW).mean()), "high_frac": _r((v >= HIGH).mean()),
            "mid_frac": _r(((v >= LOW) & (v < HIGH)).mean()),
            "bins": np.bincount(bins, minlength=100 // HIST_BIN).tolist()}


def _median(xs):
    xs = [x for x in xs if x is not None]
    return _r(float(np.median(xs))) if xs else None


def landscape(pooled: dict, u: dict) -> dict:
    """Per unit and per study: how gene-body and promoter methylation are distributed."""
    per_unit = []
    for c, _ in CLASSES:
        pct = _pct(u["matrices"][(c, "meth")], u["matrices"][(c, "calls")], u["matrices"][(c, "cpgs")])
        for j, m in enumerate(u["units"]):
            if c == CLASSES[0][0]:
                per_unit.append({"unit_id": m["unit_id"], "study_id": m["study_id"],
                                 "libraries": len(m["libraries"]), "tissue": m["tissue"],
                                 "condition": m["condition"], "cpg_meth_pct": m["cpg_meth_pct"]})
            per_unit[j][c] = _shape(pct[:, j])
    studies = []
    for st in pooled["studies"]:
        rows = [p for p in per_unit if p["study_id"] == st["study_id"]]
        row = {"study_id": st["study_id"], "units": len(rows), "tissues": st["tissues"]}
        for c, _ in CLASSES:
            row[c] = {k: _median([p[c].get(k) for p in rows]) for k in ("median", "low_frac", "high_frac", "mid_frac")}
            row[c]["genes"] = _median([p[c].get("genes") for p in rows])
        row["cpg_meth_pct"] = _median([p["cpg_meth_pct"] for p in rows])
        studies.append(row)
    return {"units": per_unit, "studies": studies}


def _study_means(pct: np.ndarray, cols: list[int]) -> np.ndarray:
    """A study's value per gene: the mean over its units with the gene measured, kept
    when at least half of its units measure it."""
    sub = pct[:, cols]
    ok = ~np.isnan(sub)
    n_ok = ok.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(ok, sub, 0.0).sum(axis=1) / n_ok
    mean[n_ok < max(1, (len(cols) + 1) // 2)] = np.nan
    return mean


def concordance(pooled: dict, u: dict) -> dict:
    """Rank agreement of gene-body methylation between units and between studies."""
    c = "gene_body"
    pct = _pct(u["matrices"][(c, "meth")], u["matrices"][(c, "calls")], u["matrices"][(c, "cpgs")])
    meta = u["units"]
    n = len(meta)
    rho = [[None] * n for _ in range(n)]
    n_genes = [[0] * n for _ in range(n)]
    for i in range(n):
        rho[i][i], n_genes[i][i] = 1.0, int((~np.isnan(pct[:, i])).sum())
        for j in range(i + 1, n):
            r, k = spearman(pct[:, i], pct[:, j])
            rho[i][j] = rho[j][i] = _r(r)
            n_genes[i][j] = n_genes[j][i] = k

    sids = [st["study_id"] for st in pooled["studies"]]
    cols = {sid: [j for j, m in enumerate(meta) if m["study_id"] == sid] for sid in sids}
    blocks = []
    for a_i, a in enumerate(sids):
        for b in sids[a_i:]:
            vals = [rho[i][j] for i in cols[a] for j in cols[b] if i != j and (a != b or i < j)
                    and rho[i][j] is not None]
            blocks.append({"a": a, "b": b, "pairs": len(vals),
                           "mean_rho": _r(np.mean(vals)) if vals else None,
                           "min_rho": _r(min(vals)) if vals else None,
                           "max_rho": _r(max(vals)) if vals else None})

    means = {sid: _study_means(pct, cols[sid]) for sid in sids}
    tissues = {st["study_id"]: set(st["tissues"]) for st in pooled["studies"]}
    pairs = []
    for a_i, a in enumerate(sids):
        for b in sids[a_i + 1:]:
            x, y = means[a], means[b]
            r, k = spearman(x, y)
            ok = ~np.isnan(x) & ~np.isnan(y)
            xa, ya = x[ok], y[ok]
            step = 100 / DENSITY_BINS
            grid = np.zeros((DENSITY_BINS, DENSITY_BINS), dtype=int)
            np.add.at(grid, (np.minimum((ya // step).astype(int), DENSITY_BINS - 1),
                             np.minimum((xa // step).astype(int), DENSITY_BINS - 1)), 1)
            both_low = int(((xa < LOW) & (ya < LOW)).sum())
            both_high = int(((xa >= HIGH) & (ya >= HIGH)).sum())
            opposite = int((((xa < LOW) & (ya >= HIGH)) | ((xa >= HIGH) & (ya < LOW))).sum())
            pairs.append({"a": a, "b": b, "rho": _r(r), "genes": k,
                          "shared_tissue": sorted(tissues[a] & tissues[b]),
                          "both_low": both_low, "both_high": both_high, "opposite": opposite,
                          "other": k - both_low - both_high - opposite,
                          "mean_abs_diff": _r(np.abs(xa - ya).mean(), 2) if k else None,
                          "density": grid.tolist()})
    return {"class": c, "units": [{"unit_id": m["unit_id"], "study_id": m["study_id"], "tissue": m["tissue"],
                                   "condition": m["condition"]} for m in meta],
            "rho": rho, "min_genes": min((n_genes[i][j] for i in range(n) for j in range(n) if i != j), default=0),
            "blocks": blocks, "study_pairs": pairs}


def design(pooled: dict) -> dict:
    """Which studies each tissue, life stage, chemistry and sample unit rests on, and
    the comparisons that cannot be separated from study."""
    studies = pooled["studies"]
    out, caveats = {}, []
    for key, label in (("tissues", "tissue"), ("life_stages", "life stage"), ("chemistry", "chemistry"),
                       ("sample_unit", "sample unit")):
        levels: dict[str, list[str]] = {}
        for st in studies:
            for v in st[key] if isinstance(st[key], list) else [st[key]]:
                levels.setdefault(v or "unrecorded", []).append(st["study_id"])
        out[key] = [{"level": k, "studies": v} for k, v in sorted(levels.items())]
        if len(levels) > 1:
            for k, v in sorted(levels.items()):
                if len(v) == 1:
                    caveats.append(f"{label.capitalize()} '{k}' comes from one study ({v[0]}), so its "
                                   f"difference from other {label}s cannot be separated from study.")
        elif len(studies) > 1:
            caveats.append(f"Every study has the same {label} ({next(iter(levels))}); nothing here speaks to "
                           f"{label} differences.")
    return {"levels": out, "caveats": caveats}


# --------------------------------------------------------------------------- build / write

def headline(data: dict) -> list[str]:
    """Plain statements, each with the number of studies and samples it rests on."""
    out = []
    n_st, n_u = len(data["studies"]), sum(s["units"] for s in data["studies"])
    if not n_st:
        return out
    ls = data["landscape"]["studies"]
    gb = [s["gene_body"]["high_frac"] for s in ls if s["gene_body"]["high_frac"] is not None]
    lo = [s["gene_body"]["low_frac"] for s in ls if s["gene_body"]["low_frac"] is not None]
    if gb:
        out.append(f"Across {n_st} {'study' if n_st == 1 else 'studies'} ({n_u} samples), a median "
                   f"{100 * float(np.median(gb)):.0f}% of genes have gene-body methylation >= {HIGH}% and "
                   f"{100 * float(np.median(lo)):.0f}% under {LOW}% (study medians "
                   f"{100 * min(gb):.0f}-{100 * max(gb):.0f}% methylated).")
    blocks = data["concordance"]["blocks"]
    within = [b["mean_rho"] for b in blocks if b["a"] == b["b"] and b["mean_rho"] is not None]
    between = [b["mean_rho"] for b in blocks if b["a"] != b["b"] and b["mean_rho"] is not None]
    if within and between:
        out.append(f"Rank agreement of gene-body methylation (mean Spearman rho between samples): "
                   f"{min(within):.2f}-{max(within):.2f} within a study, {min(between):.2f}-{max(between):.2f} "
                   f"between studies ({n_st} studies, {n_u} samples).")
    return out


def build(repo: Path, *, repo_url: str = "", now: str | None = None) -> tuple[dict, dict]:
    """(atlas.json data, pooled tables)."""
    pooled = pool(repo)
    u = units(pooled)
    data = {
        "generated_at": now or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "commit": _git(Path(repo), "rev-parse", "--short", "HEAD"),
        "repo_url": repo_url,
        "genome_assembly": pooled["genome_assembly"],
        "genes": len(pooled["genes"]),
        "libraries": len(pooled["libraries"]),
        "min_cpgs": GENE_MIN_CPGS, "bin": HIST_BIN, "low": LOW, "high": HIGH,
        "studies": pooled["studies"],
        "excluded": pooled["excluded"],
        "design": design(pooled),
        "landscape": landscape(pooled, u),
        "concordance": concordance(pooled, u),
        "tables": [f"tables/{c}.{m}.tsv.gz" for c, _ in CLASSES for m in (*METRICS, "meth_pct")]
                  + ["tables/sample_metadata.tsv", "tables/sample_qc.tsv"],
    }
    data["headline"] = headline(data)
    return data, pooled


def _write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _write_matrix(path: Path, genes: list[str], libs: list[str], values: np.ndarray, fmt) -> None:
    # mtime=0: the same tables give the same bytes.
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="\n")
        text.write("\t".join(["gene_id", *libs]) + "\n")
        for g, row in zip(genes, values, strict=True):
            text.write(g + "\t" + "\t".join(fmt(v) for v in row) + "\n")
        text.flush()
        text.detach()


def write(data: dict, pooled: dict, out_dir: Path) -> Path:
    """atlas.json and tables/ under out_dir (site/atlas/data by default)."""
    import json

    out_dir = Path(out_dir)
    tables = out_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    libs = [m["library_id"] for m in pooled["libraries"]]
    for c, _ in CLASSES:
        for metric in METRICS:
            _write_matrix(tables / f"{c}.{metric}.tsv.gz", pooled["genes"], libs,
                          pooled["matrices"][(c, metric)], str)
        mx = pooled["matrices"]
        with np.errstate(invalid="ignore", divide="ignore"):
            pct = 100.0 * mx[(c, "meth")] / mx[(c, "calls")]
        _write_matrix(tables / f"{c}.meth_pct.tsv.gz", pooled["genes"], libs, pct,
                      lambda v: "" if np.isnan(v) else f"{v:.3f}")
    _write_tsv(tables / "sample_metadata.tsv", METADATA_COLUMNS, pooled["libraries"])
    qc_cols = list(dict.fromkeys(k for r in pooled["qc"] for k in r)) or ["study_id"]
    _write_tsv(tables / "sample_qc.tsv", qc_cols, pooled["qc"])
    path = out_dir / "atlas.json"
    path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    return path


def summary_text(data: dict) -> str:
    lines = [f"{len(data['studies'])} studies, {data['libraries']} libraries, "
             f"{sum(s['units'] for s in data['studies'])} samples, {data['genes']:,} genes"]
    lines += [f"  excluded {e['study_id']}: {e['reason']}" for e in data["excluded"]]
    lines += [f"  {h}" for h in data["headline"]]
    return "\n".join(lines)

