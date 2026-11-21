"""The read-only questions `whale mcp` answers, as plain functions returning JSON-ready dicts.

Every answer names where its numbers come from: a `sources` list with each study's
study_id, Nextflow run name, SNP mask and genome assembly. Long answers are
paginated (`limit`, `offset`) under a hard cap that is reported with them.
"""
from __future__ import annotations

import re
import warnings
from contextlib import contextmanager

import numpy as np

from whale.atlas import _study_means
from whale.query.store import GENE_MIN_CPGS, REGIONS, Store, nan_to_none

MAX_LIMIT = 200          # rows per page for any list
MAX_GENES_DETAILED = 10  # genes per get_gene call
LAMBDA_WARN = 1.3        # genomic inflation above which significant counts are not trusted (docs/replication.md)
VERDICTS = ("replicated", "partly", "not_replicated", "inconclusive", "not_testable", "pending")


class QueryError(ValueError):
    """A question the data cannot answer as asked (unknown study, bad region, ...)."""


# --------------------------------------------------------------------------- helpers

def _page(rows: list, limit: int, offset: int, cap: int = MAX_LIMIT) -> dict:
    limit = max(1, min(int(limit), cap))
    offset = max(0, int(offset))
    chunk = rows[offset: offset + limit]
    end = offset + len(chunk)
    return {"rows": chunk, "page": {"total": len(rows), "offset": offset, "returned": len(chunk),
                                    "limit": limit, "cap": cap, "truncated": end < len(rows),
                                    "next_offset": end if end < len(rows) else None}}


def _sources(store: Store, study_ids) -> list[dict]:
    """Where a study's numbers come from: its recorded run, SNP mask and assembly."""
    summaries = store.summaries()
    out = []
    for sid in dict.fromkeys(study_ids):
        for s in summaries.get(sid) or [{}]:
            out.append({"study_id": sid, "part": s.get("part"),
                        "nextflow_run": s.get("nextflow_run"),
                        "snp_mask": (s.get("snp_mask") or {}).get("mask_id"),
                        "genome_assembly": s.get("genome_assembly"),
                        "summary_file": s.get("_path")})
    return out


def _check_study(store: Store, study_id: str) -> None:
    known = set(store.registry()) | set(store.summaries())
    if study_id not in known:
        raise QueryError(f"Unknown study_id {study_id!r}. Known: {', '.join(sorted(known))}")


def _region(region: str) -> str:
    if region not in REGIONS:
        raise QueryError(f"region must be one of {', '.join(REGIONS)}, not {region!r}")
    return region


def _gene_record(store: Store, gene_id: str) -> dict:
    gi = store.gene_info.get(gene_id, {})
    rec = {"gene_id": gene_id, "symbol": gi.get("symbol", ""), "description": gi.get("description", "")}
    if gi.get("chrom"):
        rec["location"] = f"{gi['chrom']}:{int(gi['start']) + 1}-{gi['end']}({gi.get('strand', '')})"
    return rec


def _match_genes(store: Store, query: str) -> list[str]:
    """GeneIDs matching an NCBI GeneID, a symbol (exact, case-insensitive) or, failing those,
    every word of a product description."""
    q = (query or "").strip()
    if not q:
        raise QueryError("give a GeneID, a symbol or words from a product description")
    info = store.gene_info
    gid = re.sub(r"^(GeneID:|gene-)", "", q, flags=re.IGNORECASE)
    if gid in info:
        return [gid]
    if not info:
        idx = store.atlas().gene_index
        return [gid] if gid in idx else []
    low = q.lower()
    exact = [g for g, r in info.items() if r.get("symbol", "").lower() == low]
    if exact:
        return exact
    words = low.split()
    return [g for g, r in info.items()
            if all(w in (r.get("description", "") + " " + r.get("symbol", "")).lower() for w in words)]


def _norm(v: str) -> str:
    return (v or "").strip().lower()


def _select_units(store: Store, *, tissue=None, life_stage=None, chemistry=None, study_ids=None):
    """Unit columns that pass the filters, per study, and which filters left units out
    because the attribute is not recorded."""
    at = store.atlas()
    wanted = set(study_ids) if study_ids else None
    if wanted:
        unknown = wanted - set(at.study_ids)
        if unknown:
            raise QueryError(f"Not in the atlas: {', '.join(sorted(unknown))}. In the atlas: "
                             f"{', '.join(at.study_ids)}")
    cols: dict[str, list[int]] = {}
    unrecorded: dict[str, int] = {}
    for sid in at.study_ids:
        if wanted and sid not in wanted:
            continue
        for j in at.columns[sid]:
            m = at.units[j]
            ok = True
            for key, val in (("tissue", tissue), ("life_stage", life_stage), ("chemistry", chemistry)):
                if val is None:
                    continue
                if not m.get(key):
                    unrecorded[sid] = unrecorded.get(sid, 0) + 1
                    ok = False
                    break
                if _norm(m[key]) != _norm(val):
                    ok = False
                    break
            if ok:
                cols.setdefault(sid, []).append(j)
    notes = [f"{sid}: {n} unit(s) left out because the filtered attribute is not recorded"
             for sid, n in unrecorded.items()]
    return cols, notes


def _study_matrix(store: Store, region: str, cols: dict[str, list[int]]):
    """genes x studies: each study's mean over its selected units (kept when at least half
    of them measure the gene), and the unit counts."""
    at = store.atlas()
    pct = at.pct[region]
    sids = list(cols)
    means = np.full((len(at.genes), len(sids)), np.nan)
    within_sd = np.full((len(at.genes), len(sids)), np.nan)
    for k, sid in enumerate(sids):
        means[:, k] = _study_means(pct, cols[sid])
        sub = pct[:, cols[sid]]
        if sub.shape[1] > 1:
            with _quiet():
                within_sd[:, k] = np.nanstd(sub, axis=1, ddof=1)
    return sids, means, within_sd


@contextmanager
def _quiet():
    """nanmean and friends warn on all-NaN rows; those rows are expected (gene not measured)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def _filters(**kw) -> dict:
    return {k: v for k, v in kw.items() if v not in (None, [], "")}


# --------------------------------------------------------------------------- studies

def overview(store: Store) -> dict:
    """What is in the resource and which outputs exist."""
    reg, summ = store.registry(), store.summaries()
    at = store.atlas()
    studies = sorted(set(reg) | set(summ))
    have = {k: [s for s in studies if (store.study_dir(s) / f).exists()]
            for k, f in (("differential", "differential/tests.tsv"),
                         ("sample_curation", "sample_curation.tsv"),
                         ("replication", "replication.json"))}
    out = {
        "resource": "White Whale: public Pacific oyster (Magallana gigas) DNA-methylation studies "
                    "re-processed with one pipeline on one reference and pooled into an atlas.",
        "genome_assembly": at.genome_assembly,
        "registered_studies": len(reg),
        "studies_with_recorded_run": len(summ),
        "atlas": {"studies": at.study_ids, "units": len(at.units), "libraries": len(at.libraries),
                  "genes": len(at.genes), "min_cpgs_per_gene": GENE_MIN_CPGS,
                  "excluded": at.excluded,
                  "unit_note": "libraries of one BioSample are summed into one unit before comparison"},
        "outputs_present": have,
        "outputs_missing_note": "Studies not listed under an output do not have it yet: differential "
                                "comes from `whale differential`, sample_curation from "
                                "data/studies/<ID>/sample_curation.tsv, replication from `whale replicate`.",
        "sources": _sources(store, at.study_ids),
    }
    built = store.atlas_json()
    if built:
        out["atlas"]["headline"] = built.get("headline")
        out["atlas"]["design_caveats"] = (built.get("design") or {}).get("caveats")
        out["atlas"]["built_at"] = built.get("generated_at")
    else:
        from whale.atlas import design
        out["atlas"]["design_caveats"] = design({"studies": at.studies})["caveats"]
    return out


def _study_row(store: Store, sid: str, study: dict, summaries: list[dict]) -> dict:
    s = study.get("samples") or {}
    summ = summaries[0] if summaries else {}
    libs = sum(x.get("libraries") or 0 for x in summaries)
    passed = sum(x.get("passed") or 0 for x in summaries)
    return {
        "study_id": sid,
        "bioproject": (study.get("accessions") or {}).get("bioproject"),
        "chemistry": study.get("chemistry") or summ.get("chemistry"),
        "tissues": s.get("tissues") or [],
        "life_stages": s.get("life_stages") or [],
        "sample_unit": s.get("sample_unit"),
        "n_runs": s.get("n_runs"),
        "n_libraries": s.get("n_libraries"),
        "factors": [f.get("factor_id") for f in study.get("factors") or []],
        "paired_rnaseq": (study.get("paired_rnaseq") or {}).get("relationship"),
        "quality_flags": study.get("quality_flags") or [],
        "stage": store.status(sid).get("stage"),
        "qc": ({"libraries": libs, "passed": passed,
                "pass_rate": round(passed / libs, 3) if libs else None} if summaries else None),
        "nextflow_run": summ.get("nextflow_run"),
        "snp_mask": (summ.get("snp_mask") or {}).get("mask_id"),
    }


def list_studies(store: Store, tissue=None, life_stage=None, chemistry=None, stage=None) -> dict:
    reg, summ = store.registry(), store.summaries()
    rows = []
    for sid in sorted(set(reg) | set(summ)):
        row = _study_row(store, sid, reg.get(sid, {}), summ.get(sid, []))
        if tissue and _norm(tissue) not in map(_norm, row["tissues"]):
            continue
        if life_stage and _norm(life_stage) not in map(_norm, row["life_stages"]):
            continue
        if chemistry and _norm(chemistry) != _norm(row["chemistry"]):
            continue
        if stage and _norm(stage) != _norm(row["stage"]):
            continue
        rows.append(row)
    return {"filters": _filters(tissue=tissue, life_stage=life_stage, chemistry=chemistry, stage=stage),
            "studies": rows, "count": len(rows),
            "note": "tissues/life_stages are the registry's study-level terms; per-sample values "
                    "are in get_samples.",
            "sources": _sources(store, [r["study_id"] for r in rows if r["qc"]])}


def describe_study(store: Store, study_id: str) -> dict:
    _check_study(store, study_id)
    study = store.registry().get(study_id, {})
    summaries = store.summaries().get(study_id, [])
    row = _study_row(store, study_id, study, summaries)
    at = store.atlas()
    excluded = [e for e in at.excluded if e["study_id"] == study_id]
    runs = []
    for s in summaries:
        runs.append({k: s.get(k) for k in ("part", "nextflow_run", "generated_on", "genome_assembly",
                                           "libraries", "passed", "failed", "qc", "cpg_methylation_pct",
                                           "by_factor", "by_condition", "genes")}
                    | {"snp_mask": s.get("snp_mask"),
                       "methylomes_url": (s.get("methylomes") or {}).get("url")})
    rep = store.replication(study_id)
    diff = store.differential_summary(study_id)
    return {
        **row,
        "citation": study.get("citation"),
        "accessions": study.get("accessions"),
        "population": study.get("population"),
        "samples": study.get("samples"),
        "factors": study.get("factors"),
        "paired_rnaseq": study.get("paired_rnaseq"),
        "limitations": study.get("limitations"),
        "data_availability": study.get("data_availability"),
        "runs": runs,
        "atlas": {"included": study_id in at.study_ids,
                  "units": len(at.columns.get(study_id, [])),
                  "excluded_reason": excluded[0]["reason"] if excluded else None},
        "available": {
            "differential": bool(store.differential_tests(study_id)),
            "differential_generated_on": diff.get("generated_on"),
            "sample_curation": bool(store.curation(study_id)),
            "replication": ({"counts": rep.get("counts"), "caveats": rep.get("caveats")} if rep else None),
        },
        "sources": _sources(store, [study_id]),
    }


SAMPLE_COLUMNS = ("sample_id", "sample_accession", "condition", "label")
QC_COLUMNS = ("qc_pass", "qc_reasons", "mapping_pct", "conversion_pct", "cpg_meth_pct", "cpgs_5x", "mean_depth")


def get_samples(store: Store, study_id: str, limit: int = 100, offset: int = 0, full: bool = False) -> dict:
    """One row per library: deposit fields, factor levels, curated attributes and QC."""
    _check_study(store, study_id)
    study = store.registry().get(study_id, {})
    summaries = store.summaries().get(study_id, [])
    columns = {f["factor_id"]: f.get("column") or f["factor_id"]
               for s in summaries for f in s.get("by_factor") or []}
    for f in study.get("factors") or []:
        columns.setdefault(f["factor_id"], f["factor_id"])
    sheet = store.samplesheet(study_id)
    qc = {r["sample_id"]: r for r in store.qc(study_id)}
    cur = store.curation(study_id)
    at = store.atlas()
    unit_of = {lib: m["unit_id"] for m in at.units if m["study_id"] == study_id for lib in m["libraries"]}
    rows = []
    for r in sheet:
        sid = r.get("sample_id")
        row = {k: r.get(k, "") for k in SAMPLE_COLUMNS}
        row["factors"] = {f: r.get(col, "") for f, col in columns.items() if r.get(col, "") != ""}
        row["deposit"] = {k: r[k] for k in ("tissue", "tissue_type", "life_stage", "dev_stage", "sex", "ploidy")
                          if r.get(k)}
        c = cur.get(sid)
        if c:
            keep = ["animal_id", *(f for f in ("tissue", "life_stage", "age", "sex", "ploidy", "sample_unit",
                                               "pool_size", "population", "treatment", "confidence"))]
            if full:
                keep += ["sources", "notes"]
            row["curated"] = {k: c.get(k, "") for k in keep}
        q = qc.get(sid)
        row["qc"] = {k: q.get(k) for k in QC_COLUMNS} if q else None
        row["atlas_unit"] = unit_of.get(sid)
        rows.append(row)
    page = _page(rows, limit, offset)
    return {"study_id": study_id,
            "metadata_source": "curated (sample_curation.tsv)" if cur else
                               "deposit only (no sample_curation.tsv yet); use the registry's study-level terms",
            "samples": page["rows"], "page": page["page"],
            "note": None if full or not cur else "curated sources and notes omitted; pass full=true",
            "sources": _sources(store, [study_id])}


# --------------------------------------------------------------------------- genes

def find_genes(store: Store, query: str, limit: int = 25, offset: int = 0) -> dict:
    """Genes by GeneID, symbol or product-description words."""
    ids = _match_genes(store, query)
    at = store.atlas()
    rows = [{**_gene_record(store, g), "in_atlas": g in at.gene_index} for g in ids]
    page = _page(rows, limit, offset)
    return {"query": query, "genes": page["rows"], "page": page["page"],
            "match": "description words" if ids and len(ids) > 1 else "id or symbol"}


def _group_means(pct_row: np.ndarray, units: list[dict], cols: list[int]) -> list[dict]:
    """Mean % per factor level, per condition and per tissue (when the study has several)."""
    groups: dict[tuple[str, str], list[float]] = {}
    for j in cols:
        m = units[j]
        keys = [(f, v) for f, v in m["factor_levels"].items()]
        if m.get("condition"):
            keys.append(("condition", m["condition"]))
        if m.get("tissue"):
            keys.append(("tissue", m["tissue"]))
        for k in keys:
            groups.setdefault(k, [])
            if not np.isnan(pct_row[j]):
                groups[k].append(float(pct_row[j]))
    by_dim: dict[str, int] = {}
    for f, _ in groups:
        by_dim[f] = by_dim.get(f, 0) + 1
    out = []
    for (f, v), vals in sorted(groups.items()):
        if by_dim[f] < 2:
            continue
        out.append({"group_by": f, "level": v, "n_measured": len(vals),
                    "mean_pct": nan_to_none(np.mean(vals)) if vals else None,
                    "sd": nan_to_none(np.std(vals, ddof=1)) if len(vals) > 1 else None})
    return out


def get_gene(store: Store, gene: str, include_groups: bool = True, include_differential: bool = True,
             max_genes: int = 3) -> dict:
    """Methylation of one gene (or the few matching a symbol / description) in every atlas study."""
    ids = _match_genes(store, gene)
    if not ids:
        return {"query": gene, "genes": [], "note": "no gene matches; try find_genes with fewer words"}
    at = store.atlas()
    cap = max(1, min(int(max_genes), MAX_GENES_DETAILED))
    shown = ids[:cap]
    genes = []
    for g in shown:
        rec = _gene_record(store, g)
        i = at.gene_index.get(g)
        if i is None:
            genes.append({**rec, "note": "not in the atlas gene tables"})
            continue
        studies = []
        for sid in at.study_ids:
            cols = at.columns[sid]
            entry = {"study_id": sid}
            for r in REGIONS:
                row = at.pct[r][i]
                vals = row[cols]
                ok = vals[~np.isnan(vals)]
                entry[r] = {"mean_pct": nan_to_none(_study_means(at.pct[r][i:i + 1], cols)[0]),
                            "units_measured": int(len(ok)), "units": len(cols),
                            "min": nan_to_none(ok.min()) if len(ok) else None,
                            "max": nan_to_none(ok.max()) if len(ok) else None}
                if include_groups and r == "gene_body":
                    entry["gene_body_groups"] = _group_means(row, at.units, cols)
            studies.append(entry)
        genes.append({**rec, "studies": studies})
    out = {"query": gene, "matched": len(ids), "shown": len(shown),
           "truncated": len(ids) > len(shown), "cap": MAX_GENES_DETAILED,
           "genes": genes,
           "method": f"per unit: % methylated calls over the region's CpGs, needing >= {GENE_MIN_CPGS} "
                     "covered CpGs; study mean_pct needs at least half the study's units measured",
           "sources": _sources(store, at.study_ids)}
    if include_differential:
        out["differential"] = gene_differential(store, ",".join(shown))["results"]
    return out


def _rank(store: Store, region, tissue, life_stage, chemistry, study_ids, min_studies):
    region = _region(region)
    cols, notes = _select_units(store, tissue=tissue, life_stage=life_stage, chemistry=chemistry,
                                study_ids=study_ids)
    if not cols:
        raise QueryError("no atlas units match these filters" + (f" ({'; '.join(notes)})" if notes else ""))
    sids, means, within = _study_matrix(store, region, cols)
    n = (~np.isnan(means)).sum(axis=1)
    need = len(sids) if min_studies is None else max(1, min(int(min_studies), len(sids)))
    return region, cols, notes, sids, means, within, n, need


def _gene_row(store, g, sids, means_row, extra):
    return {**_gene_record(store, g), **extra,
            "per_study": {sid: nan_to_none(v) for sid, v in zip(sids, means_row, strict=True)}}


def rank_genes(store: Store, region: str = "gene_body", order: str = "high", by: str = "mean",
               tissue=None, life_stage=None, chemistry=None, study_ids=None, min_studies=None,
               min_study_pct=None, max_study_pct=None, limit: int = 50, offset: int = 0) -> dict:
    """Genes ranked by methylation level across the selected studies."""
    if order not in ("high", "low"):
        raise QueryError("order must be 'high' or 'low'")
    if by not in ("mean", "min", "max"):
        raise QueryError("by must be 'mean', 'min' or 'max'")
    region, cols, notes, sids, means, _, n, need = _rank(store, region, tissue, life_stage, chemistry,
                                                         study_ids, min_studies)
    at = store.atlas()
    with _quiet():
        mean, lo, hi = np.nanmean(means, axis=1), np.nanmin(means, axis=1), np.nanmax(means, axis=1)
    keep = n >= need
    if min_study_pct is not None:
        keep &= lo >= float(min_study_pct)
    if max_study_pct is not None:
        keep &= hi <= float(max_study_pct)
    score = {"mean": mean, "min": lo, "max": hi}[by]
    idx = np.where(keep)[0]
    idx = idx[np.argsort(-score[idx] if order == "high" else score[idx], kind="mergesort")]
    page = _page(list(idx), limit, offset)
    rows = [_gene_row(store, at.genes[i], sids, means[i],
                      {"mean_pct": nan_to_none(mean[i]), "min_study_pct": nan_to_none(lo[i]),
                       "max_study_pct": nan_to_none(hi[i]), "n_studies": int(n[i])})
            for i in page["rows"]]
    return {"region": region, "order": order, "ranked_by": f"{by} of study means",
            "filters": _filters(tissue=tissue, life_stage=life_stage, chemistry=chemistry, study_ids=study_ids,
                                min_study_pct=min_study_pct, max_study_pct=max_study_pct),
            "studies": [{"study_id": s, "units": len(cols[s])} for s in sids],
            "min_studies": need, "genes_passing": int(keep.sum()),
            "genes": rows, "page": page["page"], "notes": notes,
            "method": f"study value = mean % over its selected units (>= {GENE_MIN_CPGS} CpGs per unit, "
                      "half the units measured); genes measured in fewer than min_studies are left out",
            "sources": _sources(store, sids)}


def variable_genes(store: Store, region: str = "gene_body", tissue=None, life_stage=None, chemistry=None,
                   study_ids=None, min_studies=None, min_mean_pct=None, limit: int = 50,
                   offset: int = 0) -> dict:
    """Genes whose methylation differs most between studies (SD of the study means)."""
    region, cols, notes, sids, means, within, n, need = _rank(store, region, tissue, life_stage, chemistry,
                                                              study_ids, min_studies)
    if len(sids) < 2:
        raise QueryError(f"variability across studies needs at least two studies; the filters leave {sids}")
    need = max(need, 2)
    at = store.atlas()
    with _quiet():
        sd = np.nanstd(means, axis=1, ddof=1)
        mean = np.nanmean(means, axis=1)
        rng = np.nanmax(means, axis=1) - np.nanmin(means, axis=1)
        wsd = np.nanmedian(within, axis=1)
    keep = (n >= need) & ~np.isnan(sd)
    if min_mean_pct is not None:
        keep &= mean >= float(min_mean_pct)
    idx = np.where(keep)[0]
    idx = idx[np.argsort(-sd[idx], kind="mergesort")]
    page = _page(list(idx), limit, offset)
    rows = [_gene_row(store, at.genes[i], sids, means[i],
                      {"between_study_sd": nan_to_none(sd[i]), "range": nan_to_none(rng[i]),
                       "mean_pct": nan_to_none(mean[i]), "median_within_study_sd": nan_to_none(wsd[i]),
                       "n_studies": int(n[i])})
            for i in page["rows"]]
    chems = sorted({m["chemistry"] for s in sids for m in (at.units[j] for j in cols[s])})
    caveats = ["Between-study differences mix biology (tissue, stage, population) with study effects; "
               "compare between_study_sd with median_within_study_sd."]
    if len(chems) > 1:
        caveats.append(f"Selected studies span chemistries {chems}; filter by chemistry to compare within one.")
    return {"region": region, "ranked_by": "SD of study means",
            "filters": _filters(tissue=tissue, life_stage=life_stage, chemistry=chemistry, study_ids=study_ids,
                                min_mean_pct=min_mean_pct),
            "studies": [{"study_id": s, "units": len(cols[s])} for s in sids], "min_studies": need,
            "genes_passing": int(keep.sum()), "genes": rows, "page": page["page"],
            "notes": notes, "caveats": caveats, "sources": _sources(store, sids)}


# --------------------------------------------------------------------------- differential

def _lambda_lookup(tests: list[dict]) -> dict:
    return {(t["job_id"], t["analysis"], t["region"], t["coef"]): t for t in tests}


def _num(v):
    try:
        x = float(v) if v not in ("", None, "NA") else None
    except ValueError:
        return None
    return int(x) if x is not None and x.is_integer() and "." not in str(v) and "e" not in str(v) else x


def list_contrasts(store: Store, study_id: str | None = None) -> dict:
    """Every tested coefficient of the shared DSS model, with counts and genomic inflation."""
    sids = [study_id] if study_id else sorted(store.summaries())
    if study_id:
        _check_study(store, study_id)
    out, missing = [], []
    for sid in sids:
        tests = store.differential_tests(sid)
        if not tests:
            missing.append(sid)
            continue
        summ = store.differential_summary(sid)
        rows = []
        for t in tests:
            lam = _num(t.get("lambda"))
            rows.append({"job_id": t["job_id"], "analysis": t["analysis"], "region": t["region"],
                         "coef": t["coef"], "units": _num(t.get("units")), "genes_tested": _num(t.get("genes_tested")),
                         "n_sig_fdr_0.05": _num(t.get("n_sig")), "up": _num(t.get("up")), "down": _num(t.get("down")),
                         "abs_diff_ge_10": _num(t.get("large")), "lambda": lam,
                         "inflated": lam is not None and lam > LAMBDA_WARN,
                         "outliers_removed": t.get("outliers") or None, "note": t.get("note") or None})
        jobs = [{"job_id": j.get("job_id"), "formula": j.get("formula"), "stratum": j.get("stratum")}
                for j in summ.get("jobs") or []]
        out.append({"study_id": sid, "generated_on": summ.get("generated_on"),
                    "gene_tables_from_run": summ.get("gene_tables_from_run"), "snp_mask": summ.get("snp_mask"),
                    "jobs": jobs, "skipped": summ.get("skipped") or [], "tests": rows})
    return {"studies": out, "without_differential": missing,
            "reading": f"lambda > {LAMBDA_WARN} means test statistics are inflated: significant counts "
                       "are not findings. 'without_pca_outliers' repeats each test without PCA outliers.",
            "sources": _sources(store, [s["study_id"] for s in out])}


def _coefs(header: list[str]) -> list[str]:
    return [h[2:] for h in header if h.startswith("p_")]


def get_differential(store: Store, study_id: str, job_id: str | None = None, region: str = "gene_body",
                     coef: str | None = None, max_fdr: float | None = 0.05, direction: str | None = None,
                     limit: int = 50, offset: int = 0) -> dict:
    """Genes from one study's contrast, most significant first."""
    _check_study(store, study_id)
    region = _region(region)
    tables = store.differential_tables(study_id)
    if not tables:
        return {"study_id": study_id, "available": False,
                "note": "no differential results for this study yet (`whale differential` writes "
                        "data/studies/<ID>/differential/)", "sources": _sources(store, [study_id])}
    jobs = sorted({j for j, _, _ in tables})
    job_id = job_id or ("main" if "main" in jobs else jobs[0])
    path = next((p for j, r, p in tables if j == job_id and r == region), None)
    if path is None:
        raise QueryError(f"no {region} table for job {job_id!r}; jobs: {', '.join(jobs)}")
    header, rows = store.differential_table(path)
    coefs = _coefs(header)
    coef = coef or coefs[0]
    if coef not in coefs:
        raise QueryError(f"coef {coef!r} not in job {job_id}; coefficients: {', '.join(coefs)}")
    h = {name: k for k, name in enumerate(header)}
    pc, fc, dc, sc = h[f"p_{coef}"], h[f"fdr_{coef}"], h.get(f"diff_{coef}"), h.get(f"stat_{coef}")
    mean_cols = [name for name in header if name.startswith("mean_")]
    picked = []
    for r in rows:
        p, fdr = _num(r[pc]), _num(r[fc])
        if p is None:
            continue
        if max_fdr is not None and (fdr is None or fdr > float(max_fdr)):
            continue
        stat = _num(r[sc]) if sc is not None else None
        diff = _num(r[dc]) if dc is not None else None
        sign = diff if diff is not None else stat
        if direction == "up" and not (sign is not None and sign > 0):
            continue
        if direction == "down" and not (sign is not None and sign < 0):
            continue
        picked.append((p, r, diff, stat, fdr))
    picked.sort(key=lambda x: x[0])
    page = _page(picked, limit, offset)
    genes = [{"gene_id": r[0], "symbol": r[h["symbol"]] if "symbol" in h else "",
              "description": r[h["description"]] if "description" in h else "",
              **{c: _num(r[h[c]]) for c in mean_cols},
              "diff_pp": diff, "stat": stat, "p": p, "fdr": fdr}
             for p, r, diff, stat, fdr in page["rows"]]
    lam = _lambda_lookup(store.differential_tests(study_id))
    t = lam.get((job_id, "all_units", region, coef), {})
    s = lam.get((job_id, "without_pca_outliers", region, coef))
    lam_v = _num(t.get("lambda"))
    return {"study_id": study_id, "job_id": job_id, "region": region, "coef": coef,
            "coefficients": coefs, "jobs": jobs, "analysis": "all_units",
            "filters": _filters(max_fdr=max_fdr, direction=direction),
            "genes_tested": len(rows), "genes_passing": len(picked),
            "lambda": lam_v, "inflated": lam_v is not None and lam_v > LAMBDA_WARN,
            "sensitivity_without_pca_outliers": ({"n_sig": _num(s.get("n_sig")), "lambda": _num(s.get("lambda")),
                                                  "outliers": s.get("outliers")} if s else None),
            "genes": genes, "page": page["page"],
            "columns": "diff_pp = difference in raw group means (percentage points, level vs reference); "
                       "stat = DSS Wald statistic; fdr = Benjamini-Hochberg within job, region and coefficient",
            "file": str(path.relative_to(store.repo)),
            "sources": _sources(store, [study_id])}


def gene_differential(store: Store, gene: str) -> dict:
    """One gene (or a comma-separated list of GeneIDs) in every study's shared-model tables."""
    ids = []
    for part in [p.strip() for p in gene.split(",") if p.strip()]:
        ids += _match_genes(store, part)[:MAX_GENES_DETAILED]
    ids = list(dict.fromkeys(ids))[:MAX_GENES_DETAILED]
    if not ids:
        return {"query": gene, "results": [], "note": "no gene matches"}
    results, without = [], []
    for sid in sorted(store.summaries()):
        tables = store.differential_tables(sid)
        if not tables:
            without.append(sid)
            continue
        lam = _lambda_lookup(store.differential_tests(sid))
        for job, region, path in tables:
            header, rows = store.differential_rows_for(path, set(ids))
            h = {name: k for k, name in enumerate(header)}
            for r in rows:
                for c in _coefs(header):
                    lam_v = _num(lam.get((job, "all_units", region, c), {}).get("lambda"))
                    results.append({"gene_id": r[0], "study_id": sid, "job_id": job, "region": region,
                                    "coef": c, "mean_pct": _num(r[h["mean_pct"]]) if "mean_pct" in h else None,
                                    "diff_pp": _num(r[h[f"diff_{c}"]]) if f"diff_{c}" in h else None,
                                    "p": _num(r[h[f"p_{c}"]]), "fdr": _num(r[h[f"fdr_{c}"]]),
                                    "lambda": lam_v, "inflated": lam_v is not None and lam_v > LAMBDA_WARN})
    results.sort(key=lambda x: (x["gene_id"], x["p"] if x["p"] is not None else 2))
    return {"query": gene, "genes": [_gene_record(store, g) for g in ids], "results": results,
            "studies_without_differential": without,
            "note": "a gene absent from a study's table was not tested there (fewer than the minimum "
                    "CpGs or calls in some unit)",
            "sources": _sources(store, sorted({r["study_id"] for r in results}))}


# --------------------------------------------------------------------------- replication

def get_replication(store: Store, study_id: str | None = None, verdict: str | None = None,
                    kind: str | None = None, text: str | None = None, limit: int = 50,
                    offset: int = 0) -> dict:
    """Replication verdicts for published claims, per study and per claim."""
    if verdict and verdict not in VERDICTS:
        raise QueryError(f"verdict must be one of {', '.join(VERDICTS)}")
    if study_id:
        _check_study(store, study_id)
    sids = [study_id] if study_id else sorted(set(store.registry()) | set(store.summaries()))
    studies, claims, missing = [], [], []
    for sid in sids:
        rep = store.replication(sid)
        if not rep:
            missing.append(sid)
            continue
        studies.append({"study_id": sid, "generated_on": rep.get("generated_on"), "counts": rep.get("counts"),
                        "based_on": rep.get("based_on"), "caveats": rep.get("caveats"),
                        "papers": [{"id": s.get("id"), "citation": s.get("citation"), "doi": s.get("doi")}
                                   for s in rep.get("sources") or []]})
        for c in rep.get("claims") or []:
            if verdict and c.get("verdict") != verdict:
                continue
            if kind and c.get("kind") != kind:
                continue
            if text and text.lower() not in (c.get("statement") or "").lower():
                continue
            claims.append({"study_id": sid, **c})
    page = _page(claims, limit, offset)
    return {"filters": _filters(study_id=study_id, verdict=verdict, kind=kind, text=text),
            "studies": studies, "claims": page["rows"], "page": page["page"],
            "studies_without_replication": missing,
            "reading": "Each verdict states the rule that produced it and the numbers it rests on "
                       "(docs/replication.md). 'partly' is also the cap when the paper's contrast is "
                       "narrower than the shared model or test statistics are inflated.",
            "sources": _sources(store, [s["study_id"] for s in studies])}
