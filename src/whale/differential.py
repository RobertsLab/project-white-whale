"""Re-derived differential methylation for one study (PLAN.md §6.3): every study's
registered contrasts run through one shared model, so "did it replicate" can be
asked of the same statistic everywhere.

Inputs, all committed on main or on the study branch:

* registry/studies/<ID>.yaml                 factors (levels, reference level) and an
                                             optional `differential` block
* data/studies/<ID>/study_summary.yaml       assembly and SNP mask of the gene tables
* data/studies/<ID>/samplesheet.tsv          the factor columns of every library
* data/studies/<ID>/genes/{class}.{metric}.tsv.gz   genes x libraries (meth, calls, cpgs)

`plan` turns these into jobs. A job is one model fitted to one set of biological
units: libraries of one BioSample (unless the study is flagged shared_biosample)
are summed into one unit first, as in the atlas, so technical replicates never
count as animals. Without a `differential` block the plan is:

* one additive model over every testable factor (a reference level, and at least
  two levels with at least two units each), testing each of them;
* when exactly two factors are testable and every combination has at least two
  units, the full model as well, testing their interaction.

A study whose animals were sampled more than once, or whose design is nested,
names its models in the block instead (see the schema and docs/differential.md).

`write_inputs` writes the plan (design.tsv, jobs.tsv) into
data/studies/<ID>/differential/, and `run` fits it there with
workflows/differential/dss.R in the pinned DSS container: DSS's beta-binomial
regression, one "site" per gene holding its summed calls, Wald tests per
coefficient, Benjamini-Hochberg within each test and region. `write_summary` then
turns the R output into summary.yaml.
"""
from __future__ import annotations

import csv
import itertools
import shutil
import subprocess
from collections import Counter
from datetime import date
from pathlib import Path

import yaml

from whale.registry import load_studies
from whale.results import _factor_column

DSS_IMAGE = "quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0"
SCRIPT = "workflows/differential/dss.R"
OUT_NAME = "differential"
REGIONS = ("gene_body", "promoter")
MIN_UNITS = 2           # per tested level (and per cell, for an interaction)
INDIVIDUAL = "individual"


class DifferentialError(Exception):
    pass


def _tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _libraries(study_dir: Path) -> list[str]:
    """Libraries in the gene tables, i.e. those that passed the QC gate."""
    import gzip
    path = study_dir / "genes" / "gene_body.meth.tsv.gz"
    if not path.exists():
        raise DifferentialError(f"{path} is missing: the study has no gene tables yet")
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")[1:]


def _terms(formula: str) -> list[str]:
    """The variables a formula like '~ a * b + c' uses, in order."""
    body = formula.split("~", 1)[-1]
    for ch in "+*:()":
        body = body.replace(ch, " ")
    return list(dict.fromkeys(body.split()))


def _interactions(formula: str) -> list[str]:
    """The interaction terms a formula implies, as 'a:b' (a * b implies a:b)."""
    out = []
    for part in formula.split("~", 1)[-1].split("+"):
        names = [n.strip() for n in part.replace("*", ":").split(":") if n.strip()]
        if len(names) > 1:
            out.append(":".join(names))
    return out


def design(repo: Path, study_id: str) -> dict:
    """Units x design columns for a study, with the factors' resolved columns.

    Returns {"study", "summary", "factors" [{id, column, levels, reference}],
    "units" [{unit_id, libraries, <factor id>: level, individual?}], "skipped" [reasons]}."""
    repo = Path(repo)
    study = load_studies(repo).get(study_id)
    if study is None:
        raise DifferentialError(f"registry/studies/{study_id}.yaml is missing")
    d = repo / "data" / "studies" / study_id
    if (d / "parts").is_dir():
        raise DifferentialError(f"{study_id} is processed in parts; differential runs on whole studies only")
    summary_path = d / "study_summary.yaml"
    if not summary_path.exists():
        raise DifferentialError(f"{summary_path.relative_to(repo)} is missing: record the run first")
    summary = yaml.safe_load(summary_path.read_text()) or {}
    if not (summary.get("snp_mask") or {}).get("mask_id"):
        raise DifferentialError(f"{study_id}'s gene tables are not SNP-masked; re-run with --snp_mask first")

    sheet = _tsv(d / "samplesheet.tsv")
    columns = list(sheet[0]) if sheet else []
    by_id = {r["sample_id"]: r for r in sheet}
    libs = _libraries(d)
    shared = "shared_biosample" in (study.get("quality_flags") or [])
    spec = study.get("differential") or {}
    ind_cols = spec.get("individual") or []
    for c in ind_cols:
        if c not in columns:
            raise DifferentialError(f"differential.individual: {c!r} is not a sample-sheet column")

    factors, skipped = [], []
    for f in study.get("factors") or []:
        col = _factor_column(f, columns, sheet)
        levels = [str(v) for v in f.get("levels") or []]
        if col is None or not levels or not {r.get(col, "") for r in sheet} & set(levels):
            skipped.append(f"factor {f['factor_id']}: no sample-sheet column holds its levels")
            continue
        ref = f.get("reference_level")
        factors.append({"id": f["factor_id"], "column": col, "levels": levels,
                        "reference": None if ref is None else str(ref)})

    units: dict[str, dict] = {}
    for lib in libs:
        r = by_id.get(lib)
        if r is None:
            raise DifferentialError(f"library {lib} is in the gene tables but not in the sample sheet")
        uid = lib if shared or not r.get("sample_accession") else r["sample_accession"]
        row = {"unit_id": uid, **{f["id"]: r.get(f["column"], "") for f in factors}}
        if ind_cols:
            row[INDIVIDUAL] = "_".join(r[c] for c in ind_cols)
        if uid in units:
            other = {k: v for k, v in units[uid].items() if k != "libraries"}
            if other != row:
                raise DifferentialError(f"libraries of unit {uid} disagree on the design: {other} vs {row}")
            units[uid]["libraries"].append(lib)
        else:
            units[uid] = {**row, "libraries": [lib]}
    return {"study": study, "summary": summary, "factors": factors,
            "units": list(units.values()), "skipped": skipped}


def _testable(f: dict, units: list[dict]) -> bool:
    n = Counter(u[f["id"]] for u in units)
    return f["reference"] is not None and sum(1 for lv in f["levels"] if n[lv] >= MIN_UNITS) >= 2 \
        and n[f["reference"]] >= MIN_UNITS


def _default_models(factors: list[dict], units: list[dict], skipped: list[str]) -> list[dict]:
    testable = []
    for f in factors:
        if _testable(f, units):
            testable.append(f["id"])
        elif f["reference"] is None:
            skipped.append(f"factor {f['id']}: no reference level, so not tested by default")
        else:
            skipped.append(f"factor {f['id']}: fewer than two levels with >= {MIN_UNITS} units")
    if not testable:
        return []
    models = [{"id": "main", "formula": "~ " + " + ".join(testable), "within": [], "test": testable}]
    if len(testable) == 2:
        a, b = (next(f for f in factors if f["id"] == t) for t in testable)
        n = Counter((u[a["id"]], u[b["id"]]) for u in units)
        if all(n[(x, y)] >= MIN_UNITS for x in a["levels"] for y in b["levels"]):
            models.append({"id": "interaction", "formula": f"~ {a['id']} * {b['id']}", "within": [],
                           "test": [f"{a['id']}:{b['id']}"]})
        else:
            skipped.append(f"interaction {a['id']}:{b['id']}: not every combination has >= {MIN_UNITS} units")
    return models


def plan(repo: Path, study_id: str) -> dict:
    """The jobs for one study, with the reasons anything was left out.

    Returns design() plus {"models", "jobs" [{job_id, model_id, stratum, formula,
    tests, units}]}."""
    out = design(repo, study_id)
    factors, units, skipped = out["factors"], out["units"], out["skipped"]
    by_id = {f["id"]: f for f in factors}
    spec = out["study"].get("differential") or {}
    if spec.get("models"):
        models = []
        for m in spec["models"]:
            used = _terms(m["formula"]) + list(m.get("within") or [])
            for t in used:
                if t == INDIVIDUAL and not spec.get("individual"):
                    raise DifferentialError(f"model {m['id']}: uses {INDIVIDUAL} but differential.individual is not set")
                if t != INDIVIDUAL and t not in by_id:
                    raise DifferentialError(f"model {m['id']}: {t!r} is not a resolved factor of the study")
            default_tests = [t for t in _terms(m["formula"]) if t != INDIVIDUAL and by_id[t]["reference"] is not None]
            tests = list(m.get("test") or default_tests + _interactions(m["formula"]))
            models.append({"id": m["id"], "formula": m["formula"], "within": list(m.get("within") or []),
                           "test": tests})
    else:
        if not factors:
            skipped.append("no factors registered")
        models = _default_models(factors, units, skipped)

    jobs = []
    for m in models:
        strata = [by_id[w]["levels"] for w in m["within"]]
        for combo in itertools.product(*strata) if strata else [()]:
            sel = [u for u in units if all(u[w] == lv for w, lv in zip(m["within"], combo, strict=True))]
            job_id = ".".join([m["id"], *combo])
            stratum = ",".join(f"{w}={lv}" for w, lv in zip(m["within"], combo, strict=True))
            main = [t for t in m["test"] if ":" not in t]
            short = [t for t in main if not _testable(by_id[t], sel)]
            if short:
                skipped.append(f"job {job_id}: {', '.join(short)} has fewer than two levels with >= {MIN_UNITS} units")
                continue
            if len(sel) <= len(_terms(m["formula"])) + 1:
                skipped.append(f"job {job_id}: {len(sel)} units are too few for {m['formula']}")
                continue
            jobs.append({"job_id": job_id, "model_id": m["id"], "stratum": stratum, "formula": m["formula"],
                         "tests": m["test"], "units": [u["unit_id"] for u in sel]})
    return {**out, "models": models, "jobs": jobs}


def out_dir(repo: Path, study_id: str) -> Path:
    return Path(repo) / "data" / "studies" / study_id / OUT_NAME


def write_inputs(p: dict, dest: Path) -> None:
    """design.tsv (one row per unit; levels in the order R should use, reference
    first, are in jobs.tsv) and jobs.tsv for dss.R. Clears earlier output."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    cols = ["unit_id", "libraries"] + [f["id"] for f in p["factors"]]
    if any(INDIVIDUAL in u for u in p["units"]):
        cols.append(INDIVIDUAL)
    with open(dest / "design.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for u in p["units"]:
            w.writerow([";".join(u[c]) if c == "libraries" else u[c] for c in cols])
    levels = {f["id"]: [f["reference"]] + [lv for lv in f["levels"] if lv != f["reference"]]
              if f["reference"] is not None else f["levels"] for f in p["factors"]}
    with open(dest / "jobs.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["job_id", "model_id", "stratum", "formula", "tests", "levels", "units"])
        for j in p["jobs"]:
            used = [t for t in _terms(j["formula"]) if t != INDIVIDUAL]
            w.writerow([j["job_id"], j["model_id"], j["stratum"], j["formula"], ",".join(j["tests"]),
                        ";".join(f"{t}=" + ",".join(levels[t]) for t in used), ",".join(j["units"])])


def command(repo: Path, study_id: str, engine: str) -> list[str]:
    """The command that runs dss.R for a study, from the repository root."""
    rel = f"data/studies/{study_id}/{OUT_NAME}"
    args = ["Rscript", SCRIPT, study_id, rel]
    if engine == "docker":
        return ["docker", "run", "--rm", "--platform", "linux/amd64", "-v", f"{Path(repo).resolve()}:/repo",
                "-w", "/repo", DSS_IMAGE, *args]
    if engine == "apptainer":
        return ["apptainer", "exec", "--bind", f"{Path(repo).resolve()}:/repo", "--pwd", "/repo",
                f"docker://{DSS_IMAGE}", *args]
    if engine == "local":
        return args
    raise DifferentialError(f"unknown engine {engine!r}: docker, apptainer or local")


def run(repo: Path, study_id: str, engine: str = "docker") -> None:
    res = subprocess.run(command(repo, study_id, engine), cwd=repo)
    if res.returncode != 0:
        raise DifferentialError(f"dss.R failed for {study_id} (exit {res.returncode})")


def write_summary(p: dict, dest: Path, *, today: str | None = None) -> dict:
    """summary.yaml from the plan and dss.R's tests.tsv / run_info.tsv."""
    tests = _tsv(dest / "tests.tsv") if (dest / "tests.tsv").exists() else []
    info = {r["key"]: r["value"] for r in _tsv(dest / "run_info.tsv")} if (dest / "run_info.tsv").exists() else {}
    jobs = []
    for j in p["jobs"]:
        rows = [r for r in tests if r["job_id"] == j["job_id"]]
        entry = {"job_id": j["job_id"], "model": j["model_id"], "formula": j["formula"],
                 "stratum": j["stratum"] or None, "units": len(j["units"])}
        failed = next((r["note"] for r in rows if r.get("note")), None)
        if failed:
            entry["note"] = failed
        for analysis in ("all_units", "without_pca_outliers"):
            sub = [r for r in rows if r["analysis"] == analysis and not r.get("note")]
            if not sub:
                continue
            block = {"units": int(sub[0]["units"]), "regions": {}}
            for r in sub:
                reg = block["regions"].setdefault(r["region"], {"genes_tested": int(r["genes_tested"]), "tests": {}})
                t = {"fdr_lt_0.05": int(r["n_sig"]), "lambda": float(r["lambda"])}
                if r["up"] != "":
                    t.update({"up": int(r["up"]), "down": int(r["down"]), "and_abs_diff_ge_10": int(r["large"])})
                reg["tests"][r["coef"]] = t
            if analysis == "without_pca_outliers":
                block["outliers"] = sub[0]["outliers"].split(",")
            entry[analysis] = block
        jobs.append(entry)
    s = p["summary"]
    out = {
        "study_id": s.get("study_id"),
        "generated_on": today or date.today().isoformat(),
        "genome_assembly": s.get("genome_assembly"),
        "snp_mask": (s.get("snp_mask") or {}).get("mask_id"),
        "gene_tables_from_run": s.get("nextflow_run"),
        "method": {"test": "DSS beta-binomial regression (DMLfit.multiFactor, Wald test per coefficient)",
                   "dss_version": info.get("dss_version"), "container": DSS_IMAGE,
                   "site": "one gene: its summed methylated and total calls in the region",
                   "filter": {"min_cpgs": int(info.get("min_cpgs", 5)), "min_calls": int(info.get("min_calls", 10)),
                              "every_unit": True},
                   "fdr": "Benjamini-Hochberg within each job, region and coefficient; significant at < 0.05",
                   "units": "libraries of one BioSample summed into one unit",
                   "sensitivity": "repeated without units more than 3 SD from the mean on gene-body PC1 or PC2"},
        "quality_flags": p["study"].get("quality_flags") or [],
        "units": len(p["units"]),
        "skipped": p["skipped"],
        "jobs": jobs,
    }
    (dest / "summary.yaml").write_text(
        "# Written by `whale differential` (src/whale/differential.py, workflows/differential/dss.R).\n"
        + yaml.safe_dump(out, sort_keys=False))
    return out
