"""Replication reports: does a study's published finding hold when its data are
re-derived here? `whale replicate`.

Each study may have data/studies/<STUDY_ID>/claims.yaml: the paper's claims,
written from the paper without sight of our results (docs/replication.md and
docs/replication_claims_brief.md). Every claim has a `kind` and a `check`; this
module tests each check against the study's committed results and writes
data/studies/<STUDY_ID>/replication.md and replication.json.

Inputs, all committed:

* qc_summary.tsv, samplesheet.tsv           per-library CpG methylation, factor levels
* genes/{gene_body,promoter}.{meth,calls,cpgs}.tsv.gz   per-gene calls (SNP-masked)
* differential/ (summary.yaml, jobs.tsv, design.tsv, pca.tsv, <job>.<region>.tsv.gz)
  from `whale differential`; claims about differential methylation are
  `pending` until it has run
* data/reference/gene_info/<ASSEMBLY>.tsv.gz  symbols and descriptions

Each claim gets one verdict:

  replicated      the re-derived result agrees, by the rule for its kind
  partly          agrees in direction or presence but not in size or significance
  not_replicated  the re-derived result disagrees
  inconclusive    tested, but too little signal to decide
  not_testable    our results cannot address it (the claim says why)
  pending         testable once a missing input exists (e.g. the differential model)

The rules are deliberately simple and stated in each row, so a reader can
disagree with a verdict from the numbers beside it.
"""
from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from datetime import date
from math import comb
from pathlib import Path

import numpy as np
import yaml

from whale.registry import load_studies
from whale.results import _factor_column

VERDICTS = ["replicated", "partly", "not_replicated", "inconclusive", "not_testable", "pending"]
KINDS = {"global_level", "group_difference", "region_contrast", "methylated_fraction", "dm_count",
         "dm_direction", "dm_genes", "clustering", "not_testable"}
REGIONS = ("gene_body", "promoter")
MIN_CPGS = 5            # a gene's region counts in a unit when it has >= MIN_CPGS covered CpGs
LEVEL_TOL = (2.0, 5.0)  # global level: |difference| in points for replicated / partly
FRACTION_TOL = (5.0, 15.0)
MIN_SIG_FOR_DIRECTION = 5
COUNT_RATIO = 3.0       # a gene count within this factor of the paper's (either way) replicates it
LAMBDA_INFLATED = 1.3   # genomic inflation above which significant-gene counts are not trusted
PERMUTATIONS = 20000
SEED = 20261121


class ReplicationError(RuntimeError):
    pass


# --------------------------------------------------------------------------- inputs

def _tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


class Study:
    """The committed results of one study, loaded on demand."""

    def __init__(self, repo: Path, study_id: str):
        self.repo, self.id = Path(repo), study_id
        self.dir = self.repo / "data" / "studies" / study_id
        self.registry = load_studies(self.repo).get(study_id) or {}
        self.summary = yaml.safe_load((self.dir / "study_summary.yaml").read_text()) \
            if (self.dir / "study_summary.yaml").exists() else {}
        self.sheet = _tsv(self.dir / "samplesheet.tsv")
        self.qc = _tsv(self.dir / "qc_summary.tsv")
        self.diff_dir = self.dir / "differential"
        self.diff = yaml.safe_load((self.diff_dir / "summary.yaml").read_text()) \
            if (self.diff_dir / "summary.yaml").exists() else None
        self._genes = None
        self._info = None

    # factor levels per library --------------------------------------------------
    def factor_column(self, factor: str) -> str | None:
        columns = list(self.sheet[0]) if self.sheet else []
        for f in self.registry.get("factors") or []:
            if f["factor_id"] == factor:
                return _factor_column(f, columns, self.sheet)
        return factor if factor in columns else None

    def units(self) -> dict[str, list[str]]:
        """unit id -> QC-passing libraries (libraries of one BioSample are one unit)."""
        passed = [r["sample_id"] for r in self.qc if r.get("qc_pass") == "true"]
        shared = "shared_biosample" in (self.registry.get("quality_flags") or [])
        acc = {r["sample_id"]: r.get("sample_accession") or r["sample_id"] for r in self.sheet}
        out: dict[str, list[str]] = {}
        for lib in passed:
            out.setdefault(lib if shared else acc.get(lib, lib), []).append(lib)
        return out

    def unit_level(self, factor: str) -> dict[str, str]:
        col = self.factor_column(factor)
        if col is None:
            raise ReplicationError(f"factor {factor!r}: no sample-sheet column holds it")
        by_lib = {r["sample_id"]: str(r.get(col, "")) for r in self.sheet}
        return {u: by_lib[libs[0]] for u, libs in self.units().items()}

    def unit_global(self) -> dict[str, float]:
        """unit -> mean CpG methylation (%) of its libraries."""
        pct = {r["sample_id"]: _float(r.get("cpg_meth_pct")) for r in self.qc}
        out = {}
        for u, libs in self.units().items():
            vals = [pct[lib] for lib in libs if pct.get(lib) is not None]
            if vals:
                out[u] = float(np.mean(vals))
        return out

    # genes -----------------------------------------------------------------------
    def genes(self) -> dict:
        """{(region, metric): (gene ids, units, genes x units array)}, libraries summed per unit."""
        if self._genes is None:
            units = self.units()
            data = {}
            for region in REGIONS:
                for metric in ("meth", "calls", "cpgs"):
                    rows = _tsv(self.dir / "genes" / f"{region}.{metric}.tsv.gz")
                    if not rows:
                        raise ReplicationError(f"{self.id} has no gene tables")
                    genes = [r["gene_id"] for r in rows]
                    cols = []
                    for libs in units.values():
                        m = np.array([[_float(r.get(lib)) or 0.0 for lib in libs] for r in rows])
                        cols.append(m.max(axis=1) if metric == "cpgs" else m.sum(axis=1))
                    data[(region, metric)] = (genes, list(units), np.column_stack(cols))
            self._genes = data
        return self._genes

    def gene_info(self) -> list[dict]:
        if self._info is None:
            asm = self.summary.get("genome_assembly") or "GCF_963853765.1"
            self._info = _tsv(self.repo / "data" / "reference" / "gene_info" / f"{asm}.tsv.gz")
        return self._info

    # differential ----------------------------------------------------------------
    def job_for(self, factor: str, stratum: dict | None) -> tuple[dict | None, bool]:
        """The differential job testing `factor` within `stratum` (main model preferred), and
        whether it is only approximate: when no job has exactly that stratum, the job whose
        stratum is the largest part of it (it pools over the levels it leaves out)."""
        want = {str(k): str(v) for k, v in (stratum or {}).items()}
        exact, broader = [], []
        for j in _tsv(self.diff_dir / "jobs.tsv"):
            if factor not in j["tests"].split(","):
                continue
            got = dict(kv.split("=", 1) for kv in j["stratum"].split(",")) if j["stratum"] else {}
            if got == want:
                exact.append(j)
            elif got.items() <= want.items():
                broader.append((len(got), j))
        exact.sort(key=lambda j: (j["model_id"] != "main", j["job_id"]))
        if exact:
            return exact[0], False
        broader.sort(key=lambda x: (-x[0], x[1]["model_id"] != "main", x[1]["job_id"]))
        return (broader[0][1], True) if broader else (None, False)

    def job_summary(self, job_id: str) -> dict | None:
        for j in (self.diff or {}).get("jobs") or []:
            if j["job_id"] == job_id:
                return j
        return None


# --------------------------------------------------------------------------- statistics

def permutation_p(a: list[float], b: list[float], *, seed: int = SEED) -> float | None:
    """Two-sided permutation p for a difference in means (exact when small)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return None
    pooled = np.concatenate([a, b])
    obs = abs(a.mean() - b.mean())
    rng = np.random.default_rng(seed)
    n = len(a)
    hits = 0
    for _ in range(PERMUTATIONS):
        perm = rng.permutation(pooled)
        if abs(perm[:n].mean() - perm[n:].mean()) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (PERMUTATIONS + 1)


def r_squared(values: list[float], groups: list[str]) -> float | None:
    """Share of variance in `values` explained by `groups` (one-way), adjusted for the number of
    groups so that a factor with many levels and few samples is not credited by chance."""
    v = np.asarray(values, float)
    n, k = len(v), len(set(groups))
    if k < 2 or n <= k or v.var() == 0:
        return None
    total = ((v - v.mean()) ** 2).sum()
    within = sum(((v[[g == x for g in groups]] - v[[g == x for g in groups]].mean()) ** 2).sum()
                 for x in set(groups))
    return float(max(0.0, 1 - (within / (n - k)) / (total / (n - 1))))


def _group_units(st: Study, groups) -> tuple[set[str] | None, str]:
    """Units selected by a claim's `groups` ({factor: level} or 'all'), and a label."""
    if isinstance(groups, dict) and groups:
        sel = set(st.units())
        for factor, level in groups.items():
            lv = st.unit_level(factor)
            sel &= {u for u, x in lv.items() if x == str(level)}
        return sel, ", ".join(f"{f} = {v}" for f, v in groups.items())
    return None, "all samples"


def _gene_levels(st: Study, region: str, units: set[str] | None):
    """Per-gene mean methylation (%) over the selected units, for genes with >= MIN_CPGS
    covered CpGs in every one of them; (levels, number of units)."""
    genes, names, meth = st.genes()[(region, "meth")]
    _, _, calls = st.genes()[(region, "calls")]
    _, _, cpgs = st.genes()[(region, "cpgs")]
    cols = [i for i, u in enumerate(names) if units is None or u in units]
    m, c, k = meth[:, cols], calls[:, cols], cpgs[:, cols]
    ok = (k >= MIN_CPGS).all(axis=1) & (c > 0).all(axis=1)
    return (100 * m[ok] / c[ok]).mean(axis=1), len(cols)


def _r(x, d=2):
    return None if x is None else round(float(x), d)


# --------------------------------------------------------------------------- checks

def _out(verdict: str, ours: str, rule: str, **extra) -> dict:
    assert verdict in VERDICTS, verdict
    return {"verdict": verdict, "re_derived": ours, "rule": rule, **extra}


def check_global_level(st: Study, c: dict) -> dict:
    """Whole-library CpG methylation, or with `region` the mean per-gene level of that region."""
    sel, label = _group_units(st, c.get("groups", "all"))
    region = c.get("region")
    if region in REGIONS:
        levels, n = _gene_levels(st, region, sel)
        if not len(levels):
            return _out("not_testable", "no genes covered in every selected sample", "")
        ours = float(levels.mean())
        text = (f"{ours:.2f}% (mean {region} methylation of {len(levels)} genes with >= {MIN_CPGS} covered "
                f"CpGs in each of {n} samples, {label})")
    else:
        vals = {u: v for u, v in st.unit_global().items() if sel is None or u in sel}
        if not vals:
            return _out("not_testable", "no QC-passing samples in that group", "")
        ours = float(np.mean(list(vals.values())))
        text = f"{ours:.2f}% (mean over {len(vals)} samples, {label}; CpG calls pooled per library)"
    theirs = _float(c.get("level_pct"))
    rule = f"|difference| <= {LEVEL_TOL[0]} points: replicated; <= {LEVEL_TOL[1]}: partly"
    if theirs is None:
        return _out("inconclusive", text, "the paper gives no number", value=_r(ours))
    d = abs(ours - theirs)
    verdict = "replicated" if d <= LEVEL_TOL[0] else "partly" if d <= LEVEL_TOL[1] else "not_replicated"
    return _out(verdict, f"{text}; paper {theirs}%", rule, value=_r(ours), difference=_r(ours - theirs))


def check_group_difference(st: Study, c: dict) -> dict:
    factor = c["factor"]
    sel, where = _group_units(st, c.get("stratum") or "all")
    lv = st.unit_level(factor)
    vals = {u: v for u, v in st.unit_global().items() if sel is None or u in sel}
    if c.get("expect") == "no_difference" or c.get("higher") is None:
        levels = sorted({lv[u] for u in vals})
        if len(levels) != 2:
            return _out("not_testable", f"{factor} has {len(levels)} levels here; the test needs two", "")
        a = [v for u, v in vals.items() if lv[u] == levels[0]]
        b = [v for u, v in vals.items() if lv[u] == levels[1]]
        p = permutation_p(a, b)
        if p is None:
            return _out("not_testable", f"{len(a)} vs {len(b)} samples", "needs >= 2 samples per group")
        floor = 2 / comb(len(a) + len(b), len(a))
        text = (f"{factor} {levels[0]} {np.mean(a):.2f}% vs {levels[1]} {np.mean(b):.2f}% "
                f"(difference {np.mean(b) - np.mean(a):+.2f} points, permutation p = {p:.3g})")
        rule = ("the paper reports no difference: p >= 0.05 replicated; inconclusive when the group sizes "
                "cannot reach p < 0.05")
        if floor >= 0.05:
            return _out("inconclusive", f"{text}; {len(a)} vs {len(b)} samples cannot show a difference", rule,
                        p=_r(p, 4))
        return _out("replicated" if p >= 0.05 else "not_replicated", text, rule, p=_r(p, 4))
    hi, lo = str(c["higher"]), str(c["lower"])
    a = [v for u, v in vals.items() if lv.get(u) == hi]
    b = [v for u, v in vals.items() if lv.get(u) == lo]
    if len(a) < 2 or len(b) < 2:
        return _out("not_testable", f"{len(a)} vs {len(b)} samples", "needs >= 2 samples per group")
    diff, p = float(np.mean(a) - np.mean(b)), permutation_p(a, b)
    floor = 2 / comb(len(a) + len(b), len(a))      # the smallest two-sided p these group sizes allow
    separated = min(a) > max(b)
    text = (f"{factor} {hi} {np.mean(a):.2f}% vs {lo} {np.mean(b):.2f}% (difference {diff:+.2f} points, "
            f"n = {len(a)} vs {len(b)}, permutation p = {p:.3g}"
            + (f", the smallest possible with these group sizes is {floor:.2g}" if floor >= 0.05 else "")
            + (", groups completely separated" if separated else "")
            + ("" if sel is None else f", within {where}") + ")")
    rule = ("same direction and p < 0.05 (or, when group sizes cannot reach 0.05, complete separation): "
            "replicated; same direction: partly; opposite: not replicated")
    if diff > 0 and (p < 0.05 or (floor >= 0.05 and separated)):
        v = "replicated"
    elif diff > 0:
        v = "partly"
    else:
        v = "not_replicated"
    return _out(v, text, rule, difference=_r(diff), p=_r(p, 4))


def _region_pooled(st: Study, region: str) -> dict[str, float]:
    genes, units, meth = st.genes()[(region, "meth")]
    _, _, calls = st.genes()[(region, "calls")]
    tot = calls.sum(axis=0)
    return {u: float(100 * meth[:, i].sum() / tot[i]) for i, u in enumerate(units) if tot[i]}


def check_region_contrast(st: Study, c: dict) -> dict:
    hi, lo = c["higher"], c["lower"]
    level = {r: _region_pooled(st, r) for r in REGIONS}
    level["genome"] = st.unit_global()

    def mean(region):
        if region not in level:
            raise ReplicationError(f"region {region!r}: use gene_body, promoter or genome")
        return float(np.mean(list(level[region].values())))
    a, b = mean(hi), mean(lo)
    text = f"{hi} {a:.2f}% vs {lo} {b:.2f}% (CpG calls pooled over all genes, mean over samples)"
    return _out("replicated" if a > b else "not_replicated", text, "direction only", higher=_r(a), lower=_r(b))


def check_methylated_fraction(st: Study, c: dict) -> dict:
    if c.get("unit") == "cpgs":
        return _out("not_testable", "our tables are per gene, not per CpG", "", reason="needs_cpg_level")
    region = c.get("region") or "gene_body"
    if region not in REGIONS:
        return _out("not_testable", f"no gene-level table for {region}", "", reason="needs_other_regions")
    sel, label = _group_units(st, c.get("groups", "all"))
    pct, n = _gene_levels(st, region, sel)
    if not len(pct):
        return _out("not_testable", "no genes covered in every selected sample", "")
    theirs = _float(c.get("fraction_pct"))
    thr = _float(c.get("threshold_pct"))
    if thr is None:
        # A paper that calls a gene methylated by counting methylated CpGs has no level cut-off
        # to compare with; report the share at two cut-offs and leave the verdict open.
        f10, f50 = (float(100 * (pct >= x).mean()) for x in (10, 50))
        text = (f"{f10:.1f}% of {len(pct)} genes have mean {region} methylation >= 10%, {f50:.1f}% >= 50% "
                f"({n} samples, {label})" + (f"; paper {theirs}%" if theirs is not None else ""))
        return _out("inconclusive", text, "the paper's definition has no methylation cut-off to compare with",
                    ge10=_r(f10), ge50=_r(f50))
    ours = float(100 * (pct >= thr).mean())
    text = (f"{ours:.1f}% of {len(pct)} genes with >= {MIN_CPGS} covered CpGs in every sample have "
            f"mean {region} methylation >= {thr:g}% ({n} samples, {label})")
    rule = f"|difference| <= {FRACTION_TOL[0]} points: replicated; <= {FRACTION_TOL[1]}: partly"
    if theirs is None:
        return _out("inconclusive", text, "the paper gives no number", value=_r(ours))
    d = abs(ours - theirs)
    v = "replicated" if d <= FRACTION_TOL[0] else "partly" if d <= FRACTION_TOL[1] else "not_replicated"
    return _out(v, f"{text}; paper {theirs}%", rule, value=_r(ours), difference=_r(ours - theirs))


def _diff_context(st: Study, c: dict):
    """(job, coefficient, job summary) for a differential claim, or a verdict dict."""
    if st.diff is None:
        return _out("pending", "the shared differential model has not been run for this study", "")
    factor, level = c.get("factor"), c.get("level")
    if not factor or level is None:
        return _out("not_testable", "the claim names no registered factor and level", "", reason="other_data")
    job, approx = st.job_for(factor, c.get("stratum"))
    if job is None:
        stratum = c.get("stratum") or {}
        return _out("not_testable", f"no differential job tests {factor} within {stratum or 'all samples'}", "",
                    reason="other_data")
    coef = f"{factor}{level}"
    js = st.job_summary(job["job_id"]) or {}
    block = js.get("all_units") or {}
    tests = {r: ((block.get("regions") or {}).get(r) or {}).get("tests", {}).get(coef) for r in REGIONS}
    if not any(tests.values()):
        return _out("not_testable", f"job {job['job_id']} has no coefficient {coef}", "", reason="other_data")
    if approx:
        left = {k: v for k, v in (c.get("stratum") or {}).items() if f"{k}={v}" not in job["stratum"].split(",")}
        job = {**job, "approx": "approximate: our job pools over " + ", ".join(f"{k} (paper: {v} only)"
                                                                         for k, v in left.items())}
    return job, coef, tests


def _cap(res: dict) -> dict:
    if res["verdict"] == "replicated":
        res["verdict"] = "partly"
    return res


def _approx(job: dict, res: dict) -> dict:
    """A job that pools levels the paper kept apart tests a different contrast: at most partly."""
    if job.get("approx"):
        res["rule"] = f"{res['rule']}; {job['approx']}, so at most partly"
        res["approximate"] = True
        _cap(res)
    return res


def _inflated(tests: dict, res: dict) -> dict:
    """Inflated test statistics overstate significant-gene counts: at most partly."""
    lams = [t["lambda"] for t in tests.values() if t and t.get("lambda") is not None]
    if lams and max(lams) > LAMBDA_INFLATED:
        res["rule"] = (f"{res['rule']}; test statistics are inflated (lambda {max(lams)} > {LAMBDA_INFLATED}), "
                       "so the count is not trusted: at most partly")
        res["inflated"] = True
        _cap(res)
    return res


def check_dm_count(st: Study, c: dict) -> dict:
    ctx = _diff_context(st, c)
    if isinstance(ctx, dict):
        return ctx
    job, coef, tests = ctx
    parts, total = [], 0
    for r in REGIONS:
        t = tests[r]
        if t:
            total += t["fdr_lt_0.05"]
            parts.append(f"{r}: {t['fdr_lt_0.05']} genes at FDR < 0.05 (lambda {t['lambda']})")
    n = _float(c.get("n"))
    text = f"job {job['job_id']}, {coef}: " + "; ".join(parts)
    unit = str(c.get("unit") or "features")
    if unit == "DMGs":
        rule = (f"genes vs genes: within {COUNT_RATIO:g}x of the paper's count replicated; fewer or more "
                "but some: partly; none: not replicated")
    else:
        rule = (f"units differ (the paper counts {unit}, we count genes), so only presence is compared: "
                "both some partly, both none replicated, otherwise not replicated")
    if n is None:
        return _inflated(tests, _approx(job, _out("inconclusive", text, rule, value=total)))
    if n == 0 or total == 0:
        v = "replicated" if n == total == 0 else "not_replicated"
    elif unit == "DMGs":
        v = "replicated" if 1 / COUNT_RATIO <= total / n <= COUNT_RATIO else "partly"
    else:
        v = "partly"
    res = _out(v, f"{text}; paper {int(n)} {unit}", rule, value=total)
    return _inflated(tests, _approx(job, res))


def check_dm_direction(st: Study, c: dict) -> dict:
    ctx = _diff_context(st, c)
    if isinstance(ctx, dict):
        return ctx
    job, coef, tests = ctx
    up = sum((tests[r] or {}).get("up") or 0 for r in REGIONS)
    down = sum((tests[r] or {}).get("down") or 0 for r in REGIONS)
    text = f"job {job['job_id']}, {coef}: {up} higher and {down} lower in {c['level']} (FDR < 0.05, both regions)"
    rule = f"majority direction of >= {MIN_SIG_FOR_DIRECTION} significant genes"
    if up + down < MIN_SIG_FOR_DIRECTION:
        return _approx(job, _out("inconclusive", text, rule, up=up, down=down))
    ours = "hyper" if up > down else "hypo" if down > up else "balanced"
    return _approx(job, _out("replicated" if ours == c.get("direction") else "not_replicated", text, rule,
                             up=up, down=down))


def _norm(text: str) -> str:
    import re
    t = (text or "").lower()
    t = re.sub(r"\(.*?\)", " ", t)
    t = re.sub(r"\b(putative|probable|uncharacterized|isoform x\d+)\b", " ", t)
    t = re.sub(r"-like\b", " ", t)
    return " ".join(t.replace(",", " ").split())


def _match_genes(st: Study, name: str) -> list[dict]:
    """Annotation genes for a published name. Names may give alternatives ("bak1 / Bcl-2 ...");
    a symbol or a whole product description must match (ignoring case, "-like", "putative");
    a description containing the name is used only when nothing matches whole."""
    alts = [a.strip() for a in (name or "").split("/") if a.strip()]
    info = st.gene_info()
    whole, part = [], []
    for a in alts:
        q, qn = a.lower(), _norm(a)
        for g in info:
            if g["symbol"].lower() == q or (qn and _norm(g["description"]) == qn):
                whole.append(g)
            elif len(qn) >= 8 and qn in _norm(g["description"]):
                part.append(g)
    seen, out = set(), []
    for g in whole or part[:50]:
        if g["gene_id"] not in seen:
            seen.add(g["gene_id"])
            out.append(g)
    return out


def check_dm_genes(st: Study, c: dict) -> dict:
    ctx = _diff_context(st, c)
    if isinstance(ctx, dict):
        return ctx
    job, coef, _ = ctx
    tables = {r: {row["gene_id"]: row for row in _tsv(st.diff_dir / f"{job['job_id']}.{r}.tsv.gz")} for r in REGIONS}
    per_gene, counts = [], Counter()
    for g in c.get("genes") or []:
        matches = _match_genes(st, g.get("name"))
        regions = [g["region"]] if g.get("region") in REGIONS else list(REGIONS)
        best = None
        for m in matches:
            for r in regions:
                row = tables[r].get(m["gene_id"])
                p = _float((row or {}).get(f"p_{coef}"))
                if p is not None and (best is None or p < best[0]):
                    best = (p, r, m, row)
        if not matches:
            status, detail = "not_found", "no gene with this symbol or description in the annotation"
        elif best is None:
            status, detail = "not_tested", f"{len(matches)} matching gene(s), none passed the coverage filter"
        else:
            p, r, m, row = best
            fdr, dif = _float(row.get(f"fdr_{coef}")), _float(row.get(f"diff_{coef}"))
            k = len(matches) * len(regions)          # the best of k tests was taken: correct for it
            p_adj, fdr_adj = min(1.0, p * k), (None if fdr is None else min(1.0, fdr * k))
            want = g.get("direction")
            same = want is None or dif is None or (dif > 0) == (want == "hyper")
            if fdr_adj is not None and fdr_adj < 0.05 and same:
                status = "replicated"
            elif p_adj < 0.05 and same:
                status = "partly"
            else:
                status = "not_replicated"
            detail = (f"{m['symbol']} ({m['description']}), {r}: diff {dif:+.2f} points, p {p:.3g}, FDR {fdr:.3g}"
                      + (f"; best of {k} tests ({len(matches)} matching genes), so p and FDR are multiplied "
                         f"by {k} for the verdict" if k > 1 else "")
                      + ("" if same else f"; opposite to the paper's {want}"))
        counts[status] += 1
        per_gene.append({"name": g.get("name"), "status": status, "detail": detail})
    tested = counts["replicated"] + counts["partly"] + counts["not_replicated"]
    text = (f"job {job['job_id']}, {coef}: {counts['replicated']} replicated, {counts['partly']} partly, "
            f"{counts['not_replicated']} not, of {tested} tested ({counts['not_found']} not found in the "
            f"annotation, {counts['not_tested']} not tested)")
    rule = ("per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for "
            "the number of matching genes and regions); claim: most tested genes replicated, any "
            "replicated or partly: partly")
    if tested == 0:
        v = "not_testable"
    elif counts["replicated"] * 2 > tested:
        v = "replicated"
    elif counts["replicated"] or counts["partly"]:
        v = "partly"
    else:
        v = "not_replicated"
    return _approx(job, _out(v, text, rule, genes=per_gene))


def check_clustering(st: Study, c: dict) -> dict:
    if st.diff is None:
        return _out("pending", "the shared differential model (which computes the PCA) has not been run", "")
    comp_raw = str(c.get("component") or "PC1").upper()
    comps = [x for x in ("PC1", "PC2") if x in comp_raw]
    if not comps:
        return _out("not_testable", f"only a PCA's PC1 and PC2 are kept, not {c.get('component')}", "",
                    reason="other_data")
    absent = c.get("expect") in ("no_separation", "no_difference") or any(
        w in str(c.get("statement") or "").lower() for w in ("no separation", "do not separate", "not separate"))
    factor = c["factor"]
    pca = _tsv(st.diff_dir / "pca.tsv")
    design = {r["unit_id"]: r for r in _tsv(st.diff_dir / "design.tsv")}
    if factor not in next(iter(design.values()), {}):
        return _out("not_testable", f"{factor} is not in the differential design", "", reason="other_data")
    # the largest PCA in which the factor takes at least two levels
    by_job: dict[str, list[dict]] = {}
    for r in pca:
        by_job.setdefault(r["job_id"], []).append(r)
    usable = [(len(rs), j) for j, rs in by_job.items() if len({design[r["unit_id"]][factor] for r in rs}) >= 2]
    if not usable:
        return _out("not_testable", f"no PCA in which {factor} varies", "", reason="other_data")
    job = max(usable, key=lambda x: (x[0], x[1] == "main"))[1]
    rows = by_job[job]
    groups = [design[r["unit_id"]][factor] for r in rows]
    found = []
    for comp in comps:
        r2 = r_squared([float(r[comp.lower()]) for r in rows], groups)
        found.append((r2 or 0.0, comp, rows[0][f"{comp.lower()}_pct"]))
    r2, comp, pct = max(found)
    text = (f"{factor} explains {100 * r2:.0f}% (adjusted R^2) of {comp} ({pct}% of variance)"
            + (f", the larger of {' and '.join(comps)}" if len(comps) > 1 else "") + f"; PCA of job {job}")
    if absent:
        rule = "the paper reports no separation: adjusted R^2 < 0.2 replicated; < 0.5 partly"
        v = "replicated" if r2 < 0.2 else "partly" if r2 < 0.5 else "not_replicated"
    else:
        rule = "adjusted R^2 >= 0.5 replicated; >= 0.2 partly"
        v = "replicated" if r2 >= 0.5 else "partly" if r2 >= 0.2 else "not_replicated"
    return _out(v, text, rule, r2=_r(r2, 3))


CHECKS = {"global_level": check_global_level, "group_difference": check_group_difference,
          "region_contrast": check_region_contrast, "methylated_fraction": check_methylated_fraction,
          "dm_count": check_dm_count, "dm_direction": check_dm_direction, "dm_genes": check_dm_genes,
          "clustering": check_clustering}


# --------------------------------------------------------------------------- report

def evaluate(repo: Path, study_id: str) -> dict:
    path = Path(repo) / "data" / "studies" / study_id / "claims.yaml"
    if not path.exists():
        raise ReplicationError(f"{path} is missing")
    doc = yaml.safe_load(path.read_text()) or {}
    st = Study(repo, study_id)
    results = []
    for claim in doc.get("claims") or []:
        kind = claim.get("kind")
        if kind not in KINDS:
            raise ReplicationError(f"{claim.get('id')}: unknown kind {kind!r}")
        if kind == "not_testable":
            chk = claim.get("check") or {}
            res = _out("not_testable", chk.get("note") or "", "", reason=chk.get("reason"))
        else:
            try:
                res = CHECKS[kind](st, claim.get("check") or {})
            except (ReplicationError, KeyError) as exc:
                res = _out("not_testable", f"could not test: {exc}", "", reason="other_data")
        results.append({**{k: claim.get(k) for k in ("id", "source", "where", "statement", "kind")}, **res})
    s = st.summary
    return {
        "study_id": study_id,
        "generated_on": date.today().isoformat(),
        "sources": doc.get("sources") or [],
        "caveats": doc.get("caveats") or [],
        "based_on": {"nextflow_run": s.get("nextflow_run"), "genome_assembly": s.get("genome_assembly"),
                     "snp_mask": (s.get("snp_mask") or {}).get("mask_id"),
                     "differential": (st.diff or {}).get("generated_on"),
                     "samples": len(st.units())},
        "counts": {v: sum(r["verdict"] == v for r in results) for v in VERDICTS},
        "claims": results,
    }


def _cell(text) -> str:
    return str(text if text is not None else "").replace("|", "\\|").replace("\n", " ")


def render(rep: dict) -> str:
    b = rep["based_on"]
    lines = [f"# Replication: {rep['study_id']}", "",
             "Written by `whale replicate` (src/whale/replication.py) from `claims.yaml` and the study's "
             "committed results. Do not edit by hand; edit the claims or re-run.", "",
             "## Sources", ""]
    for s in rep["sources"]:
        lines.append(f"- **{s.get('id')}**: {s.get('citation')}" + (f" doi:{s['doi']}" if s.get("doi") else "")
                     + (f" — read: {s['read']}" if s.get("read") else ""))
    lines += ["", "## Re-derived from", "",
              f"- Pipeline run `{b['nextflow_run']}` on {b['genome_assembly']}, SNP mask `{b['snp_mask']}`, "
              f"{b['samples']} QC-passing samples",
              f"- Shared differential model: {'run ' + b['differential'] if b['differential'] else 'not run yet'}",
              "", *(["## Caveats", "", *[f"- {c}" for c in rep.get("caveats") or []], ""] if rep.get("caveats") else []),
              "## Summary", "",
              "| " + " | ".join(VERDICTS) + " |", "|" + "---|" * len(VERDICTS),
              "| " + " | ".join(str(rep["counts"][v]) for v in VERDICTS) + " |", "",
              "## Claims", "",
              "| id | where | claim (paraphrased) | verdict | re-derived | rule |", "|---|---|---|---|---|---|"]
    for r in rep["claims"]:
        lines.append("| " + " | ".join(_cell(x) for x in (
            r["id"], f"{r.get('source') or ''} {r.get('where') or ''}".strip(), r["statement"],
            f"**{r['verdict']}**", r["re_derived"], r["rule"])) + " |")
    genes = [r for r in rep["claims"] if r.get("genes")]
    if genes:
        lines += ["", "## Named genes", ""]
        for r in genes:
            lines += [f"**{r['id']}**", "", "| gene (as published) | status | re-derived |", "|---|---|---|"]
            lines += ["| " + " | ".join(_cell(x) for x in (g["name"], g["status"], g["detail"])) + " |"
                      for g in r["genes"]]
            lines.append("")
    lines += ["", "## Reading this", "",
              "Differential results here are per gene (gene body or promoter, summed calls, DSS "
              "beta-binomial model, FDR < 0.05), on GCF_963853765.1. Papers count CpGs or regions on other "
              "assemblies, so counts are compared by presence, and named genes are matched by symbol or "
              "product description. See docs/replication.md.", ""]
    return "\n".join(lines)


def write(repo: Path, study_id: str) -> dict:
    rep = evaluate(repo, study_id)
    d = Path(repo) / "data" / "studies" / study_id
    (d / "replication.json").write_text(json.dumps(rep, indent=1) + "\n")
    (d / "replication.md").write_text(render(rep))
    return rep


def studies_with_claims(repo: Path) -> list[str]:
    return sorted(p.parent.name for p in (Path(repo) / "data" / "studies").glob("*/claims.yaml"))
