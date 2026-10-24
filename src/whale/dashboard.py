"""Build site/data.json, the data behind the GitHub Pages dashboard (site/index.html).

Nothing on the page is typed in. It is assembled from:

* registry/backlog.tsv                      every known study
* registry/studies/*.yaml                   studies registered on main
* study branches' status.yaml               claims: machine, profile, stage (via `list_claims`)
* data/studies/*/status.yaml on main        stage of merged studies
* raw/<ACCESSION>/mirror_manifest.tsv       FASTQ mirrored on gannet (data_host.yaml url_base)
* data/studies/*/study_summary.yaml         findings, one per study (and part) whose run has been
                                            recorded (`whale record-run`), on main or its branch
* data/studies/*/genes/*.tsv.gz             gene x sample matrices: each study's gene section, with
                                            names from data/reference/gene_info/<ASSEMBLY>.tsv.gz

The Pages workflow (.github/workflows/pages.yml) runs `whale build-dashboard`
on every push to main and every few hours, so claims made on study branches
show up without a merge.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from whale.paths import SUMMARY_NAME
from whale.reference import gene_info_path, read_gene_info
from whale.registry import load_backlog, load_studies
from whale.vocab import load_vocab

DONE_STAGES = {"summarized", "merged"}
CLOSED_STAGES = {"excluded", "released"}


def _float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _summary_paths(base: Path, study_id: str) -> list[tuple[str | None, Path]]:
    d = base / "data" / "studies" / study_id
    out = [(None, d / SUMMARY_NAME)]
    if (d / "parts").is_dir():
        out += [(p.name, p / SUMMARY_NAME) for p in sorted((d / "parts").iterdir())]
    return [(part, path) for part, path in out if path.exists()]


def summaries_on_main(repo: Path) -> list[dict]:
    """Study summaries committed on main, each with ref = None."""
    out = []
    base = Path(repo) / "data" / "studies"
    for d in sorted(base.iterdir()) if base.is_dir() else []:
        for _part, path in _summary_paths(Path(repo), d.name):
            out.append({"summary": yaml.safe_load(path.read_text()) or {}, "ref": None,
                        "path": str(path.relative_to(repo))})
    return out


def summaries_on_branches(repo: Path, claims: list[dict] | None, remote: str = "origin") -> list[dict]:
    """Study summaries on the claim branches (already fetched by list_claims)."""
    out = []
    for c in claims or []:
        rel = f"data/studies/{c['study_id']}/" + (f"parts/{c['part']}/" if c.get("part") else "") + SUMMARY_NAME
        text = _git(Path(repo), "show", f"{remote}/{c['branch']}:{rel}")
        if text:
            out.append({"summary": yaml.safe_load(text) or {}, "ref": c["branch"], "path": rel})
    return out


def _pct(x) -> str:
    return "—" if x is None else f"{x:g}%"


def finding(entry: dict, study: dict | None, repo_url: str) -> dict:
    """One study's finding, stated from its study_summary.yaml and nothing else."""
    s = entry["summary"]
    label = s["study_id"] + (f" ({s['part']})" if s.get("part") else "")
    qc = s.get("qc") or {}
    meth = s.get("cpg_methylation_pct") or {}
    cpgs = (qc.get("cpgs_5x") or {}).get("median")
    sentences = [
        f"{s['passed']} of {s['libraries']} libraries passed the QC gate (median mapping "
        f"{_pct((qc.get('mapping_pct') or {}).get('median'))}, bisulfite conversion "
        f"{_pct((qc.get('conversion_pct') or {}).get('median'))}"
        + (f", {cpgs / 1e6:.1f} M CpGs at ≥5×" if cpgs else "") + ")."
    ]
    if meth:
        sentences.append(f"CpG methylation per library ranged {_pct(meth['min'])}–{_pct(meth['max'])} "
                         f"(median {_pct(meth['median'])}).")
    for f in s.get("by_factor") or []:
        parts = [f"{lv['level']} {_pct(lv.get('mean'))}" + (f" ± {lv['sd']:g}" if lv.get("sd") is not None else "")
                 + f" (n={lv.get('n', 0)})" for lv in f["levels"]]
        sentences.append(f"By {f['factor_id']}: " + ", ".join(parts) + ".")
    g = s.get("genes") or {}
    for cls, kind in (("gene_body", "Gene bodies"), ("promoter", "Promoters")):
        c = g.get(cls)
        if c and c.get("meth_pct"):
            sentences.append(f"{kind}: {c['covered_in_all']:,} of {g['universe']:,} genes have ≥{g['min_cpgs']} "
                             f"covered CpGs in every library; median methylation {_pct(c['meth_pct']['median'])}, "
                             f"{c['methylated_ge_50pct']:,} at ≥50%.")
    if s.get("failed"):
        sentences.append("Excluded: " + "; ".join(f"{x['sample_id']} ({x['reasons']})" for x in s["failed"]) + ".")
    if s.get("snp_mask"):
        sentences.append(f"Gene values leave out {s['snp_mask']['cpgs_masked']:,} CpGs on C/T SNPs "
                         f"(mask {s['snp_mask']['mask_id']}).")
    sentences.append("Descriptive; no differential test yet.")
    ref = entry["ref"] or "main"
    goals = (study or {}).get("goals") or []
    return {
        "date": str(s.get("generated_on") or ""),
        "title": f"{label}: median CpG methylation {_pct(meth.get('median'))} across {s['passed']} libraries"
                 if meth else f"{label}: no library passed the QC gate",
        "summary": " ".join(sentences),
        "goal": ", ".join(str(g) for g in goals) or None,
        "studies": [s["study_id"]],
        "status": "provisional: study PR not merged" if entry["ref"] else "merged",
        "link": f"{repo_url}/blob/{ref}/{entry['path']}",
    }


def build_findings(repo: Path, registered: dict, branch_entries: list[dict] | None, repo_url: str) -> list[dict]:
    seen = set()
    out = []
    # A branch's summary is newer than main's for the same study (a re-run before merging).
    for e in [*(branch_entries or []), *summaries_on_main(repo)]:
        key = (e["summary"].get("study_id"), e["summary"].get("part"))
        if not key[0] or key in seen:
            continue
        seen.add(key)
        out.append(finding(e, registered.get(key[0]), repo_url))
    return sorted(out, key=lambda f: (f["date"], f["title"]), reverse=True)


GENE_MIN_CPGS = 5        # a library counts towards a gene's mean only with >= 5 covered CpGs
HIST_BIN = 2             # % methylation per histogram bin
GENE_CLASSES = [("gene_body", "Gene body"), ("promoter", "Promoter")]


def _read(repo: Path, ref: str | None, rel: str, remote: str) -> bytes | None:
    """A file from the working tree (ref None: main as checked out) or from a remote branch."""
    if ref is None:
        path = Path(repo) / rel
        return path.read_bytes() if path.exists() else None
    proc = subprocess.run(["git", "show", f"{remote}/{ref}:{rel}"], cwd=repo, capture_output=True)
    return proc.stdout if proc.returncode == 0 else None


def _matrix(raw: bytes) -> tuple[list[str], list[str], list[list[str]]]:
    rows = list(csv.reader(io.StringIO(gzip.decompress(raw).decode()), delimiter="\t"))
    return rows[0][1:], [r[0] for r in rows[1:]], [r[1:] for r in rows[1:]]


def _mean(values: list[float]) -> float | None:
    return round(statistics.fmean(values), 1) if values else None


def gene_view(repo: Path, entry: dict, repo_url: str, remote: str = "origin") -> dict | None:
    """A study's gene section: per-gene means for the table, histograms, and links to the matrices.
    None if the study has no gene matrices (recorded before FEATURE_SUMMARY existed)."""
    s = entry["summary"]
    base = entry["path"].rsplit("/", 1)[0]
    raw = {}
    for cls, _ in GENE_CLASSES:
        for metric in ("meth_pct", "cpgs"):
            raw[(cls, metric)] = _read(repo, entry["ref"], f"{base}/genes/{cls}.{metric}.tsv.gz", remote)
    if any(v is None for v in raw.values()):
        return None
    sheet_raw = _read(repo, entry["ref"], f"data/studies/{s['study_id']}/samplesheet.tsv", remote)
    sheet = {r["sample_id"]: r for r in csv.DictReader(io.StringIO(sheet_raw.decode()), delimiter="\t")} \
        if sheet_raw else {}
    info_path = gene_info_path(repo, s.get("genome_assembly") or "")
    info = read_gene_info(info_path) if info_path.exists() else {}

    factors = [(f["factor_id"], f["column"], [lv["level"] for lv in f["levels"]]) for f in s.get("by_factor") or []]
    per_class, samples, genes = {}, None, None
    for cls, _ in GENE_CLASSES:
        smp, gns, pct = _matrix(raw[(cls, "meth_pct")])
        _, _, cpg = _matrix(raw[(cls, "cpgs")])
        samples, genes = samples or smp, genes or gns
        per_class[cls] = (pct, cpg)

    columns = ["gene_id", "symbol", "description", "biotype"]
    for cls, _ in GENE_CLASSES:
        columns += [cls, f"{cls}_n"]
    for fid, _, levels in factors:
        columns += [f"gene_body:{fid}={lv}" for lv in levels]
    rows, hists = [], {}
    for cls, _ in GENE_CLASSES:
        hists[cls] = {"bins": [0] * (100 // HIST_BIN), "genes": 0}
    for i, g in enumerate(genes):
        gi = info.get(g, {})
        row = [g, gi.get("symbol", ""), gi.get("description", ""), gi.get("biotype", "")]
        for cls, _ in GENE_CLASSES:
            pct, cpg = per_class[cls]
            ok = [float(pct[i][j]) for j in range(len(samples)) if pct[i][j] and int(cpg[i][j]) >= GENE_MIN_CPGS]
            m = _mean(ok)
            row += [m, len(ok)]
            if m is not None and len(ok) == len(samples):
                h = hists[cls]
                h["bins"][min(int(m // HIST_BIN), len(h["bins"]) - 1)] += 1
                h["genes"] += 1
        pct, cpg = per_class["gene_body"]
        for _, col, levels in factors:
            for lv in levels:
                row.append(_mean([float(pct[i][j]) for j, smp in enumerate(samples)
                                  if str(sheet.get(smp, {}).get(col, "")) == lv
                                  and pct[i][j] and int(cpg[i][j]) >= GENE_MIN_CPGS]))
        rows.append(row)

    ref = entry["ref"] or "main"
    return {
        "study_id": s["study_id"], "part": s.get("part"), "libraries": len(samples),
        "min_cpgs": GENE_MIN_CPGS, "bin": HIST_BIN, "genes": len(genes),
        "named": bool(info), "status": "provisional: study PR not merged" if entry["ref"] else "merged",
        "hist": [{"id": cls, "label": label, **hists[cls]} for cls, label in GENE_CLASSES],
        "factors": [{"id": fid, "levels": levels} for fid, _, levels in factors],
        "tables": [{"label": f"{label} {metric.replace('_', ' ')}",
                    "url": f"{repo_url}/raw/{ref}/{base}/genes/{cls}.{metric}.tsv.gz"}
                   for cls, label in GENE_CLASSES for metric in ("meth_pct", "cpgs", "meth", "calls")],
        "snp_mask": (s.get("snp_mask") or {}).get("mask_id"),
        "per_sample_url": (s.get("methylomes") or {}).get("url"),
        "columns": columns,
        "rows": rows,
    }


def build_gene_views(repo: Path, branch_entries: list[dict] | None, repo_url: str,
                     remote: str = "origin") -> list[dict]:
    seen, out = set(), []
    for e in [*(branch_entries or []), *summaries_on_main(repo)]:
        key = (e["summary"].get("study_id"), e["summary"].get("part"))
        if not key[0] or key in seen:
            continue
        seen.add(key)
        v = gene_view(repo, e, repo_url, remote)
        if v:
            out.append(v)
    return sorted(out, key=lambda v: v["study_id"])


def mirror_summary(ledger: dict[str, dict] | None) -> dict | None:
    if not ledger:
        return None
    runs = {r.get("run_accession") for r in ledger.values()}
    sources: dict[str, int] = {}
    for r in ledger.values():
        sources[r.get("source", "")] = sources.get(r.get("source", ""), 0) + 1
    return {"files": len(ledger), "runs": len(runs),
            "gb": round(sum(int(r.get("bytes") or 0) for r in ledger.values()) / 1e9, 2),
            "sources": dict(sorted(sources.items())),
            "last_mirrored": max((r.get("mirrored_on") or "" for r in ledger.values()), default="")}


def _main_status(repo: Path, study_id: str) -> dict | None:
    path = Path(repo) / "data" / "studies" / study_id / "status.yaml"
    return yaml.safe_load(path.read_text()) if path.exists() else None


def _git(repo: Path, *args: str) -> str:
    proc = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 else ""


def build(repo: Path, *, claims: list[dict] | None, mirrors: dict[str, dict | None] | None,
          repo_url: str, now: str | None = None, branch_summaries: list[dict] | None = None,
          remote: str = "origin") -> dict:
    """The dashboard's data. `claims` from list_claims (None if unavailable), `mirrors`
    maps accession -> ledger (None if unreadable), `branch_summaries` from
    summaries_on_branches."""
    repo = Path(repo)
    stages = load_vocab(repo, "pipeline_stage")
    stage_order = [s for s, v in sorted(stages.items(), key=lambda kv: kv[1]["order"]) if s not in CLOSED_STAGES]
    registered = load_studies(repo)
    by_accession = {(s.get("accessions") or {}).get("bioproject"): sid for sid, s in registered.items()}
    claims_by_acc: dict[str, list[dict]] = {}
    for c in claims or []:
        claims_by_acc.setdefault(c["accession"], []).append(c)

    rows = []
    for b in load_backlog(repo):
        acc = b["accession"]
        study_id = by_accession.get(acc)
        study = registered.get(study_id) if study_id else None
        held = claims_by_acc.get(acc, [])
        main_status = _main_status(repo, study_id) if study_id else None
        if held:
            current = held[0] if len(held) == 1 else max(held, key=lambda c: stages.get(c["stage"], {}).get("order", -1))
            stage, machine, profile, updated = current["stage"], current["machine"], current["profile"], current["updated"]
            state = "in progress"
        elif main_status:
            stage, machine, profile = main_status["stage"], main_status["machine"], main_status["profile"]
            updated = (main_status.get("history") or [{}])[-1].get("date", "")
            state = "done" if stage in DONE_STAGES else ("excluded" if stage == "excluded" else "registered")
        else:
            stage, machine, profile, updated = None, None, None, None
            state = "registered" if study_id else "not started"
        rows.append({
            "accession": acc,
            "study_id": study_id or (b["suggested_study_id"] if b["suggested_study_id"] != "(name on claim)" else None),
            "registered": bool(study_id),
            "chemistry": (study or {}).get("chemistry") or b["chemistry"],
            "runs": int(b["runs"]) if b["runs"].isdigit() else None,
            "fastq_gb": (study or {}).get("data_availability", {}).get("fastq_size_gb") or _float(b["fastq_gb"]),
            "paired_rnaseq": ((study or {}).get("paired_rnaseq") or {}).get("relationship") or b["paired_rnaseq"],
            "goals": (study or {}).get("goals") or [],
            "tissues": ((study or {}).get("samples") or {}).get("tissues") or [],
            "state": state,
            "stage": stage,
            "stage_index": stage_order.index(stage) if stage in stage_order else None,
            "machine": machine,
            "profile": profile or b["suggested_profile"],
            "updated": updated,
            "parts": [c["part"] for c in held if c["part"]],
            "branch_url": f"{repo_url}/tree/{held[0]['branch']}" if held else None,
            "mirror": mirror_summary((mirrors or {}).get(acc)),
            "notes": b["notes"],
        })

    def gb(key):
        return round(sum(r[key] or 0 for r in rows), 1)

    commit = _git(repo, "rev-parse", "--short", "HEAD")
    return {
        "generated_at": now or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "commit": commit,
        "repo_url": repo_url,
        "claims_available": claims is not None,
        "mirror_available": mirrors is not None,
        "stages": [{"id": s, "label": stages[s]["label"]} for s in stage_order],
        "totals": {
            "backlog": len(rows),
            "registered": sum(r["registered"] for r in rows),
            "in_progress": sum(r["state"] == "in progress" for r in rows),
            "done": sum(r["state"] == "done" for r in rows),
            "backlog_gb": gb("fastq_gb"),
            "mirrored_gb": round(sum((r["mirror"] or {}).get("gb", 0) for r in rows), 1),
            "mirrored_studies": sum(bool(r["mirror"]) for r in rows),
            "by_stage": {s: sum(r["stage"] == s for r in rows) for s in stage_order},
        },
        "studies": rows,
        "findings": build_findings(repo, registered, branch_summaries, repo_url),
        "genes": build_gene_views(repo, branch_summaries, repo_url, remote),
    }


def write(data: dict, out: Path) -> Path:
    """data.json, and each study's gene table as genes/<STUDY_ID>[.<PART>].json beside it
    (several MB each, so the page loads one only when asked)."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    data = dict(data)
    views = []
    for v in data.get("genes") or []:
        v = dict(v)
        rows = v.pop("rows", None)
        if rows is not None:
            name = v["study_id"] + (f".{v['part']}" if v.get("part") else "") + ".json"
            (out.parent / "genes").mkdir(exist_ok=True)
            (out.parent / "genes" / name).write_text(
                json.dumps({"columns": v["columns"], "rows": rows}, separators=(",", ":")))
            v["file"] = f"genes/{name}"
        views.append(v)
    if "genes" in data:
        data["genes"] = views
    out.write_text(json.dumps(data, indent=2) + "\n")
    return out
