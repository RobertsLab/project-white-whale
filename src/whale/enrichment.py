"""Record a finished enrichment-tier run (MeDIP, MBD-seq) on its study branch:
`whale record-enrichment`. The counterpart of `whale record-run` (whale.results)
for workflows/enrichment; docs/enrichment.md explains the method.

It checks the run is complete, then:

1. copies qc_summary.tsv and run_manifest.json into the study's data folder;
2. writes methylome_manifest.tsv (every per-sample file's size, SHA-256 and
   data-host URL: tiles, gene features, QC reports) and uploads those files
   to the data host under methylomes/<STUDY_ID>/<sample>/;
3. builds gene x sample matrices from the per-sample feature tables, under
   enrichment/ (not genes/, which holds base-resolution methylation):
   <class>.{count,cpm,obs_exp}.tsv.gz for gene_body and promoter, and
   cpgs.tsv.gz (each gene's CpGs, the same for every library);
4. writes study_summary.yaml with `tier: enrichment` and a statement of what
   the values can be compared with; no CpG methylation percentage is ever
   written for these studies;
5. moves status.yaml through aligned, called, qc_passed and summarized,
   commits and pushes, as record-run does.

`whale build-atlas` leaves every enrichment-tier study out of the pooled
base-resolution tables (`is_enrichment`).
"""
from __future__ import annotations

import csv
import gzip
import json
import shutil
import statistics
import subprocess
from datetime import date
from pathlib import Path

import yaml

from whale.checksums import sha256_file
from whale.claims import (
    ClaimError,
    _run,
    advance_status,
    current_branch,
    dump_status,
    machine_name,
    parse_branch,
)
from whale.paths import SUMMARY_NAME, status_path, studies_dir
from whale.results import (
    MANIFEST_COLUMNS,
    METHYLOME_MANIFEST,
    RECORDED_STAGES,
    DataHost,
    RecordError,
    _data_host,
    _float,
    _samplesheet,
    _stats,
    _write_gz,
    read_tsv,
)
from whale.vocab import load_vocab, stage_order

TIER = "enrichment"
FEATURE_SUFFIX = ".enrichment_features.tsv.gz"
TILE_SUFFIX = ".tiles_1kb.tsv.gz"
MATRIX_DIR = "enrichment"
CLASSES = ["gene_body", "promoter"]
METRICS = ["count", "cpm", "obs_exp"]
QC_KEYS = ["mapping_pct", "duplication_pct", "fragments", "relh", "goge", "cpgs_1x_pct", "saturation_est",
           "gene_body_obs_exp", "intergenic_obs_exp"]
COMPARABILITY = (
    "Enrichment tier (MeDIP/MBD): relative read density over methylated DNA, normalised for library size "
    "and CpG count. Comparable within this study (ranks, fold differences between groups) and, as ranks, "
    "with other enrichment studies; never pooled with or read as base-resolution CpG methylation "
    "percentages (PLAN.md §4.1)."
)


def chemistry_tier(repo: Path, chemistry: str | None) -> str | None:
    """The assay_chemistry tier of a term (base_resolution, enrichment), or None."""
    if not chemistry:
        return None
    try:
        return (load_vocab(Path(repo), "assay_chemistry").get(chemistry) or {}).get("tier")
    except FileNotFoundError:
        return None


def is_enrichment(repo: Path, summary: dict) -> bool:
    """A study summary from an enrichment-tier run: by its own tier, or by its chemistry's."""
    return summary.get("tier") == TIER or chemistry_tier(repo, summary.get("chemistry")) == TIER


# --------------------------------------------------------------------------- tables

def sample_files(results: Path, sample_id: str) -> list[Path]:
    files = [results / "features" / f"{sample_id}{FEATURE_SUFFIX}"]
    tiles = results / "tiles" / f"{sample_id}{TILE_SUFFIX}"
    if tiles.exists():
        files.append(tiles)
    qc = results / "qc"
    if qc.is_dir():
        files += sorted(p for p in qc.iterdir()
                        if p.is_file() and (p.name.startswith(f"{sample_id}_") or p.name.startswith(f"{sample_id}.")))
    return files


def gene_matrices(results: Path, samples: list[str], dest: Path) -> dict | None:
    """Gene x sample matrices under dest/enrichment/ and the summary's `genes` block.

    Rows are every gene of the reference bundle in its order; columns `samples`.
    obs_exp is empty where a gene's region has too few CpGs (enrichment_min_cpgs)."""
    if not samples:
        return None
    genes: list[str] | None = None
    cpgs: dict[str, list[str]] = {}
    data = {c: {m: [] for m in METRICS} for c in CLASSES}
    for sid in samples:
        path = results / "features" / f"{sid}{FEATURE_SUFFIX}"
        cols = {c: {m: [] for m in METRICS} for c in CLASSES}
        order: dict[str, list[str]] = {c: [] for c in CLASSES}
        n_cpgs: dict[str, list[str]] = {c: [] for c in CLASSES}
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                c = r["feature"]
                if c not in cols:
                    continue
                order[c].append(r["gene_id"])
                n_cpgs[c].append(r["cpgs"])
                for m in METRICS:
                    cols[c][m].append(r[m])
        for c in CLASSES:
            if genes is None:
                genes = order[c]
            if order[c] != genes:
                raise RecordError(f"{path.name}: its {c} genes differ from the other libraries' "
                                  "(made against another reference bundle?)")
            if c in cpgs and cpgs[c] != n_cpgs[c]:
                raise RecordError(f"{path.name}: its {c} CpG counts differ from the other libraries'")
            cpgs[c] = n_cpgs[c]
            for m in METRICS:
                data[c][m].append(cols[c][m])
    out = dest / MATRIX_DIR
    out.mkdir(exist_ok=True)
    _write_gz(out / "cpgs.tsv.gz", [["gene_id", *CLASSES],
                                    *([g, *(cpgs[c][i] for c in CLASSES)] for i, g in enumerate(genes))])
    block: dict = {"universe": len(genes), "values": "obs_exp: reads per CpG over the genome's reads per CpG"}
    for c in CLASSES:
        for m in METRICS:
            columns = data[c][m]
            _write_gz(out / f"{c}.{m}.tsv.gz",
                      [["gene_id", *samples], *([g, *(col[i] for col in columns)] for i, g in enumerate(genes))])
        # Genes with an obs_exp in every library: the study-wide mean per gene.
        means = []
        for i in range(len(genes)):
            vals = [_float(col[i]) for col in data[c]["obs_exp"]]
            if all(v is not None for v in vals):
                means.append(statistics.fmean(vals))
        block[c] = {"measured_in_all": len(means), "obs_exp": _stats(means),
                    "enriched_gt_1": sum(v > 1 for v in means)}
    return block


def summarize(study: dict, qc: list[dict], manifest: dict, *, part: str | None, run_name: str | None,
              today: str) -> dict:
    passed = [r for r in qc if r.get("qc_pass") == "true"]
    params = manifest.get("params", {})
    return {
        "study_id": study.get("study_id") or params.get("study_id"),
        "part": part,
        "tier": TIER,
        "comparability": COMPARABILITY,
        "generated_on": today,
        "nextflow_run": run_name or params.get("run_name"),
        "commit": params.get("commit"),
        "workflow": "enrichment",
        "genome_assembly": params.get("genome_assembly"),
        "chemistry": params.get("chemistry") or study.get("chemistry"),
        "enrichment_method": params.get("enrichment_method"),
        "window_size": params.get("window_size"),
        "libraries": len(qc),
        "passed": len(passed),
        "failed": [{"sample_id": r["sample_id"], "reasons": r.get("qc_reasons", "")}
                   for r in qc if r.get("qc_pass") != "true"],
        "qc": {k: _stats([_float(r.get(k)) for r in qc]) for k in QC_KEYS},
        # Per passing library: gene bodies over intergenic DNA, the first thing to look at in a
        # mosaic (invertebrate) methylome, where gene bodies carry most methylation.
        "gene_body_over_intergenic": _stats([
            _float(r.get("gene_body_obs_exp")) / _float(r.get("intergenic_obs_exp"))
            for r in passed if _float(r.get("gene_body_obs_exp")) and _float(r.get("intergenic_obs_exp"))]),
    }


# --------------------------------------------------------------------------- record

def record_enrichment_run(repo, results, *, run_name: str | None = None, upload: bool = True,
                          upload_target: str | None = None, commit: bool = True, push: bool = True,
                          remote: str = "origin", today: str | None = None, mark_ready: bool = True) -> dict:
    repo, results = Path(repo), Path(results)
    today = today or date.today().isoformat()
    branch = current_branch(repo)
    parsed = parse_branch(branch)
    if parsed is None:
        raise ClaimError(f"{branch!r} is not a study branch; run this inside the study's worktree.")
    _, study_id, part = parsed
    spath = status_path(repo, study_id, part)
    if not spath.exists():
        raise ClaimError(f"{spath.relative_to(repo)} is missing on this branch.")
    status = yaml.safe_load(spath.read_text())

    for name in ("run_manifest.json", "qc_summary.tsv"):
        if not (results / name).exists():
            raise RecordError(f"{results / name} is missing; did the run finish?")
    manifest = json.loads((results / "run_manifest.json").read_text())
    params = manifest.get("params", {})
    if params.get("tier") != TIER and chemistry_tier(repo, params.get("chemistry")) != TIER:
        raise RecordError(f"{results} is not an enrichment-tier run (chemistry {params.get('chemistry')!r}); "
                          "record it with `whale record-run`.")
    qc = read_tsv(results / "qc_summary.tsv")
    if not qc:
        raise RecordError(f"{results / 'qc_summary.tsv'} lists no libraries.")
    missing = [r["sample_id"] for r in qc if not (results / "features" / f"{r['sample_id']}{FEATURE_SUFFIX}").exists()]
    if missing:
        raise RecordError(f"no enrichment tables in {results / 'features'} for {len(missing)} librar(ies): "
                          f"{', '.join(missing[:5])}")
    run_name = run_name or params.get("run_name")

    dest = spath.parent
    shutil.copyfile(results / "qc_summary.tsv", dest / "qc_summary.tsv")
    shutil.copyfile(results / "run_manifest.json", dest / "run_manifest.json")

    host = _data_host(repo)
    url_base = (host.get("url_base") or "").rstrip("/")
    prefix = f"methylomes/{study_id}" + (f"/parts/{part}" if part else "")
    per_sample = {r["sample_id"]: sample_files(results, r["sample_id"]) for r in qc}
    rows = [{"sample_id": sid, "file": f.name, "bytes": f.stat().st_size, "sha256": sha256_file(f),
             "url": f"{url_base}/{prefix}/{sid}/{f.name}" if url_base else ""}
            for sid, files in per_sample.items() for f in files]
    with open(dest / METHYLOME_MANIFEST, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    uploaded, upload_error = False, None
    target = upload_target or host.get("upload_target")
    if upload and target:
        store = DataHost(target)
        try:
            for sid, files in per_sample.items():
                store.put(files, f"{prefix}/{sid}")
            uploaded = True
        except (RecordError, OSError) as exc:
            upload_error = str(exc)
    elif upload:
        upload_error = "data/data_host.yaml has no upload_target"

    study_yaml = studies_dir(repo) / f"{study_id}.yaml"
    study = yaml.safe_load(study_yaml.read_text()) if study_yaml.exists() else {}
    summary = summarize(study, qc, manifest, part=part, run_name=run_name, today=today)
    passed_ids = [r["sample_id"] for r in qc if r.get("qc_pass") == "true"]
    genes_block = gene_matrices(results, passed_ids, dest)
    if genes_block:
        summary["enrichment"] = genes_block
    summary["per_sample_files"] = {"files": len(rows), "bytes": sum(r["bytes"] for r in rows),
                                   "on_data_host": uploaded, "url": f"{url_base}/{prefix}/" if url_base else None}
    sheet = _samplesheet(repo, study_id, part)
    conditions = sorted({r.get("condition") for r in sheet if r.get("condition")})
    if conditions:
        summary["conditions"] = conditions
    (dest / SUMMARY_NAME).write_text(
        "# Written by `whale record-enrichment` from this study's enrichment run; read by the dashboard.\n"
        "# tier: enrichment -- never pooled with base-resolution methylation (whale build-atlas).\n"
        + yaml.safe_dump(summary, sort_keys=False, allow_unicode=True))

    n, n_pass = len(qc), summary["passed"]
    mapping = (summary["qc"]["mapping_pct"] or {}).get("median")
    notes = {
        "aligned": f"{n} enrichment libraries aligned (bowtie2); median mapping {mapping}%",
        "called": f"enrichment counted for {n} libraries (run {run_name})",
        "qc_passed": f"{n_pass} of {n} libraries pass the enrichment QC gate"
                     + (f"; failed: {', '.join(f['sample_id'] for f in summary['failed'])}" if summary["failed"] else ""),
        "summarized": f"enrichment summaries committed; {len(rows)} files on {host.get('host') or target}",
    }
    order = stage_order(repo)
    machine = machine_name()
    added = []
    for stage in RECORDED_STAGES:
        if order[stage] <= order.get(status["stage"], -1):
            continue
        if stage == "qc_passed" and n_pass == 0:
            break
        if stage == "summarized" and not uploaded:
            break
        status = advance_status(status, stage, repo=repo, machine=machine, today=today,
                                note=notes[stage], nextflow_run=run_name)
        added.append(stage)
    spath.write_text(dump_status(status))

    result = {"study_id": study_id, "part": part, "branch": branch, "stage": status["stage"], "added": added,
              "libraries": n, "passed": n_pass, "uploaded": uploaded, "upload_error": upload_error,
              "committed": False, "pushed": False, "pr_ready": False}
    if commit:
        files = [spath, dest / "qc_summary.tsv", dest / "run_manifest.json", dest / METHYLOME_MANIFEST,
                 dest / SUMMARY_NAME]
        if genes_block:
            files += sorted((dest / MATRIX_DIR).glob("*.tsv.gz"))
        _run(["git", "add", *[str(f.relative_to(repo)) for f in files]], repo)
        if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=repo).returncode != 0:
            label = f"{study_id}{f' part {part}' if part else ''}"
            msg = (f"{label}: {status['stage']} (record-enrichment {run_name})" if added
                   else f"{label}: enrichment summaries from {run_name}")
            _run(["git", "commit", "--quiet", "-m", msg], repo)
            result["committed"] = True
        if push:
            _run(["git", "push", "--quiet", remote, branch], repo)
            result["pushed"] = True
    if mark_ready and push and status["stage"] == "summarized" and shutil.which("gh"):
        proc = subprocess.run(["gh", "pr", "ready", branch], cwd=repo, capture_output=True, text=True)
        result["pr_ready"] = proc.returncode == 0
    return result


# --------------------------------------------------------------------------- dashboard

def _num(x, fmt="{:g}") -> str:
    return "—" if x is None else fmt.format(x)


def finding(entry: dict, study: dict | None, repo_url: str) -> dict:
    """The dashboard finding for an enrichment-tier study, from its study_summary.yaml."""
    s = entry["summary"]
    label = s["study_id"] + (f" ({s['part']})" if s.get("part") else "")
    qc = s.get("qc") or {}

    def med(k):
        return (qc.get(k) or {}).get("median")

    chem = {"medip": "MeDIP", "mbd": "MBD-seq"}.get(s.get("chemistry"), s.get("chemistry") or "enrichment")
    sentences = [
        f"Enrichment tier ({chem}): {s['passed']} of {s['libraries']} libraries passed the QC gate (median mapping "
        f"{_num(med('mapping_pct'))}%, duplication {_num(med('duplication_pct'))}%, CpG enrichment relH "
        f"{_num(med('relh'))}, saturation {_num(med('saturation_est'))})."
    ]
    ratio = (s.get("gene_body_over_intergenic") or {}).get("median")
    if ratio is not None:
        sentences.append(f"Gene bodies are {ratio:g}x as enriched per CpG as intergenic DNA (median over libraries).")
    e = s.get("enrichment") or {}
    gb = e.get("gene_body") or {}
    if gb.get("obs_exp"):
        sentences.append(f"Gene bodies: {gb['measured_in_all']:,} of {e['universe']:,} genes measured in every "
                         f"library; {gb['enriched_gt_1']:,} enriched above the genome's average CpG.")
    if s.get("failed"):
        sentences.append("Excluded: " + "; ".join(f"{x['sample_id']} ({x['reasons']})" for x in s["failed"]) + ".")
    sentences.append("Relative enrichment, not % methylation; kept out of the pooled atlas.")
    from whale.dashboard import excluded_summary

    stats = [
        {"label": "Passed QC", "value": f"{s['passed']} / {s['libraries']}"},
        {"label": "Mapping", "value": f"{_num(med('mapping_pct'))}%"},
        {"label": "Duplication", "value": f"{_num(med('duplication_pct'))}%"},
        {"label": "relH", "value": _num(med("relh"))},
        {"label": "Saturation", "value": _num(med("saturation_est"))},
    ]
    details = [f"{chem}: relative enrichment, not % methylation; kept out of the pooled atlas"]
    if ratio is not None:
        details.append(f"Gene bodies {ratio:g}× as enriched per CpG as intergenic DNA")
    if gb.get("obs_exp"):
        details.append(f"{gb['enriched_gt_1']:,} of {gb['measured_in_all']:,} genes enriched above the genome's average CpG")
    ref = entry["ref"] or "main"
    goals = (study or {}).get("goals") or []
    return {
        "date": str(s.get("generated_on") or ""),
        "title": (f"{label}: enrichment tier, {s['passed']} libraries"
                  + (f", gene bodies {ratio:g}x intergenic" if ratio is not None else "")
                  if s.get("passed") else f"{label}: no library passed the enrichment QC gate"),
        "summary": " ".join(sentences),
        "stats": stats,
        "details": details,
        "excluded": excluded_summary(s.get("failed")),
        "goal": ", ".join(str(g) for g in goals) or None,
        "studies": [s["study_id"]],
        "tier": TIER,
        "status": "provisional: study PR not merged" if entry["ref"] else "merged",
        "link": f"{repo_url}/blob/{ref}/{entry['path']}",
    }
