"""Tests for `whale record-run` (src/whale/results.py), against a local remote and a local data host."""
from __future__ import annotations

import csv
import gzip
import json
import shutil

import pytest
import yaml
from conftest import REPO, git

from whale.checksums import sha256_file
from whale.claims import list_claims, set_status, start_study
from whale.dashboard import finding, summaries_on_branches
from whale.results import RecordError, record_run, summarize
from whale.scope import check_scope
from whale.validate import validate_status

# Fictional, as in test_claims.py.
ACC, SID = "PRJNA999903", "TESTER2099_RECORD_WGBS"
URL = "https://github.com/RobertsLab/project-white-whale"
QC_HEADER = ["sample_id", "condition", "mapping_pct", "duplication_pct", "cpg_meth_pct", "conversion_pct",
             "cpgs_5x", "mean_depth", "qc_pass", "qc_reasons"]
SAMPLES = [  # sample, ploidy, ph, mapping, cpg meth, pass
    ("S1", "diploid", "8.2", 30.0, 11.0, True),
    ("S2", "diploid", "7.7", 31.0, 12.0, True),
    ("S3", "triploid", "8.2", 29.0, 13.0, True),
    ("S4", "triploid", "7.7", 12.0, 14.0, False),
]


@pytest.fixture
def study(claim_env):
    """A claimed, registered study at fastq_staged, and a finished run's outputs."""
    clone, wt = claim_env["clone"], claim_env["tmp"] / "wt"
    (clone / "data").mkdir(exist_ok=True)
    shutil.copy2(REPO / "data" / "data_host.yaml", clone / "data" / "data_host.yaml")
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", "data host", cwd=clone)
    git("push", "--quiet", "origin", "main", cwd=clone)
    start_study(ACC, SID, repo=claim_env["clone"], profile="klone", worktree_dir=wt,
                registered_by="Tester", machine="klone", open_pr=False)
    yml = wt / f"registry/studies/{SID}.yaml"
    s = yaml.safe_load(yml.read_text())
    s["goals"] = [1]
    s["factors"] = [
        {"factor_id": "ploidy", "factor_original": "", "factor_standardized": "genotype_or_population",
         "levels": ["diploid", "triploid"], "reference_level": "diploid"},
        {"factor_id": "ph", "factor_original": "", "factor_standardized": "ocean_acidification",
         "levels": ["8.2", "7.7"], "reference_level": "8.2"},
    ]
    yml.write_text(yaml.safe_dump(s, sort_keys=False))
    sheet = wt / f"data/studies/{SID}/samplesheet.tsv"
    sheet.write_text("sample_id\tcondition\tploidy\twater_ph\n" + "".join(
        f"{sid}\t{p}_{ph.replace('.', '_')}\t{p}\t{ph}\n" for sid, p, ph, *_ in SAMPLES))
    git("add", "-A", cwd=wt)
    git("commit", "--quiet", "-m", "register", cwd=wt)
    set_status("fastq_staged", repo=wt, machine="klone")

    res = claim_env["tmp"] / "methylomes" / SID
    (res / "cpg").mkdir(parents=True)
    (res / "qc").mkdir()
    with open(res / "qc_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(QC_HEADER)
        for sid, p, ph, mapping, meth, ok in SAMPLES:
            w.writerow([sid, f"{p}_{ph}", mapping, 10.0, meth, 98.6, 2_000_000, 12.0,
                        "true" if ok else "false", "" if ok else "mapping 12.0% < 20.0%"])
            with gzip.open(res / "cpg" / f"{sid}.CpG_report.merged_CpG_evidence.cov.gz", "wt") as g:
                g.write("chr1\t0\t2\t50\t1\t1\n")
            (res / "qc" / f"{sid}_PE_report.txt").write_text("report")
            _features(res, sid, meth)
    (res / "run_manifest.json").write_text(json.dumps(
        {"params": {"study_id": SID, "run_name": "happy_turing", "commit": "abc123", "genome_assembly": "GCF_X",
                    "chemistry": "wgbs"}}))
    return {"wt": wt, "results": res, "host": claim_env["tmp"] / "gannet", **claim_env}


GENES = ["G1", "G2", "G3"]


def _features(res, sid, meth):
    """A FEATURE_SUMMARY table: G1 well covered, G2 thinly, G3 not at all."""
    (res / "features").mkdir(exist_ok=True)
    with gzip.open(res / "features" / f"{sid}.features.tsv.gz", "wt") as g:
        g.write("feature\tgene_id\tcpgs\tcpgs_5x\tmeth\tunmeth\tmeth_pct\tmean_meth_5x_pct\n")
        for cls in ("gene_body", "promoter", "exon", "intron"):
            for gene, n in zip(GENES, (10, 2, 0), strict=True):
                m = int(meth) * n if cls == "gene_body" else 0
                u = (100 - int(meth)) * n if cls == "gene_body" else 5 * n
                pct = f"{100 * m / (m + u):.3f}" if n else ""
                g.write(f"{cls}\t{gene}\t{n}\t{n}\t{m}\t{u}\t{pct}\t{pct}\n")


def _status(wt):
    return yaml.safe_load((wt / f"data/studies/{SID}/status.yaml").read_text())


def test_record_run_summarizes_uploads_and_pushes(study):
    wt, res, host = study["wt"], study["results"], study["host"]
    out = record_run(wt, res, upload_target=str(host), today="2026-10-24")
    assert out["added"] == ["aligned", "called", "qc_passed", "summarized"]
    assert (out["passed"], out["libraries"], out["uploaded"], out["pushed"]) == (3, 4, True, True)

    status = _status(wt)
    assert status["stage"] == "summarized" and status["nextflow_run"] == "happy_turing"
    assert [h["stage"] for h in status["history"]][-4:] == ["aligned", "called", "qc_passed", "summarized"]
    assert "3 of 4 libraries pass" in status["history"][-2]["note"] and "S4" in status["history"][-2]["note"]
    assert validate_status(status, wt) == []

    d = wt / f"data/studies/{SID}"
    with open(d / "methylome_manifest.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 12 and {r["sample_id"] for r in rows} == {"S1", "S2", "S3", "S4"}
    assert (host / "methylomes" / SID / "S2" / "S2.features.tsv.gz").exists()

    # Gene x sample matrices: every gene, the libraries that passed QC.
    with gzip.open(d / "genes" / "gene_body.meth_pct.tsv.gz", "rt") as fh:
        body = list(csv.reader(fh, delimiter="\t"))
    assert body[0] == ["gene_id", "S1", "S2", "S3"] and [r[0] for r in body[1:]] == GENES
    assert body[1][1:] == ["11.000", "12.000", "13.000"] and body[3][1:] == ["", "", ""]
    with gzip.open(d / "genes" / "promoter.cpgs.tsv.gz", "rt") as fh:
        assert [r[1] for r in csv.reader(fh, delimiter="\t")] == ["S1", "10", "2", "0"]
    assert sorted(p.name for p in (d / "genes").iterdir()) == sorted(
        f"{c}.{m}.tsv.gz" for c in ("gene_body", "promoter") for m in ("meth_pct", "cpgs", "calls"))
    cpg = next(r for r in rows if r["sample_id"] == "S1" and r["file"].endswith(".cov.gz"))
    assert cpg["url"] == f"https://gannet.fish.washington.edu/v1_web/project-white-whale/methylomes/{SID}/S1/{cpg['file']}"
    uploaded = host / "methylomes" / SID / "S1" / cpg["file"]
    assert sha256_file(uploaded) == cpg["sha256"]

    summary = yaml.safe_load((d / "study_summary.yaml").read_text())
    assert summary["cpg_methylation_pct"]["median"] == 12.0 and summary["failed"][0]["sample_id"] == "S4"
    ploidy, ph = summary["by_factor"]
    assert ploidy["column"] == "ploidy" and ph["column"] == "water_ph"
    assert [(lv["level"], lv["n"], lv.get("mean")) for lv in ploidy["levels"]] == [("diploid", 2, 11.5), ("triploid", 1, 13.0)]
    assert [(lv["level"], lv["n"]) for lv in ph["levels"]] == [("8.2", 2), ("7.7", 1)]
    genes = summary["genes"]
    assert genes["universe"] == 3 and genes["gene_body"]["covered_in_all"] == 1
    assert genes["gene_body"]["meth_pct"]["median"] == 12.0 and genes["gene_body"]["methylated_ge_50pct"] == 0

    # Committed and pushed, inside the study's own files.
    git("fetch", "--quiet", "origin", cwd=wt)
    assert git("status", "--porcelain", cwd=wt).strip() == ""
    changed = git("diff", "--name-only", "origin/main...HEAD", cwd=wt).split()
    branch = git("rev-parse", "--abbrev-ref", "HEAD", cwd=wt).strip()
    study_yaml = yaml.safe_load((wt / f"registry/studies/{SID}.yaml").read_text())
    assert check_scope(branch, changed, study_yaml, {}) == ([], [])

    # Another machine (the dashboard) reads the finding off the branch.
    other = study["tmp"] / "machine-b"
    git("clone", "--quiet", str(study["remote"]), str(other), cwd=study["tmp"])
    claims = list_claims(other)
    assert claims[0]["stage"] == "summarized"
    entries = summaries_on_branches(other, claims)
    f = finding(entries[0], study_yaml, URL)
    assert f["title"] == f"{SID}: median CpG methylation 12% across 3 libraries"
    assert "3 of 4 libraries passed the QC gate" in f["summary"] and "By ploidy: diploid 11.5%" in f["summary"]
    assert "Gene bodies: 1 of 3 genes have ≥5 covered CpGs in every library; median methylation 12%" in f["summary"]
    assert "Excluded: S4" in f["summary"] and f["status"].startswith("provisional") and f["goal"] == "1"
    assert f["link"] == f"{URL}/blob/{branch}/data/studies/{SID}/study_summary.yaml"

    # Running it again changes nothing.
    again = record_run(wt, res, upload_target=str(host), today="2026-10-24")
    assert again["added"] == [] and not again["committed"]


def test_failed_upload_stops_below_summarized_and_resumes(study):
    wt, res = study["wt"], study["results"]
    blocker = study["tmp"] / "not-a-dir"
    blocker.write_text("")
    out = record_run(wt, res, upload_target=str(blocker), today="2026-10-24")
    assert out["stage"] == "qc_passed" and not out["uploaded"] and out["upload_error"]
    assert yaml.safe_load((wt / f"data/studies/{SID}/study_summary.yaml").read_text())["methylomes"]["on_data_host"] is False
    out = record_run(wt, res, upload_target=str(study["host"]), today="2026-10-25")
    assert out["added"] == ["summarized"] and _status(wt)["history"][-1]["date"] == "2026-10-25"


def test_incomplete_run_is_refused(study):
    wt, res = study["wt"], study["results"]
    (res / "cpg" / "S2.CpG_report.merged_CpG_evidence.cov.gz").unlink()
    with pytest.raises(RecordError, match="no CpG calls"):
        record_run(wt, res, upload=False)
    assert _status(wt)["stage"] == "fastq_staged"


def test_no_library_passing_stops_at_called(study):
    wt, res = study["wt"], study["results"]
    text = (res / "qc_summary.tsv").read_text().replace("\ttrue\t", "\tfalse\t")
    (res / "qc_summary.tsv").write_text(text)
    out = record_run(wt, res, upload=False, push=False)
    assert out["added"] == ["aligned", "called"]


def test_feature_tables_missing_for_some_libraries(study):
    wt, res = study["wt"], study["results"]
    (res / "features" / "S2.features.tsv.gz").unlink()
    with pytest.raises(RecordError, match="feature tables missing"):
        record_run(wt, res, upload=False, push=False)


def test_runs_without_feature_tables_still_record(study):
    wt, res = study["wt"], study["results"]
    shutil.rmtree(res / "features")
    out = record_run(wt, res, upload=False, push=False)
    assert out["added"] == ["aligned", "called", "qc_passed"]
    assert "genes" not in yaml.safe_load((wt / f"data/studies/{SID}/study_summary.yaml").read_text())


def test_summarize_without_factors():
    qc = [{"sample_id": "A", "qc_pass": "true", "cpg_meth_pct": "10", "mapping_pct": "50"}]
    s = summarize({}, [], qc, {"params": {"study_id": "X2026_Y"}}, part=None, run_name=None, today="2026-01-01")
    assert s["study_id"] == "X2026_Y" and s["by_factor"] == [] and s["cpg_methylation_pct"]["sd"] is None
