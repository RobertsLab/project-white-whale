"""Tests for the RNA-seq sample sheet and `whale record-rnaseq` (src/whale/rnaseq.py)."""
from __future__ import annotations

import csv
import gzip
import json
import shutil

import pytest
import yaml
from click.testing import CliRunner
from conftest import REPO, git

import whale.rnaseq as rnaseq
from whale.checksums import sha256_file
from whale.claims import ClaimError
from whale.cli import main
from whale.results import RecordError
from whale.rnaseq import GENE_MATRICES, MANIFEST_NAME, SUMMARY_NAME, record_rnaseq

SID = "TESTER2099_RNASEQ_WGBS"
QC_HEADER = ["sample_id", "condition", "layout", "reads_in", "mapping_pct", "decoy_pct", "fragments_mapped",
             "strandedness", "genes_tpm_ge1", "qc_pass", "qc_reasons"]
SAMPLES = [("R1", "PE", 60.0, "reverse", True), ("R2", "SE", 55.0, "reverse", True),
           ("R3", "PE", 12.0, "unstranded", False)]


def _gz(path, text):
    with open(path, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(text.encode())


def make_results(root, study_id=SID, workflow="rnaseq"):
    """A finished workflows/rnaseq run's output folder."""
    res = root / "methylomes" / study_id / "rnaseq"
    (res / "matrices").mkdir(parents=True)
    (res / "qc").mkdir()
    ids = [s[0] for s in SAMPLES]
    with open(res / "qc_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(QC_HEADER)
        for sid, layout, mapping, strand, ok in SAMPLES:
            w.writerow([sid, "c", layout, 2_000_000, mapping, 3.0, 1_500_000, strand, 9000,
                        "true" if ok else "false", "" if ok else "mapping 12.0% < 30.0%"])
            s = res / "salmon" / sid
            (s / "aux_info").mkdir(parents=True)
            (s / "logs").mkdir()
            (s / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nT1\t100\t80\t1e6\t5\n")
            (s / "aux_info" / "meta_info.json").write_text("{}")
            (s / "lib_format_counts.json").write_text("{}")
            (s / "cmd_info.json").write_text("{}")
            (s / "logs" / "salmon_quant.log").write_text("log")
            (res / "qc" / f"{sid}.fastp.json").write_text("{}")
            (res / "qc" / f"{sid}_1_fastqc.zip").write_text("zip")
    header = "gene_id\t" + "\t".join(ids) + "\n"
    _gz(res / "matrices" / "gene_tpm.tsv.gz", header + "1\t5\t2\t0\n2\t0.5\t3\t0\n")
    _gz(res / "matrices" / "gene_counts.tsv.gz", header + "1\t50\t20\t0\n2\t5\t30\t0\n")
    _gz(res / "matrices" / "gene_length.tsv.gz", header + "1\t900\t900\t900\n2\t400\t410\t400\n")
    _gz(res / "matrices" / "transcript_tpm.tsv.gz", "transcript_id\tR1\tR2\tR3\nT1\t1\t1\t1\n")
    _gz(res / "matrices" / "transcript_counts.tsv.gz", "transcript_id\tR1\tR2\tR3\nT1\t1\t1\t1\n")
    (res / "matrices" / "matrices_report.json").write_text("{}")
    (res / "multiqc_report.html").write_text("<html></html>")
    (res / "run_manifest.json").write_text(json.dumps(
        {"params": {"study_id": study_id, "workflow": workflow, "run_name": "brave_salmon", "commit": "abc123",
                    "genome_assembly": "GCF_X"},
         "salmon_index": {"salmon": "1.10.3", "decoys": True, "kmer": 31}}))
    return res


@pytest.fixture
def env(claim_env):
    clone = claim_env["clone"]
    (clone / "data").mkdir(exist_ok=True)
    shutil.copy2(REPO / "data" / "data_host.yaml", clone / "data" / "data_host.yaml")
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", "data host", cwd=clone)
    git("push", "--quiet", "origin", "main", cwd=clone)
    return {**claim_env, "results": make_results(claim_env["tmp"]), "host": claim_env["tmp"] / "gannet"}


def test_record_rnaseq_commits_gene_matrices_and_uploads_the_rest(env):
    clone, res, host = env["clone"], env["results"], env["host"]
    git("switch", "--quiet", "-c", f"rnaseq/{SID}", cwd=clone)
    out = record_rnaseq(clone, res, upload_target=str(host), today="2026-11-21")
    assert (out["libraries"], out["passed"], out["uploaded"], out["committed"], out["pushed"]) == (3, 2, True, True, True)

    d = clone / "data" / "studies" / SID / "rnaseq"
    for name in ("qc_summary.tsv", "run_manifest.json", MANIFEST_NAME, SUMMARY_NAME, *GENE_MATRICES):
        assert (d / name).exists(), name
    assert not (d / "transcript_tpm.tsv.gz").exists()
    tracked = git("ls-tree", "-r", "--name-only", f"origin/rnaseq/{SID}", cwd=clone).split()
    assert f"data/studies/{SID}/rnaseq/gene_tpm.tsv.gz" in tracked
    # The methylation stage is not touched.
    assert not (clone / "data" / "studies" / SID / "status.yaml").exists()

    with open(d / MANIFEST_NAME) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["sample_id"] for r in rows} == {"R1", "R2", "R3", "_study"}
    r1 = [r["file"] for r in rows if r["sample_id"] == "R1"]
    assert r1 == ["quant.sf", "meta_info.json", "lib_format_counts.json", "cmd_info.json", "salmon_quant.log",
                  "R1.fastp.json", "R1_1_fastqc.zip"]
    study = {r["file"] for r in rows if r["sample_id"] == "_study"}
    assert study == {"transcript_tpm.tsv.gz", "transcript_counts.tsv.gz", "matrices_report.json", "multiqc_report.html"}
    q = next(r for r in rows if r["sample_id"] == "R2" and r["file"] == "quant.sf")
    assert q["url"] == f"https://gannet.fish.washington.edu/v1_web/project-white-whale/rnaseq/{SID}/R2/quant.sf"
    assert sha256_file(host / "rnaseq" / SID / "R2" / "quant.sf") == q["sha256"]
    assert (host / "rnaseq" / SID / "_study" / "multiqc_report.html").exists()

    summary = yaml.safe_load((d / SUMMARY_NAME).read_text())
    assert summary["passed"] == 2 and summary["failed"] == [{"sample_id": "R3", "reasons": "mapping 12.0% < 30.0%"}]
    assert summary["mixed_strandedness"] is True and summary["strandedness"] == {"reverse": 2, "unstranded": 1}
    assert summary["quantifier"] == {"salmon": "1.10.3", "decoys": True, "kmer": 31}
    # Gene 1 has TPM >= 1 in both passing libraries; gene 2 only in R2.
    assert summary["genes"] == {"universe": 2, "tpm_ge1_in_every_passing_library": 1, "matrices": GENE_MATRICES}
    assert summary["files"]["on_data_host"] is True and summary["qc"]["mapping_pct"]["median"] == 55.0

    # Running it again changes nothing, so commits nothing.
    again = record_rnaseq(clone, res, upload_target=str(host), today="2026-11-21")
    assert again["committed"] is False


def test_record_rnaseq_on_the_study_branch(env):
    clone = env["clone"]
    git("switch", "--quiet", "-c", f"study/PRJNA999904/{SID}", cwd=clone)
    out = record_rnaseq(clone, env["results"], upload=False, push=False, today="2026-11-21")
    assert out["committed"] and out["branch"] == f"study/PRJNA999904/{SID}" and out["upload_error"] is None


def test_record_rnaseq_refuses_other_branches_unless_not_committing(env):
    clone = env["clone"]
    with pytest.raises(ClaimError, match="neither"):
        record_rnaseq(clone, env["results"], upload=False)
    git("switch", "--quiet", "-c", "study/PRJNA999904/OTHER2099_STUDY_WGBS", cwd=clone)
    with pytest.raises(ClaimError, match="neither"):
        record_rnaseq(clone, env["results"], upload=False)
    out = record_rnaseq(clone, env["results"], upload=False, commit=False, push=False)
    assert out["committed"] is False and (clone / "data" / "studies" / SID / "rnaseq" / "gene_tpm.tsv.gz").exists()


def test_record_rnaseq_checks_the_run_is_complete_and_is_rnaseq(env, tmp_path):
    clone = env["clone"]
    git("switch", "--quiet", "-c", f"rnaseq/{SID}", cwd=clone)
    meth = make_results(tmp_path / "m", workflow=None)
    with pytest.raises(RecordError, match="not from workflows/rnaseq"):
        record_rnaseq(clone, meth, upload=False)
    (env["results"] / "salmon" / "R2" / "quant.sf").unlink()
    with pytest.raises(RecordError, match="no quant.sf .* R2"):
        record_rnaseq(clone, env["results"], upload=False)
    (env["results"] / "matrices" / "gene_counts.tsv.gz").unlink()
    with pytest.raises(RecordError, match="gene_counts.tsv.gz is missing"):
        record_rnaseq(clone, env["results"], upload=False)


def test_record_rnaseq_cli(env):
    clone = env["clone"]
    git("switch", "--quiet", "-c", f"rnaseq/{SID}", cwd=clone)
    r = CliRunner().invoke(main, ["record-rnaseq", "--results", str(env["results"]), "--no-push",
                                  "--upload-target", str(env["host"])], env={"WHALE_REPO": str(clone)})
    assert r.exit_code == 0, r.output
    assert "2 of 3 libraries pass" in r.output and "differ in strandedness" in r.output


# --- the sample sheet ---------------------------------------------------------

def _study(repo, block):
    s = yaml.safe_load((REPO / "registry/studies/_TEMPLATE.yaml").read_text())
    s["study_id"] = SID
    s["paired_rnaseq"] = block
    (repo / "registry" / "studies").mkdir(parents=True, exist_ok=True)
    (repo / "registry" / "studies" / f"{SID}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))


def test_fetch_rnaseq_samplesheet_reads_the_paired_bioproject(tmp_path, monkeypatch):
    _study(tmp_path, {"relationship": "same_animals", "bioproject": "PRJNA407831", "notes": None})
    seen = {}

    def fake(bioproject, **kw):
        seen.update(bioproject=bioproject, **kw)
        return {"bioproject": bioproject}

    monkeypatch.setattr(rnaseq, "build_samplesheet", fake)
    out = rnaseq.fetch_rnaseq_samplesheet(tmp_path, SID, condition_attributes=["treatment"])
    assert out == {"bioproject": "PRJNA407831", "relationship": "same_animals"}
    assert seen["out_dir"] == tmp_path / "data" / "studies" / SID / "rnaseq"
    assert seen["library_strategies"] == ["RNA-Seq"] and seen["study_id"] == SID


@pytest.mark.parametrize("block", [{"relationship": "none", "bioproject": None},
                                   {"relationship": "unknown", "bioproject": None}])
def test_fetch_rnaseq_samplesheet_needs_a_registered_bioproject(tmp_path, block):
    _study(tmp_path, block)
    with pytest.raises(RecordError, match="names no BioProject"):
        rnaseq.fetch_rnaseq_samplesheet(tmp_path, SID, condition_attributes=["x"])
    with pytest.raises(RecordError, match="register the study first"):
        rnaseq.fetch_rnaseq_samplesheet(tmp_path, "NOPE2099_X_WGBS", condition_attributes=["x"])


def test_every_registered_paired_bioproject_resolves():
    """Studies that name paired RNA-seq give fetch-rnaseq-samplesheet a BioProject."""
    for path in sorted((REPO / "registry" / "studies").glob("*.yaml")):
        study = yaml.safe_load(path.read_text())
        block = study.get("paired_rnaseq") or {}
        if block.get("bioproject"):
            bioproject, _ = rnaseq.paired_bioproject(study)
            assert bioproject.startswith(("PRJNA", "PRJEB", "PRJDB")), path.name
