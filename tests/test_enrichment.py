"""Tests for the enrichment tier: the counting script (workflows/enrichment/bin/
enrichment_counts.py), `whale record-enrichment` (src/whale/enrichment.py) and
the launch command. No network access; the data host is a local folder."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import shutil

import pytest
import yaml
from conftest import REPO, git

from whale.claims import config_name, nextflow_command, set_status, start_study
from whale.dashboard import finding
from whale.enrichment import is_enrichment, record_enrichment_run
from whale.results import RecordError, record_run
from whale.validate import validate_status

spec = importlib.util.spec_from_file_location("enrichment_counts",
                                              REPO / "workflows/enrichment/bin/enrichment_counts.py")
ec = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ec)

URL = "https://github.com/RobertsLab/project-white-whale"


# --------------------------------------------------------------------------- counting script

def test_fragment_geometry():
    # paired: one record per template, from the mate with TLEN > 0
    assert ec.fragment(["99", "c", "101", "50M", "180"], 200, 10_000) == (100, 280)
    assert ec.fragment(["147", "c", "231", "50M", "-180"], 200, 10_000) is None
    # single-end, forward: extended to the fragment length
    assert ec.fragment(["0", "c", "101", "50M", "0"], 200, 10_000) == (100, 300)
    # single-end, reverse: extended leftwards from the read's 3' end on the reference
    assert ec.fragment(["16", "c", "501", "40M2D10M", "0"], 200, 10_000) == (352, 552)
    # clipped to the sequence
    assert ec.fragment(["16", "c", "1", "50M", "0"], 200, 10_000) == (0, 50)
    assert ec.fragment(["0", "c", "9901", "50M", "0"], 200, 10_000) == (9900, 10_000)


def test_cpg_positions_and_obs_exp():
    assert list(ec.cpg_positions("ACGCGTTCG")) == [1, 3, 7]
    assert ec.obs_exp(10, 5, 100, 50, 5) == pytest.approx((10 / 100) / (5 / 50))
    assert ec.obs_exp(10, 4, 100, 50, 5) is None          # too few CpGs
    assert ec.obs_exp(0, 5, 0, 50, 5) is None              # empty library


@pytest.fixture
def tiny(tmp_path):
    """chrA: 3 kb; CpGs at 100..109 (5) in tile 0 and 2100..2119 (10) in tile 2; tile 1 has none."""
    seq = ["A"] * 3000
    for p in [*range(100, 110, 2), *range(2100, 2120, 2)]:
        seq[p], seq[p + 1] = "C", "G"
    (tmp_path / "genome.fa").write_text(">chrA desc\n" + "".join(seq).lower() + "\n")
    feats = tmp_path / "features"
    feats.mkdir()
    (feats / "genes.bed").write_text("chrA\t2000\t3000\tGB\tprotein_coding\t+\nchrA\t0\t500\tGA\tlncRNA\t+\n")
    # two promoter intervals for GB, overlapping: merged, so its CpGs count once
    (feats / "promoters.bed").write_text("chrA\t1900\t2110\tGB\t.\t+\nchrA\t2050\t2120\tGB\t.\t+\n"
                                         "chrA\t500\t600\tGA\t.\t+\n")
    (feats / "intergenic.bed").write_text("chrA\t500\t2000\tintergenic\t.\t.\n")
    # 8 paired fragments with midpoints in tile 2 (gene GB), 2 in tile 0 (gene GA), 2 single-end in tile 1
    rows = [f"99\tchrA\t{s + 1}\t50M\t100\n" for s in (2050,) * 8]
    rows += [f"163\tchrA\t{s + 1}\t50M\t100\n" for s in (50, 60)]
    rows += ["0\tchrA\t1201\t50M\t0\n", "16\tchrA\t1801\t50M\t0\n", "99\tchrZ\t1\t50M\t100\n"]
    with gzip.open(tmp_path / "frags.tsv.gz", "wt") as fh:
        fh.writelines(rows)
    return tmp_path


def _rows(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_count_tiles_features_and_qc(tiny):
    out = tiny / "out"
    qc = ec.count("X", tiny / "frags.tsv.gz", tiny / "genome.fa", tiny / "features", out,
                  window=1000, frag_len=200, min_cpgs=5)
    assert qc["fragments"] == 12 and qc["genome_cpgs"] == 15     # chrZ is not in the genome
    tiles = _rows(out / "X.tiles_1kb.tsv.gz")
    assert [(t["start"], t["cpgs"], t["count"]) for t in tiles] == [("0", "5", "2"), ("1000", "0", "2"),
                                                                    ("2000", "10", "8")]
    assert float(tiles[2]["cpm"]) == pytest.approx(8e6 / 12)
    assert float(tiles[2]["obs_exp"]) == pytest.approx((8 / 12) / (10 / 15))
    assert tiles[1]["obs_exp"] == ""                              # no CpGs
    feats = _rows(out / "X.enrichment_features.tsv.gz")
    assert [(f["feature"], f["gene_id"]) for f in feats] == [("gene_body", "GB"), ("gene_body", "GA"),
                                                             ("promoter", "GB"), ("promoter", "GA")]
    gb_promoter = feats[2]
    assert (gb_promoter["cpgs"], gb_promoter["count"]) == ("10", "8")   # 1900-2120 merged; mids at 2100
    assert feats[3]["obs_exp"] == "" and feats[3]["count"] == "0"
    # QC: two random halves; gene bodies vs intergenic; coverage of CpGs by whole fragments
    assert float(qc["gene_body_obs_exp"]) == pytest.approx((10 / 12) / (15 / 15), abs=1e-3)
    assert qc["intergenic_obs_exp"] == ""                         # no CpGs in 500-2000
    assert float(qc["relh"]) > 1                                  # fragments are CpG-richer than the genome
    assert qc["cpgs_1x_pct"] == "100.00" and qc["cpgs_5x_pct"] == "66.67"
    with open(out / "X.enrichment_qc.tsv") as fh:
        assert next(csv.DictReader(fh, delimiter="\t"))["fragments"] == "12"


def test_count_is_deterministic(tiny):
    a = ec.count("X", tiny / "frags.tsv.gz", tiny / "genome.fa", tiny / "features", tiny / "a")
    b = ec.count("X", tiny / "frags.tsv.gz", tiny / "genome.fa", tiny / "features", tiny / "b")
    assert a == b


def test_qc_table_gates(tmp_path):
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample_id\trun_accession\tfastq_1\tfastq_2\tcondition\nA\tR1\ta\t\tc\nB\tR2\tb\t\tt\n")
    r = tmp_path / "r"
    r.mkdir()
    (r / "A.bowtie2.log").write_text("1000 reads; of these:\n  ...\n91.50% overall alignment rate\n")
    (r / "B.bowtie2.log").write_text("1000 reads; of these:\n40.00% overall alignment rate\n")
    (r / "A.markdup.txt").write_text("READ: 900\nEXAMINED: 900\nDUPLICATE TOTAL: 90\nESTIMATED_LIBRARY_SIZE: 5000\n")
    for sid, frags, relh in (("A", 800, "1.800"), ("B", 10, "0.900")):
        (r / f"{sid}.enrichment_qc.tsv").write_text(
            "\t".join(ec.QC_COLUMNS) + "\n" + "\t".join(
                {"sample_id": sid, "fragments": str(frags), "relh": relh}.get(c, "") for c in ec.QC_COLUMNS) + "\n")
    rows = ec.qc_table(sorted(r.iterdir()), sheet, tmp_path / "qc.tsv", min_mapping=50.0, min_fragments=100,
                       min_relh=1.0)
    a, b = rows
    assert (a["mapping_pct"], a["duplication_pct"], a["reads"], a["qc_pass"]) == ("91.50", "10.00", 1000, "true")
    assert a["condition"] == "c" and a["relh"] == "1.800"
    assert b["qc_pass"] == "false"
    assert b["qc_reasons"] == "mapping 40.0% < 50.0%; fragments 10 < 100; CpG enrichment relH 0.9 < 1.0"


# --------------------------------------------------------------------------- record-enrichment

ACC, SID = "PRJNA999907", "TESTER2099_RECORD_MEDIP"
SAMPLES = [("E1", "ctl", True), ("E2", "trt", True), ("E3", "trt", False)]
GENES = ["G1", "G2", "G3"]


def _features(res, sid, scale):
    (res / "features").mkdir(exist_ok=True)
    with gzip.open(res / "features" / f"{sid}.enrichment_features.tsv.gz", "wt") as g:
        g.write("feature\tgene_id\tcpgs\tcount\tcpm\tobs_exp\n")
        for cls in ("gene_body", "promoter"):
            for gene, cpgs, n in zip(GENES, (50, 20, 2), (100, 10, 1), strict=True):
                count = n * scale if cls == "gene_body" else 1
                oe = f"{2.0 * scale * n / 100:.4f}" if cpgs >= 5 else ""
                g.write(f"{cls}\t{gene}\t{cpgs}\t{count}\t{count * 10:.3f}\t{oe}\n")


@pytest.fixture
def study(claim_env):
    clone, wt = claim_env["clone"], claim_env["tmp"] / "wt"
    (clone / "data").mkdir(exist_ok=True)
    shutil.copy2(REPO / "data" / "data_host.yaml", clone / "data" / "data_host.yaml")
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", "data host", cwd=clone)
    git("push", "--quiet", "origin", "main", cwd=clone)
    start_study(ACC, SID, repo=clone, profile="desktop", worktree_dir=wt, registered_by="Tester",
                machine="mac", open_pr=False)
    (wt / f"data/studies/{SID}/samplesheet.tsv").write_text(
        "sample_id\tcondition\n" + "".join(f"{s}\t{c}\n" for s, c, _ in SAMPLES))
    git("add", "-A", cwd=wt)
    git("commit", "--quiet", "-m", "register", cwd=wt)
    set_status("fastq_staged", repo=wt, machine="mac")

    res = claim_env["tmp"] / "methylomes" / SID
    (res / "qc").mkdir(parents=True)
    (res / "tiles").mkdir()
    cols = ["sample_id", "condition", "mapping_pct", "duplication_pct", "fragments", "relh", "saturation_est",
            "gene_body_obs_exp", "intergenic_obs_exp", "qc_pass", "qc_reasons"]
    with open(res / "qc_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for k, (sid, cond, ok) in enumerate(SAMPLES):
            w.writerow([sid, cond, 90.0 + k, 5.0, 2_000_000, 1.8, 0.95, 2.0 + k, 0.5,
                        "true" if ok else "false", "" if ok else "CpG enrichment relH 0.9 < 1.0"])
            _features(res, sid, k + 1)
            with gzip.open(res / "tiles" / f"{sid}.tiles_1kb.tsv.gz", "wt") as g:
                g.write("chrom\tstart\tend\tcpgs\tcount\tcpm\tobs_exp\nchr1\t0\t1000\t5\t3\t1.5\t1.2\n")
            (res / "qc" / f"{sid}.bowtie2.log").write_text("log")
    (res / "run_manifest.json").write_text(json.dumps({"params": {
        "study_id": SID, "run_name": "tidy_pasteur", "commit": "abc", "genome_assembly": "GCF_X",
        "chemistry": "medip", "tier": "enrichment", "window_size": 1000,
        "enrichment_method": {"name": "cpg_density_normalised_window_counts", "version": "1"}}}))
    return {"wt": wt, "results": res, "host": claim_env["tmp"] / "gannet", **claim_env}


def test_record_enrichment(study):
    wt, res, host = study["wt"], study["results"], study["host"]
    out = record_enrichment_run(wt, res, upload_target=str(host), today="2026-11-21")
    assert out["added"] == ["aligned", "called", "qc_passed", "summarized"]
    assert (out["passed"], out["libraries"], out["uploaded"], out["pushed"]) == (2, 3, True, True)
    status = yaml.safe_load((wt / f"data/studies/{SID}/status.yaml").read_text())
    assert status["stage"] == "summarized" and validate_status(status, wt) == []

    d = wt / f"data/studies/{SID}"
    assert not (d / "genes").exists()                        # base-resolution folder untouched
    assert sorted(p.name for p in (d / "enrichment").iterdir()) == sorted(
        ["cpgs.tsv.gz"] + [f"{c}.{m}.tsv.gz" for c in ("gene_body", "promoter") for m in ("count", "cpm", "obs_exp")])
    with gzip.open(d / "enrichment" / "gene_body.obs_exp.tsv.gz", "rt") as fh:
        m = list(csv.reader(fh, delimiter="\t"))
    assert m[0] == ["gene_id", "E1", "E2"] and [r[0] for r in m[1:]] == GENES
    assert m[1][1:] == ["2.0000", "4.0000"] and m[3][1:] == ["", ""]
    with gzip.open(d / "enrichment" / "cpgs.tsv.gz", "rt") as fh:
        assert list(csv.reader(fh, delimiter="\t"))[1] == ["G1", "50", "50"]

    s = yaml.safe_load((d / "study_summary.yaml").read_text())
    assert s["tier"] == "enrichment" and "never pooled" in s["comparability"]
    assert "cpg_methylation_pct" not in s and "genes" not in s
    assert s["enrichment_method"]["name"] == "cpg_density_normalised_window_counts"
    assert s["failed"] == [{"sample_id": "E3", "reasons": "CpG enrichment relH 0.9 < 1.0"}]
    assert s["gene_body_over_intergenic"]["median"] == 5.0          # 2/.5 and 3/.5: passing libraries only
    assert s["enrichment"]["gene_body"]["measured_in_all"] == 2 and s["enrichment"]["gene_body"]["enriched_gt_1"] == 1
    assert s["conditions"] == ["ctl", "trt"]
    assert (host / "methylomes" / SID / "E1" / "E1.tiles_1kb.tsv.gz").exists()
    with open(d / "methylome_manifest.tsv") as fh:
        assert len(list(csv.DictReader(fh, delimiter="\t"))) == 9   # features, tiles, log x 3

    f = finding({"summary": s, "ref": None, "path": f"data/studies/{SID}/study_summary.yaml"}, None, URL)
    assert f["tier"] == "enrichment" and f["title"].startswith(f"{SID}: enrichment tier, 2 libraries")
    assert "Enrichment tier (MeDIP)" in f["summary"] and "not % methylation" in f["summary"]
    assert "Excluded: E3" in f["summary"]

    again = record_enrichment_run(wt, res, upload_target=str(host), today="2026-11-21")
    assert again["added"] == [] and not again["committed"]


def test_record_run_and_record_enrichment_refuse_each_other(study):
    wt, res = study["wt"], study["results"]
    with pytest.raises(RecordError, match="record-enrichment"):
        record_run(wt, res, upload=False, commit=False)
    m = json.loads((res / "run_manifest.json").read_text())
    m["params"].update(chemistry="wgbs", tier=None)
    (res / "run_manifest.json").write_text(json.dumps(m))
    with pytest.raises(RecordError, match="not an enrichment-tier run"):
        record_enrichment_run(wt, res, upload=False, commit=False)


def test_missing_tables_are_refused(study):
    (study["results"] / "features" / "E2.enrichment_features.tsv.gz").unlink()
    with pytest.raises(RecordError, match="E2"):
        record_enrichment_run(study["wt"], study["results"], upload=False, commit=False)


# --------------------------------------------------------------------------- tier and launch

def test_is_enrichment_by_tier_or_chemistry(repo_root):
    assert is_enrichment(repo_root, {"tier": "enrichment"})
    assert is_enrichment(repo_root, {"chemistry": "mbd"})
    assert not is_enrichment(repo_root, {"chemistry": "wgbs"})
    assert not is_enrichment(repo_root, {})


def test_launch_command_picks_the_enrichment_workflow(tmp_path, repo_root):
    shutil.copytree(repo_root / "registry" / "controlled_vocabularies", tmp_path / "registry" / "controlled_vocabularies")
    (tmp_path / "config").mkdir()
    cfg = tmp_path / config_name("X2017_DEV_MEDIP")
    cfg.write_text("params {\n    chemistry    = 'medip'\n}\n")
    assert nextflow_command("X2017_DEV_MEDIP", "raven", repo=tmp_path)[2] == "workflows/enrichment"
    cfg.write_text("params {\n    chemistry    = 'wgbs'\n}\n")
    assert nextflow_command("X2017_DEV_MEDIP", "raven", repo=tmp_path)[2] == "workflows/methylation"
    assert nextflow_command("X2017_DEV_MEDIP", "raven")[2] == "workflows/methylation"
