"""Tests for the C/T SNP mask: workflows/snp_mask/bin/{pileup_sites.awk,snp_mask.py} and `whale record-mask`."""
from __future__ import annotations

import gzip
import importlib.util
import json
import shutil
import subprocess
import sys

import pytest
import yaml
from conftest import REPO

from whale.results import RecordError
from whale.snp_mask import record_mask

BIN = REPO / "workflows/snp_mask/bin"
_spec = importlib.util.spec_from_file_location("snp_mask_bin", BIN / "snp_mask.py")
sm = importlib.util.module_from_spec(_spec)
sys.modules["snp_mask_bin"] = sm
_spec.loader.exec_module(sm)

FMT = "GT:GL1:GQ:DP:SP:AC:AF1:CV:BT"
VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    # Two samples' definite T reads at a C: alt 3 over depth 9; Y reads are ignored.
    f"c1\t11\t.\tC\tT\t9\tPASS\tCX=CG\t{FMT}\t0/1:.:.:6:C4T2Y3:6:0.3:.:.\t0/0:.:.:4:C2T1Y5:3:0.3:.:.",
    # Only one alt read in total.
    f"c1\t20\t.\tG\tA\t5\tPASS\tCX=CG\t{FMT}\t0/0:.:.:6:G5A1R3:6:0.1:.:.\t0/0:.:.:4:G4:4:0:.:.",
    # No alternative allele at all.
    f"c1\t30\t.\tC\t.\t9\tPASS\tCX=CG\t{FMT}\t0/0:.:.:6:C6:6:0:.:.\t0/0:.:.:4:C4:4:0:.:.",
    # Reference N.
    f"c1\t40\t.\tN\tA\t9\tPASS\t.\t{FMT}\t0/1:.:.:6:A6:6:1:.:.\t0/0:.:.:4:A4:4:1:.:.",
    # Two alt bases, one sample missing SP.
    f"c2\t5\t.\tA\tG\t9\tPASS\t.\t{FMT}\t0/1:.:.:6:A3G2C1:6:0.5:.:.\t./.",
]) + "\n"


def test_pileup_sites_awk():
    out = subprocess.run(["awk", "-v", "min_alt=2", "-f", str(BIN / "pileup_sites.awk")],
                         input=VCF, capture_output=True, text=True, check=True).stdout
    assert [ln.split("\t") for ln in out.splitlines()] == [
        ["c1", "11", "C", "T", "9", "3", "2"],
        ["c2", "5", "A", "GC", "6", "3", "1"],
    ]


def _genome(path):
    #            0123456789012345678901234567890
    path.write_text(">c1\nAACGTTTACGAATTCGGGCGTT\n>c2\nACGT\n")
    return path


def _sites(path, rows):
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return path


def test_build_masks_cpgs_on_snps(tmp_path):
    genome = _genome(tmp_path / "g.fa")
    sites = _sites(tmp_path / "s.tsv", [
        ("c1", 3, "C", "T", 50, 10, 4),    # the C of CpG at 2
        ("c1", 10, "G", "A", 40, 5, 2),    # the G of CpG at 8
        ("c1", 13, "T", "C", 40, 5, 2),    # not in a CpG: a SNP, nothing masked
        ("c1", 16, "C", "T", 400, 3, 1),   # alt frac 0.0075 < 0.02: not a SNP
        ("c1", 19, "C", "T", 20, 2, 1),    # 2 alt reads < 3: not a SNP
        ("c2", 2, "C", "T", 9, 9, 3),      # CpG at 1 on c2
    ])
    m = sm.build(genome, [sites], min_alt=3, min_af=0.02, params={"mask_id": "M"}, out_dir=tmp_path)
    with gzip.open(tmp_path / "cpg_mask.bed.gz", "rt") as fh:
        bed = [ln.rstrip("\n").split("\t") for ln in fh]
    assert bed == [["c1", "2", "4", "3C>T", "10", "."], ["c1", "8", "10", "10G>A", "5", "."],
                   ["c2", "1", "3", "2C>T", "9", "."]]
    c = m["counts"]
    assert (c["snps"], c["cpgs"], c["cpgs_masked"]) == (4, 5, 3)
    assert m["params"] == {"mask_id": "M", "min_alt_reads": 3, "min_alt_frac": 0.02}
    with gzip.open(tmp_path / "snps.tsv.gz", "rt") as fh:
        assert len(fh.read().splitlines()) == 1 + 4
    # Same inputs, same bytes.
    first = (tmp_path / "cpg_mask.bed.gz").read_bytes()
    sm.build(genome, [sites], min_alt=3, min_af=0.02, params={"mask_id": "M"}, out_dir=tmp_path)
    assert (tmp_path / "cpg_mask.bed.gz").read_bytes() == first


def test_sites_on_unknown_sequences_are_refused(tmp_path):
    genome = _genome(tmp_path / "g.fa")
    sites = _sites(tmp_path / "s.tsv", [("c9", 3, "C", "T", 50, 10, 4)])
    with pytest.raises(ValueError, match="absent from the genome"):
        sm.build(genome, [sites], min_alt=3, min_af=0.02, params={}, out_dir=tmp_path)


def _mask_dir(tmp_path):
    d = tmp_path / "mask"
    d.mkdir()
    genome = _genome(tmp_path / "g.fa")
    sites = _sites(tmp_path / "s.tsv", [("c1", 3, "C", "T", 50, 10, 4)])
    sm.build(genome, [sites], min_alt=3, min_af=0.02, out_dir=d, params={
        "mask_id": "GCF_X.snpmask.v1", "genome_assembly": "GCF_X", "studies": ["S2026_A"],
        "samples": [{"sample_id": "a"}, {"sample_id": "b"}], "min_mapq": 40, "min_baseq": 20,
        "tool_versions": {"biscuit": "1.10.2"}, "commit": "abc"})
    return d


def test_record_mask(tmp_path):
    repo = tmp_path / "repo"
    (repo / "data").mkdir(parents=True)
    shutil.copy2(REPO / "data" / "data_host.yaml", repo / "data" / "data_host.yaml")
    d, host = _mask_dir(tmp_path), tmp_path / "gannet"
    res = record_mask(repo, d, upload_target=str(host), today="2026-10-24")
    up = host / "reference/GCF_X/snp_mask/GCF_X.snpmask.v1"
    assert sorted(p.name for p in up.iterdir()) == ["cpg_mask.bed.gz", "mask_manifest.json", "snps.tsv.gz"]
    reg = yaml.safe_load(res["registry"].read_text())["masks"]
    assert len(reg) == 1 and reg[0]["mask_id"] == "GCF_X.snpmask.v1" and reg[0]["libraries"] == 2
    assert reg[0]["counts"]["cpgs_masked"] == 1 and reg[0]["tool_versions"] == {"biscuit": "1.10.2"}
    assert reg[0]["files"][0]["url"].endswith("/reference/GCF_X/snp_mask/GCF_X.snpmask.v1/cpg_mask.bed.gz")
    # Recording again replaces the entry rather than adding one.
    record_mask(repo, d, upload=False, today="2026-10-25")
    reg = yaml.safe_load(res["registry"].read_text())["masks"]
    assert len(reg) == 1 and reg[0]["recorded_on"] == "2026-10-25"


def test_record_mask_refuses_changed_files(tmp_path):
    repo = tmp_path / "repo"
    (repo / "data").mkdir(parents=True)
    d = _mask_dir(tmp_path)
    with gzip.open(d / "cpg_mask.bed.gz", "wt") as fh:
        fh.write("c1\t0\t2\tx\t1\t.\n")
    with pytest.raises(RecordError, match="does not match"):
        record_mask(repo, d, upload=False)
    (d / "snps.tsv.gz").unlink()
    with pytest.raises(RecordError, match="lacks snps.tsv.gz"):
        record_mask(repo, d, upload=False)


def test_manifest_is_json(tmp_path):
    d = _mask_dir(tmp_path)
    assert json.loads((d / "mask_manifest.json").read_text())["counts"]["snps"] == 1
