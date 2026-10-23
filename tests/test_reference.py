"""Offline tests for `whale fetch-reference` (src/whale/reference.py).

A two-sequence genome and a small GTF stand in for NCBI's folder; the listing,
md5checksums.txt and downloads are served by fakes. No network access.
"""
from __future__ import annotations

import gzip
import hashlib
import json

import pytest
from click.testing import CliRunner

import whale.reference as ref
from whale.cli import main
from whale.reference import ReferenceError, build_features, check_bundle, fetch_bundle

ACC, FOLDER = "GCF_963853765.1", "GCF_963853765.1_xbMagGiga1.1"
PARENT = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/"

GENOME = ">chr1 test\n" + "ACGT" * 25 + "\n" + "ACGT" * 25 + "\n>chr2\n" + "CG" * 50 + "\n"   # 200 + 100 bp


def _gtf(source="NCBI RefSeq GCF_963853765.1-RS_2024_06"):
    def line(chrom, feat, s, e, strand, gene, tx="", extra=""):
        return (f'{chrom}\tGnomon\t{feat}\t{s}\t{e}\t.\t{strand}\t.\tgene_id "LOC{gene}"; '
                f'transcript_id "{tx}"; db_xref "GeneID:{gene}"; {extra}\n')
    return (
        f"#gtf-version 2.2\n#!genome-build xbMagGiga1.1\n#!annotation-source {source}\n"
        # + strand gene, two transcripts sharing exon 1
        + line("chr1", "gene", 21, 120, "+", 1, extra='gene_biotype "protein_coding";')
        + line("chr1", "transcript", 21, 120, "+", 1, "XM_1.1")
        + line("chr1", "exon", 21, 40, "+", 1, "XM_1.1")
        + line("chr1", "exon", 81, 120, "+", 1, "XM_1.1")
        + line("chr1", "transcript", 21, 100, "+", 1, "XM_2.1")
        + line("chr1", "exon", 21, 40, "+", 1, "XM_2.1")
        + line("chr1", "exon", 61, 100, "+", 1, "XM_2.1")
        # - strand gene near the end of chr1: promoter clipped at 200
        + line("chr1", "gene", 151, 190, "-", 2, extra='gene_biotype "lncRNA";')
        + line("chr1", "transcript", 151, 190, "-", 2, "XR_3.1")
        + line("chr1", "exon", 151, 190, "-", 2, "XR_3.1")
    )


def _gz(text: str) -> bytes:
    return gzip.compress(text.encode(), mtime=0)


class FakeNCBI:
    def __init__(self, gtf_source=None, corrupt=None):
        gtf = _gtf(gtf_source) if gtf_source else _gtf()
        self.files = {f"{FOLDER}_genomic.fna.gz": _gz(GENOME), f"{FOLDER}_genomic.gtf.gz": _gz(gtf),
                      f"{FOLDER}_rna.fna.gz": _gz(">XM_1.1\nACGT\n>XM_2.1\nACGT\n>XR_3.1\nACGT\n")}
        self.md5 = "".join(f"{hashlib.md5(b).hexdigest()}  ./{n}\n" for n, b in self.files.items())
        self.corrupt = corrupt or set()
        self.downloads = []

    def get_text(self, url):
        if url == PARENT:
            return f'<a href="{FOLDER}/">{FOLDER}/</a>'
        if url == f"{PARENT}{FOLDER}/md5checksums.txt":
            return self.md5
        raise AssertionError(url)

    def download(self, url, dest):
        name = url.rsplit("/", 1)[-1]
        self.downloads.append(name)
        dest.write_bytes(b"garbage" if name in self.corrupt else self.files[name])


# --- features ----------------------------------------------------------------

def test_build_features(tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text(_gtf())
    f = build_features(gtf, [("chr1", 200), ("chr2", 100)])
    assert f["genes"] == [("chr1", 20, 120, "1", "protein_coding", "+"), ("chr1", 150, 190, "2", "lncRNA", "-")]
    # exon 1 is shared by both transcripts and appears once
    assert [(r[1], r[2]) for r in f["exons"]] == [(20, 40), (60, 100), (80, 120), (150, 190)]
    assert [(r[1], r[2]) for r in f["introns"]] == [(40, 60), (40, 80)]
    # + strand: TSS 20 -> [-980, 220); - strand: TSS 190 -> [-10, 1190); both clipped to chr1 [0, 200)
    assert f["promoters"] == [("chr1", 0, 200, "1", ".", "+"), ("chr1", 0, 200, "2", ".", "-")]
    long = build_features(gtf, [("chr1", 5000), ("chr2", 100)])["promoters"]
    assert long == [("chr1", 0, 220, "1", ".", "+"), ("chr1", 0, 1190, "2", ".", "-")]
    assert [(r[0], r[1], r[2]) for r in f["intergenic"]] == [("chr1", 0, 20), ("chr1", 120, 150),
                                                              ("chr1", 190, 200), ("chr2", 0, 100)]


def test_gtf_sequence_missing_from_genome(tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text(_gtf())
    with pytest.raises(ReferenceError, match="absent from the genome"):
        build_features(gtf, [("chr2", 100)])


# --- fetch / check -----------------------------------------------------------

def test_fetch_bundle_end_to_end(repo_root, tmp_path):
    ncbi = FakeNCBI()
    prov = fetch_bundle(repo_root, dest=tmp_path / ACC, rnaseq=True, get_text=ncbi.get_text,
                        download=ncbi.download, today="2026-10-23")
    d = tmp_path / ACC
    assert (d / "genomic.fna").read_text() == GENOME
    assert (d / "chrom.sizes").read_text() == "chr1\t200\nchr2\t100\n"
    assert (d / "tx2gene.tsv").read_text() == "XM_1.1\t1\nXM_2.1\t1\nXR_3.1\t2\n"
    assert prov["annotation_release_matches_registry"] is True
    assert prov["genome_length"] == 300
    assert prov["feature_counts"] == {"genes": 2, "exons": 4, "introns": 2, "promoters": 2, "intergenic": 4}
    assert json.loads((d / "reference_provenance.json").read_text())["files"] == prov["files"]
    assert check_bundle(d) == []

    # A re-run reuses what is already there.
    ncbi.downloads.clear()
    fetch_bundle(repo_root, dest=d, rnaseq=True, get_text=ncbi.get_text, download=ncbi.download)
    assert ncbi.downloads == []

    # Drift is reported.
    (d / "features" / "genes.bed").write_text("edited\n")
    assert check_bundle(d) == ["features/genes.bed: SHA-256 differs from reference_provenance.json"]
    (d / "chrom.sizes").unlink()
    assert "chrom.sizes: missing" in check_bundle(d)


def test_corrupt_download_is_deleted(repo_root, tmp_path):
    ncbi = FakeNCBI(corrupt={f"{FOLDER}_genomic.gtf.gz"})
    with pytest.raises(ReferenceError, match="does not match NCBI"):
        fetch_bundle(repo_root, dest=tmp_path, get_text=ncbi.get_text, download=ncbi.download)
    assert not (tmp_path / f"{FOLDER}_genomic.gtf.gz").exists()


def test_reannotated_assembly_is_flagged(repo_root, tmp_path):
    ncbi = FakeNCBI(gtf_source="NCBI RefSeq GCF_963853765.1-RS_2026_01")
    prov = fetch_bundle(repo_root, dest=tmp_path, get_text=ncbi.get_text, download=ncbi.download)
    assert prov["annotation_release_matches_registry"] is False
    assert any("does not name the registered release" in p for p in check_bundle(tmp_path))


def test_check_without_provenance(tmp_path):
    assert "missing" in check_bundle(tmp_path)[0]


# --- assemblies and NCBI layout ---------------------------------------------

def test_resolve_assembly(repo_root):
    assert ref.resolve_assembly(repo_root)["ncbi_assembly_accession"] == ACC
    assert ref.resolve_assembly(repo_root, "xbMagGiga1.1")["ncbi_assembly_accession"] == ACC
    with pytest.raises(ReferenceError, match="no NCBI RefSeq annotation"):
        ref.resolve_assembly(repo_root, "cgigas_uk_roslin_v1")
    with pytest.raises(ReferenceError, match="not in data/reference"):
        ref.resolve_assembly(repo_root, "GCF_000000001.1")


def test_ncbi_layout_helpers():
    assert ref.accession_dir_url(ACC) == PARENT
    with pytest.raises(ReferenceError):
        ref.accession_dir_url("xbMagGiga1.1")
    with pytest.raises(ReferenceError, match="expected one folder"):
        ref.find_assembly_folder("<a href='x'>", ACC)
    assert ref.parse_md5_listing("abc  ./f.gz\n\nbad line here\n") == {"f.gz": "abc"}
    assert ref.annotation_release_matches("GCF_963853765.1-RS_2024_06", "NCBI RefSeq GCF_963853765.1-RS_2024_06")
    assert not ref.annotation_release_matches("GCF_963853765.1-RS_2024_06", None)


def test_plan_needs_every_file(repo_root):
    ncbi = FakeNCBI()
    ncbi.md5 = "\n".join(line for line in ncbi.md5.splitlines() if "gtf" not in line)
    with pytest.raises(ReferenceError, match="not listed"):
        ref.plan_bundle(repo_root, get_text=ncbi.get_text)


def test_default_dir_follows_whale_reference(repo_root, tmp_path, monkeypatch):
    monkeypatch.setenv("WHALE_REFERENCE", str(tmp_path))
    assert ref.default_bundle_dir(repo_root, ACC) == tmp_path / ACC
    monkeypatch.delenv("WHALE_REFERENCE")
    assert ref.default_bundle_dir(repo_root, ACC) == repo_root / "data" / "reference" / ACC


def test_cli_plan_and_fetch(repo_root, tmp_path, monkeypatch):
    ncbi = FakeNCBI()
    monkeypatch.setenv("WHALE_REPO", str(repo_root))
    monkeypatch.setenv("WHALE_REFERENCE", str(tmp_path))
    monkeypatch.setattr(ref, "plan_bundle", _bind(ref.plan_bundle, get_text=ncbi.get_text))
    monkeypatch.setattr(ref, "fetch_bundle", _bind(ref.fetch_bundle, get_text=ncbi.get_text, download=ncbi.download))

    res = CliRunner().invoke(main, ["fetch-reference", "--plan"])
    assert res.exit_code == 0, res.output
    assert "genomic.fna" in res.output and "genomic.gtf" in res.output and "rna.fna" not in res.output
    assert ncbi.downloads == []

    res = CliRunner().invoke(main, ["fetch-reference"])
    assert res.exit_code == 0, res.output
    assert "2 sequences" in res.output
    res = CliRunner().invoke(main, ["fetch-reference", "--check"])
    assert res.exit_code == 0 and "matches its provenance" in res.output


def _bind(fn, **fixed):
    def wrapper(*args, **kwargs):
        return fn(*args, **{**fixed, **kwargs})
    return wrapper


def test_download_failure_is_a_reference_error(tmp_path, monkeypatch):
    import whale.mirror as mirror

    def gone(url, dest, **kw):
        raise mirror.MirrorError(f"{url}: gave up after 6 attempts (timed out)")
    monkeypatch.setattr(mirror, "download", gone)
    with pytest.raises(ReferenceError, match="gave up"):
        ref._http_download("https://ftp.ncbi.nlm.nih.gov/x/genomic.fna.gz", tmp_path / "genomic.fna.gz")
