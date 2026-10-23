"""Offline tests for `whale fetch-samplesheet` (src/whale/ena.py).

PRJNA173440 (Wang et al. 2014, the first oyster methylomes) is replayed from
ENA responses recorded on 2026-10-23 into tests/fixtures/ena/. The other cases
are small synthetic projects. No network access.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

import whale.ena as ena
from whale.cli import main
from whale.ena import ENAError, build_samplesheet, run_fastq_files

FIXTURES = Path(__file__).parent / "fixtures" / "ena"


def _replay(monkeypatch, run_tsv: str, sample_xml: str):
    def fake_get(url):
        if url.startswith(ena.ENA_PORTAL):
            return run_tsv
        if url.startswith(ena.ENA_XML):
            return sample_xml
        raise AssertionError(f"unexpected URL {url}")
    monkeypatch.setattr(ena, "_get", fake_get)


@pytest.fixture
def wang2014(monkeypatch):
    _replay(monkeypatch, (FIXTURES / "PRJNA173440_read_run.tsv").read_text(),
            (FIXTURES / "PRJNA173440_samples.xml").read_text())


def _read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- a real project, replayed -----------------------------------------------

def test_wang2014_sheet_manifest_and_provenance(wang2014, tmp_path):
    rep = build_samplesheet("PRJNA173440", study_id="WANG2014_FIRST_METHYLOMES_WGBS", out_dir=tmp_path,
                            condition_attributes=["sample_title"], repo=tmp_path)
    assert (rep["n_runs"], rep["n_libraries"], rep["n_biosamples"]) == (4, 4, 2)
    assert rep["group_sizes_libraries"] == {"Inbred": 2, "Wild": 2}
    assert rep["biosamples_with_several_libraries"] == {"SAMN01125790": 2, "SAMN01125791": 2}
    assert rep["library_layouts"] == ["PAIRED"]
    assert rep["runs_by_strategy_in_project"] == {"Bisulfite-Seq": 4}
    assert rep["total_fastq_bytes"] == 36392635708

    sheet = _read_tsv(tmp_path / "samplesheet.tsv")
    assert [r["sample_id"] for r in sheet] == ["SRR546471", "SRR546472", "SRR546473", "SRR546474"]
    assert list(sheet[0])[:4] == ["sample_id", "run_accession", "fastq_1", "fastq_2"]
    assert sheet[0]["fastq_1"] == "SRR546471_1.fastq.gz" and sheet[0]["fastq_2"] == "SRR546471_2.fastq.gz"
    assert {r["condition"] for r in sheet} == {"Inbred", "Wild"}
    assert sheet[0]["library_id"] == "SRX180559"

    manifest = _read_tsv(tmp_path / "fastq_manifest.tsv")
    assert len(manifest) == 8
    assert all(len(m["md5"]) == 32 and m["url"].startswith("ftp://") for m in manifest)
    assert {m["file"] for m in manifest} == {f for r in sheet for f in (r["fastq_1"], r["fastq_2"])}

    prov = json.loads((tmp_path / "ena_provenance.json").read_text())
    assert prov["outputs"]["samplesheet.tsv"] == ena.sha256_file(tmp_path / "samplesheet.tsv")
    assert "provenance_file" not in prov


def test_wang2014_missing_attribute_is_an_error(wang2014, tmp_path):
    with pytest.raises(ENAError, match="lack one of the condition attributes"):
        build_samplesheet("PRJNA173440", study_id="X2014_Y", out_dir=tmp_path,
                          condition_attributes=["treatment"])
    assert not (tmp_path / "samplesheet.tsv").exists()


def test_wang2014_labels_from_library_name(wang2014, tmp_path):
    rep = build_samplesheet("PRJNA173440", study_id="X2014_Y", out_dir=tmp_path,
                            condition_attributes=["origin"], label_field="library_name",
                            label_pattern=r"GSM\d+: (?P<origin>Wild|Inbred)$")
    assert rep["group_sizes_runs"] == {"Inbred": 2, "Wild": 2}
    assert rep["run_labels"]["field"] == "library_name"


def test_show_attributes(wang2014):
    summary = ena.sample_attribute_summary("PRJNA173440")
    assert summary["sample_title"] == ["Inbred", "Wild"]
    assert "(run) library_name" in summary


# --- synthetic projects -----------------------------------------------------

HEADER = ena.RUN_FIELDS


def _run(acc, sample, strategy="Bisulfite-Seq", layout="PAIRED", taxid="29159", **kw):
    n = 2 if layout == "PAIRED" else 1
    files = [f"ftp.sra.ebi.ac.uk/vol1/{acc}_{m}.fastq.gz" for m in (1, 2)][:n] if n == 2 else \
        [f"ftp.sra.ebi.ac.uk/vol1/{acc}.fastq.gz"]
    row = {k: "" for k in HEADER}
    row.update(run_accession=acc, sample_accession=sample, sample_alias=sample,
               experiment_accession=f"X{acc}", library_strategy=strategy, library_layout=layout,
               instrument_model="NovaSeq 6000", scientific_name="Magallana gigas", tax_id=taxid,
               read_count="1000", fastq_ftp=";".join(files), fastq_md5=";".join(["a" * 32] * n),
               fastq_bytes=";".join(["100"] * n))
    row.update(kw)
    return row


def _tsv(rows):
    return "\t".join(HEADER) + "\n" + "".join("\t".join(r[k] for k in HEADER) + "\n" for r in rows)


def _xml(samples: dict[str, dict]):
    parts = []
    for acc, attrs in samples.items():
        a = "".join(f"<SAMPLE_ATTRIBUTE><TAG>{k}</TAG><VALUE>{v}</VALUE></SAMPLE_ATTRIBUTE>" for k, v in attrs.items())
        parts.append(f'<SAMPLE accession="{acc}" alias="{acc}"><IDENTIFIERS><PRIMARY_ID>{acc}</PRIMARY_ID>'
                     f"</IDENTIFIERS><TITLE>{acc}</TITLE><SAMPLE_ATTRIBUTES>{a}</SAMPLE_ATTRIBUTES></SAMPLE>")
    return f"<SAMPLE_SET>{''.join(parts)}</SAMPLE_SET>"


def _project(monkeypatch, runs, samples):
    _replay(monkeypatch, _tsv(runs), _xml(samples))


def test_rnaseq_in_the_same_project_is_reported_not_kept(monkeypatch, tmp_path):
    runs = [_run(f"SRR{i}", f"S{i}") for i in range(4)] + \
           [_run(f"SRR{i}", f"S{i - 4}", strategy="RNA-Seq") for i in range(4, 8)]
    _project(monkeypatch, runs, {f"S{i}": {"treatment": "low" if i < 2 else "high"} for i in range(4)})
    rep = build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["treatment"])
    assert rep["n_runs"] == 4
    assert rep["runs_not_kept"] == {"RNA-Seq": 4}
    assert rep["group_sizes_libraries"] == {"high": 2, "low": 2}
    assert rep["min_group_libraries"] == 2


def test_explicit_strategy_must_exist(monkeypatch, tmp_path):
    _project(monkeypatch, [_run("SRR1", "S1")], {"S1": {"t": "a"}})
    with pytest.raises(ENAError, match="not in this project"):
        build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["t"],
                          library_strategies=["OTHER"])


def test_project_without_methylation_runs(monkeypatch, tmp_path):
    _project(monkeypatch, [_run("SRR1", "S1", strategy="RNA-Seq")], {"S1": {"t": "a"}})
    with pytest.raises(ENAError, match="no runs with library_strategy"):
        build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["t"])


def test_wrong_taxon_is_an_error(monkeypatch, tmp_path):
    _project(monkeypatch, [_run("SRR1", "S1"), _run("SRR2", "S2", taxid="6565", scientific_name="Crassostrea virginica")],
             {"S1": {"t": "a"}, "S2": {"t": "a"}})
    with pytest.raises(ENAError, match="not Magallana"):
        build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["t"])


def test_single_end_sheet_has_empty_fastq_2(monkeypatch, tmp_path):
    _project(monkeypatch, [_run("SRR1", "S1", layout="SINGLE")], {"S1": {"t": "a"}})
    build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["t"])
    row = _read_tsv(tmp_path / "samplesheet.tsv")[0]
    assert row["fastq_1"] == "SRR1.fastq.gz" and row["fastq_2"] == ""


def test_lanes_of_one_library_share_library_id(monkeypatch, tmp_path):
    runs = [_run("SRR1", "S1", experiment_accession="SRX1"), _run("SRR2", "S1", experiment_accession="SRX1")]
    _project(monkeypatch, runs, {"S1": {"t": "a"}})
    rep = build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["t"])
    assert rep["n_libraries"] == 1 and rep["libraries_with_several_runs"] == {"SRX1": 2}


def test_attribute_clashing_with_sheet_column(monkeypatch, tmp_path):
    _project(monkeypatch, [_run("SRR1", "S1")], {"S1": {"condition": "a"}})
    with pytest.raises(ENAError, match="collide"):
        build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path, condition_attributes=["condition"])


def test_design_table_join(monkeypatch, tmp_path):
    _project(monkeypatch, [_run("SRR1", "S1"), _run("SRR2", "S2")],
             {"S1": {"isolate": "A1"}, "S2": {"isolate": "B7"}})
    table = tmp_path / "design.tsv"
    table.write_text("isolate\tploidy\nA1\tdiploid\nB7\ttriploid\n")
    rep = build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path / "out",
                            condition_attributes=["ploidy"], design_table=table, design_join="isolate",
                            design_source="authors' sheet", repo=tmp_path)
    assert rep["group_sizes_runs"] == {"diploid": 1, "triploid": 1}
    assert rep["design_table"]["path"] == "design.tsv"

    table.write_text("isolate\tploidy\nA1\tdiploid\n")
    with pytest.raises(ENAError, match="no row"):
        build_samplesheet("PRJNA1", study_id="X2020_Y", out_dir=tmp_path / "out",
                          condition_attributes=["ploidy"], design_table=table, design_join="isolate")


# --- FASTQ file selection ---------------------------------------------------

def test_submitted_fastq_used_when_ena_has_none():
    run = {"library_layout": "PAIRED", "fastq_ftp": "", "submitted_format": "FASTQ;FASTQ",
           "submitted_ftp": "h/x_R2.fq.gz;h/x_R1.fq.gz", "submitted_md5": "m2;m1", "submitted_bytes": "2;1"}
    source, files, dropped = run_fastq_files(run)
    assert source == "submitted_ftp"
    assert [f[0] for f in files] == ["h/x_R1.fq.gz", "h/x_R2.fq.gz"]
    assert [f[1] for f in files] == ["m1", "m2"]


def test_bam_only_run_has_no_fastq():
    run = {"library_layout": "PAIRED", "fastq_ftp": "", "submitted_format": "BAM", "submitted_ftp": "h/x.bam"}
    assert run_fastq_files(run) == (None, [], [])


def test_unpaired_third_file_dropped():
    run = {"library_layout": "PAIRED", "fastq_ftp": "h/R.fastq.gz;h/R_1.fastq.gz;h/R_2.fastq.gz",
           "fastq_md5": "a;b;c", "fastq_bytes": "1;2;3"}
    source, files, dropped = run_fastq_files(run)
    assert [f[0] for f in files] == ["h/R_1.fastq.gz", "h/R_2.fastq.gz"] and dropped == ["h/R.fastq.gz"]


# --- CLI --------------------------------------------------------------------

def test_cli_fetch_samplesheet(wang2014, tmp_path, repo_root, monkeypatch):
    monkeypatch.setenv("WHALE_REPO", str(repo_root))
    out = tmp_path / "study"
    res = CliRunner().invoke(main, ["fetch-samplesheet", "--bioproject", "PRJNA173440",
                                    "--study", "WANG2014_FIRST_METHYLOMES_WGBS",
                                    "--condition-attribute", "sample_title", "--out", str(out)])
    assert res.exit_code == 0, res.output
    assert "4 runs, 4 libraries, 2 BioSamples" in res.output
    assert "shared_biosample" in res.output and "low_replication" in res.output
    assert (out / "samplesheet.tsv").exists()


def test_cli_requires_condition_or_show(wang2014, repo_root, monkeypatch):
    monkeypatch.setenv("WHALE_REPO", str(repo_root))
    res = CliRunner().invoke(main, ["fetch-samplesheet", "--bioproject", "PRJNA173440", "--study", "X2014_Y"])
    assert res.exit_code != 0 and "--condition-attribute is required" in res.output
    res = CliRunner().invoke(main, ["fetch-samplesheet", "--bioproject", "PRJNA173440", "--show-attributes"])
    assert res.exit_code == 0 and "sample_title" in res.output
