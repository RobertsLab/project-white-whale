"""Unit and integration tests for download_methylation_data.py.

Run from the repository root with:  python -m pytest code/tests -q

No network access is needed: NCBI and the SRA Toolkit are replaced by small
fake executables placed on PATH inside a temporary directory.
"""

import gzip
import importlib.util
import logging
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

CODE_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CODE_DIR / "download_methylation_data.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("download_methylation_data", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # The module installs a lazily-created file handler in the current
    # directory; drop it so tests never write a log into the repository.
    root = logging.getLogger()
    for handler in list(root.handlers):
        if isinstance(handler, logging.FileHandler):
            root.removeHandler(handler)
    return module


dl = _load_module()


# A runinfo CSV with the real column set. Row 1 is an on-target run whose
# LibraryName holds a quoted comma (this used to break the awk-based filter).
# Row 2 is an oyster RNA-seq run (wrong assay); row 3 is human bisulfite
# (wrong organism); row 4 uses the Magallana name and an ERR accession.
HEADER = (
    "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,"
    "AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,"
    "LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,"
    "Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,"
    "SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,"
    "g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,"
    "Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,"
    "dbgap_study_accession,Consent,RunHash,ReadHash"
)


def _row(run, libname, strategy, taxid, organism):
    return (
        f"{run},2020,2020,1,1,1,1,1,,,SRX1,{libname},{strategy},RANDOM,GENOMIC,"
        f"PAIRED,0,0,ILLUMINA,x,SRP1,PRJNA1,,1,SRS1,SAMN1,simple,{taxid},"
        f"{organism},s1,,,,,,,,,,,,C,SUB,,public,h,h"
    )


RUNINFO_CSV = "\n".join([
    HEADER,
    _row("SRR1000001", '"lib, with comma"', "Bisulfite-Seq", 29159, "Crassostrea gigas"),
    _row("SRR1000002", "lib2", "RNA-Seq", 29159, "Crassostrea gigas"),
    _row("SRR1000003", "lib3", "Bisulfite-Seq", 9606, "Homo sapiens"),
    _row("ERR2000004", "lib4", "MeDIP-Seq", 29159, "Magallana gigas"),
]) + "\n"


# --------------------------------------------------------------------------
# Pure functions
# --------------------------------------------------------------------------

def test_parse_runs_keeps_only_oyster_methylation_runs():
    assert dl._parse_runs_from_runinfo(RUNINFO_CSV) == ["SRR1000001", "ERR2000004"]


def test_parse_runs_without_validation_keeps_everything():
    assert dl._parse_runs_from_runinfo(RUNINFO_CSV, None, None) == [
        "SRR1000001", "SRR1000002", "SRR1000003", "ERR2000004",
    ]


def test_parse_runs_handles_missing_run_column():
    assert dl._parse_runs_from_runinfo("Foo,Bar\n1,2\n") == []


@pytest.mark.parametrize("acc,ok", [
    ("PRJNA682817", True),
    ("PRJEB60400", True),
    ("PRJDB1234", True),
    (" PRJNA1 ", True),
    ("PRJNA", False),
    ("SRR123", False),
    ("PRJNA1; rm -rf /", False),
    ("PRJNA1[BioProject]", False),
    ("", False),
    (None, False),
])
def test_is_valid_bioproject(acc, ok):
    assert dl.is_valid_bioproject(acc) is ok


def test_estimated_needed_gb_full_dataset_uses_upper_bound():
    info = dl.METHYLATION_DATASETS["wgbs_ph_ploidy"]  # 100-110 GB, 24 runs
    assert dl.MethylationDataDownloader._estimated_needed_gb(info) == 110
    assert dl.MethylationDataDownloader._estimated_needed_gb(info, max_runs=0) == 110


def test_estimated_needed_gb_scales_with_max_runs():
    info = dl.METHYLATION_DATASETS["wgbs_ph_ploidy"]
    # 2 of 24 runs at 110 GB -> ceil(9.17) = 10 GB
    assert dl.MethylationDataDownloader._estimated_needed_gb(info, max_runs=2, n_bioprojects=1) == 10
    # Never more than the whole dataset.
    assert dl.MethylationDataDownloader._estimated_needed_gb(info, max_runs=500) == 110


def test_estimated_needed_gb_counts_every_bioproject():
    info = dl.METHYLATION_DATASETS["wgbs_epigenomics_series"]  # 3 projects, 50 runs, 330 GB
    one = dl.MethylationDataDownloader._estimated_needed_gb(info, max_runs=2, n_bioprojects=1)
    three = dl.MethylationDataDownloader._estimated_needed_gb(info, max_runs=2)
    assert one < three <= 330


def test_catalog_totals_match_documentation():
    runs = sum(int(d["estimated_samples"]) for d in dl.METHYLATION_DATASETS.values())
    upper = sum(dl.MethylationDataDownloader._estimated_upper_gb(d) for d in dl.METHYLATION_DATASETS.values())
    projects = [bp for d in dl.METHYLATION_DATASETS.values() for bp in d["bioprojects"]]
    assert runs == 536
    assert upper == 2950
    assert len(projects) == 10 and len(set(projects)) == 10
    assert all(dl.is_valid_bioproject(bp) for bp in projects)


# --------------------------------------------------------------------------
# Downloader behaviour (no network, no SRA toolkit)
# --------------------------------------------------------------------------

@pytest.fixture
def downloader(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return dl.MethylationDataDownloader(output_dir=str(tmp_path / "out"))


def test_constructor_has_no_side_effects(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    dl.MethylationDataDownloader(output_dir=str(tmp_path / "out"))
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob("*.log"))


def test_download_dataset_rejects_unknown_dataset(downloader):
    assert downloader.download_dataset("nope", dry_run=True) is False


def test_download_dataset_rejects_malformed_bioproject(downloader, monkeypatch):
    called = []
    monkeypatch.setattr(downloader, "get_bioproject_runs", lambda bp: called.append(bp) or [])
    assert downloader.download_dataset("wgbs_ph_ploidy", bioproject="PRJNA1; echo", dry_run=True) is False
    assert called == []


def test_download_dataset_rejects_foreign_bioproject(downloader, monkeypatch):
    called = []
    monkeypatch.setattr(downloader, "get_bioproject_runs", lambda bp: called.append(bp) or [])
    assert downloader.download_dataset("wgbs_ph_ploidy", bioproject="PRJNA324546", dry_run=True) is False
    assert called == []


def test_dry_run_records_runs_and_fetches_nothing(downloader, monkeypatch):
    monkeypatch.setattr(downloader, "get_bioproject_runs", lambda bp: ["SRR1", "SRR2", "SRR3"])
    monkeypatch.setattr(downloader, "_download_run", lambda *a: pytest.fail("must not download"))
    assert downloader.download_dataset("wgbs_ph_ploidy", max_runs=2, dry_run=True) is True
    bp_dir = downloader.output_dir / "wgbs_ph_ploidy" / "PRJNA682817"
    assert (bp_dir / "runs.txt").read_text().split() == ["SRR1", "SRR2", "SRR3"]
    assert (downloader.output_dir / "wgbs_ph_ploidy" / "dataset_info.json").exists()


def test_disk_check_is_scaled_for_small_test_downloads(downloader, monkeypatch):
    """A 2-run test of a 110 GB dataset must not be blocked on a 50 GB disk."""
    class Usage:
        free = 50 * 1024 ** 3
    monkeypatch.setattr(dl.shutil, "disk_usage", lambda p: Usage)
    monkeypatch.setattr(downloader, "get_bioproject_runs", lambda bp: ["SRR1", "SRR2", "SRR3"])
    downloaded = []
    monkeypatch.setattr(downloader, "_download_run", lambda run, d: downloaded.append(run) or True)
    downloader.download_tool = "fasterq-dump"

    assert downloader.download_dataset("wgbs_ph_ploidy", max_runs=2) is True
    assert downloaded == ["SRR1", "SRR2"]
    # The full dataset is still refused without --force.
    downloaded.clear()
    assert downloader.download_dataset("wgbs_ph_ploidy") is False
    assert downloaded == []
    assert downloader.download_dataset("wgbs_ph_ploidy", force=True) is True


def _write_gz(path: Path, payload: bytes = b"@r1\nACGT\n+\nIIII\n"):
    with gzip.open(path, "wb") as fh:
        fh.write(payload)


def test_is_run_complete_trusts_marker(downloader, tmp_path):
    d = tmp_path / "bp"
    d.mkdir()
    (d / "SRR9.done").touch()
    assert downloader._is_run_complete("SRR9", d) is True


def test_is_run_complete_verifies_legacy_files_and_writes_marker(downloader, tmp_path):
    d = tmp_path / "bp"
    d.mkdir()
    _write_gz(d / "SRR9_1.fastq.gz")
    _write_gz(d / "SRR9_2.fastq.gz")
    assert downloader._is_run_complete("SRR9", d) is True
    assert (d / "SRR9.done").exists()


def test_is_run_complete_removes_truncated_gzip(downloader, tmp_path):
    d = tmp_path / "bp"
    d.mkdir()
    _write_gz(d / "SRR9_1.fastq.gz", b"x" * 5000)
    data = (d / "SRR9_1.fastq.gz").read_bytes()
    (d / "SRR9_1.fastq.gz").write_bytes(data[: len(data) // 2])  # truncate
    assert downloader._is_run_complete("SRR9", d) is False
    assert not list(d.iterdir())


def test_is_run_complete_treats_leftover_plain_fastq_as_partial(downloader, tmp_path):
    d = tmp_path / "bp"
    d.mkdir()
    _write_gz(d / "SRR9_1.fastq.gz")
    (d / "SRR9_1.fastq").write_text("@r1\nACGT\n+\nIIII\n")  # interrupted gzip
    assert downloader._is_run_complete("SRR9", d) is False
    assert not list(d.iterdir())


def _fake_tool(bin_dir: Path, name: str, body: str) -> None:
    path = bin_dir / name
    path.write_text("#!/bin/bash\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def fake_bin(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return bin_dir


FAKE_FASTERQ_DUMP = r'''
# Minimal stand-in for fasterq-dump: writes two split FASTQ files.
out="."
while [ $# -gt 1 ]; do
  case "$1" in --outdir) out="$2"; shift 2;; *) shift;; esac
done
run="$1"
[ "${FAKE_FAIL:-0}" = "1" ] && { echo "boom" >&2; printf 'partial' > "$out/${run}_1.fastq"; exit 3; }
printf '@%s.1\nACGT\n+\nIIII\n' "$run" > "$out/${run}_1.fastq"
printf '@%s.1\nTGCA\n+\nIIII\n' "$run" > "$out/${run}_2.fastq"
'''


def test_download_run_end_to_end_with_fake_toolkit(downloader, fake_bin, tmp_path):
    _fake_tool(fake_bin, "fasterq-dump", FAKE_FASTERQ_DUMP)
    downloader.download_tool = "fasterq-dump"
    d = tmp_path / "bp"
    d.mkdir()
    assert downloader._download_run("SRR42", d) is True
    names = sorted(p.name for p in d.iterdir())
    assert names == ["SRR42.done", "SRR42_1.fastq.gz", "SRR42_2.fastq.gz"]
    with gzip.open(d / "SRR42_1.fastq.gz", "rt") as fh:
        assert fh.readline().startswith("@SRR42")
    # Second call is a no-op.
    assert downloader._download_run("SRR42", d) is True


def test_download_run_failure_cleans_up(downloader, fake_bin, tmp_path, monkeypatch):
    _fake_tool(fake_bin, "fasterq-dump", FAKE_FASTERQ_DUMP)
    monkeypatch.setenv("FAKE_FAIL", "1")
    downloader.download_tool = "fasterq-dump"
    d = tmp_path / "bp"
    d.mkdir()
    assert downloader._download_run("SRR43", d) is False
    assert not list(d.iterdir())


def test_download_run_timeout_cleans_up(downloader, fake_bin, tmp_path):
    _fake_tool(fake_bin, "fasterq-dump", 'printf "x" > "$3/${@: -1}_1.fastq"; sleep 5\n')
    downloader.download_tool = "fasterq-dump"
    downloader.run_timeout = 1
    d = tmp_path / "bp"
    d.mkdir()
    assert downloader._download_run("SRR44", d) is False
    assert not list(d.iterdir())


# --------------------------------------------------------------------------
# Generated standalone shell script
# --------------------------------------------------------------------------

FAKE_CURL = r'''
# Stand-in for curl that serves canned NCBI responses. The URL is the last arg.
url="${@: -1}"
case "$url" in
  *"[BioProject]"*) echo "curl: (3) bad range in URL position" >&2; exit 3;;
  *esearch.fcgi*) echo "<eSearchResult><Count>4</Count><RetMax>1</RetMax><QueryKey>1</QueryKey><WebEnv>MCID_TEST</WebEnv></eSearchResult>";;
  *efetch.fcgi*) cat "$FAKE_RUNINFO";;
  *) exit 1;;
esac
'''


def test_generated_script_discovers_and_downloads(downloader, fake_bin, tmp_path, monkeypatch):
    runinfo = tmp_path / "runinfo.csv"
    runinfo.write_text(RUNINFO_CSV)
    _fake_tool(fake_bin, "curl", FAKE_CURL)
    _fake_tool(fake_bin, "fasterq-dump", FAKE_FASTERQ_DUMP)
    monkeypatch.setenv("FAKE_RUNINFO", str(runinfo))

    downloader.create_download_script(["wgbs_ph_ploidy", "not_a_dataset"])
    script = downloader.output_dir / "download_script.sh"
    text = script.read_text()
    # Regression checks for the URL-globbing bug and CSV parsing.
    curl_lines = [line for line in text.splitlines() if "curl " in line and "esearch" in line]
    assert curl_lines and all("-sfg" in line and "%5BBioProject%5D" in line for line in curl_lines)
    assert all("[BioProject]" not in line for line in curl_lines)
    assert "csv.DictReader" in text
    assert "not_a_dataset" not in text
    assert subprocess.run(["bash", "-n", str(script)]).returncode == 0

    env = dict(os.environ, MAX_RUNS="1")
    result = subprocess.run(["bash", str(script)], cwd=downloader.output_dir, env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    bp_dir = downloader.output_dir / "wgbs_ph_ploidy" / "PRJNA682817"
    assert (bp_dir / "runs.txt").read_text().split() == ["SRR1000001", "ERR2000004"]
    assert sorted(p.name for p in bp_dir.iterdir()) == [
        "SRR1000001.done", "SRR1000001_1.fastq.gz", "SRR1000001_2.fastq.gz", "runs.txt",
    ]

    # Re-running skips the finished run.
    result = subprocess.run(["bash", str(script)], cwd=downloader.output_dir, env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0
    assert "already downloaded, skipping" in result.stdout


def test_generated_script_reports_ncbi_failure(downloader, fake_bin, tmp_path):
    _fake_tool(fake_bin, "curl", "exit 7\n")
    _fake_tool(fake_bin, "fasterq-dump", FAKE_FASTERQ_DUMP)
    downloader.create_download_script(["medip_development"])
    script = downloader.output_dir / "download_script.sh"
    result = subprocess.run(["bash", str(script)], cwd=downloader.output_dir,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 1
    assert "could not query NCBI" in result.stderr


# --------------------------------------------------------------------------
# Entrez Direct fallback never goes through a shell
# --------------------------------------------------------------------------

def test_edirect_passes_accession_as_argument(fake_bin, tmp_path, monkeypatch):
    runinfo = tmp_path / "runinfo.csv"
    runinfo.write_text(RUNINFO_CSV)
    log = tmp_path / "esearch.log"
    _fake_tool(fake_bin, "esearch", f'printf "%s\\n" "$@" > "{log}"; echo "<ENTREZ_DIRECT/>"\n')
    _fake_tool(fake_bin, "efetch", f'cat "{runinfo}"\n')
    assert dl.fetch_runs_via_edirect("PRJNA682817") == ["SRR1000001", "ERR2000004"]
    assert log.read_text().split("\n")[:5] == ["-db", "sra", "-query", "PRJNA682817[BioProject]", ""]
    assert dl.fetch_runs_via_edirect("PRJNA1; touch pwned") is None
    assert not (tmp_path / "pwned").exists()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
