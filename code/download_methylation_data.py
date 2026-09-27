#!/usr/bin/env python3
"""
DNA Methylation Data Downloader

Downloads DNA methylation datasets for Crassostrea gigas and Magallana gigas
from NCBI SRA, based on the BioProjects documented in this repository.

How run discovery works
------------------------
For each BioProject the script asks NCBI which sequencing runs (SRR/ERR/DRR
accessions) belong to it, then downloads those runs with the SRA Toolkit.
Run discovery is done in this order:

  1. NCBI E-utilities over HTTPS (built in -- needs no extra software).
  2. NCBI Entrez Direct (`esearch`/`efetch`) if it is installed.

Every discovered run is then validated against its NCBI metadata: only
Crassostrea/Magallana gigas runs from a methylation assay (Bisulfite-Seq or
MeDIP-Seq) are kept. Anything else is logged and skipped, so a wrong or stale
BioProject accession is reported instead of being downloaded by mistake. Pass
--skip-validation to turn this filter off (not recommended).

The script never invents accession numbers. If it cannot reach NCBI it reports
the problem and stops, so you never download placeholder data by mistake.

Author: project-white-whale repository
"""

import os
import re
import sys
import csv
import io
import json
import math
import time
import shutil
import logging
import argparse
import subprocess
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional

# Configure logging. The log file is created lazily (delay=True) so that
# read-only commands such as --list and --help do not litter the working
# directory with an empty log.
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_methylation_data.log', delay=True),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# NCBI E-utilities endpoint. An API key (free, from an NCBI account) and a
# contact email raise the rate limit and are recommended but optional.
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
HTTP_TIMEOUT = 60          # seconds per HTTP request
HTTP_RETRIES = 3           # attempts before giving up on a request
RUN_ACCESSION_PREFIXES = ("SRR", "ERR", "DRR")
# BioProject accessions look like PRJNA123456 (NCBI), PRJEB12345 (ENA) or
# PRJDB1234 (DDBJ). Anything else is rejected before it reaches NCBI or a shell.
BIOPROJECT_RE = re.compile(r"^PRJ[NED][AB]\d+$")
# Marker written next to a run's FASTQ files once the download has been
# verified. Its presence is what makes a re-run skip the run.
DONE_SUFFIX = ".done"


def is_valid_bioproject(accession: str) -> bool:
    """True if ``accession`` is a well-formed BioProject accession."""
    return bool(BIOPROJECT_RE.match((accession or "").strip()))

# Runs are only accepted if they match the target species AND a methylation
# assay. This guard is what prevents the download of off-target data when a
# BioProject accession is wrong or points at unrelated samples. NCBI treats
# "Crassostrea gigas" and "Magallana gigas" as the same taxon, so both names
# appear in run metadata and both are accepted.
ALLOWED_ORGANISMS = ("crassostrea gigas", "magallana gigas")
# In SRA, both classic bisulfite and enzymatic-methyl (EM-seq) runs are
# labelled "Bisulfite-Seq"; MeDIP runs are labelled "MeDIP-Seq".
ALLOWED_STRATEGIES = ("bisulfite-seq", "medip-seq")

# Dataset catalogue.
#
# Every BioProject below was verified live against NCBI SRA: each one actually
# contains Crassostrea/Magallana gigas Bisulfite-Seq or MeDIP-Seq runs. Sample
# counts and sizes are the real run counts and summed SRA download sizes at the
# time of verification (2026-08); they can grow if submitters add runs. The
# script re-checks every accession against NCBI at run time and filters out any
# run that is not oyster methylation data, so a stale or wrong accession is
# reported rather than silently downloaded.
METHYLATION_DATASETS = {
    "wgbs_poms_adaptation": {
        "description": "Genetic/epigenetic adaptation to Pacific Oyster Mortality Syndrome (POMS)",
        "bioprojects": ["PRJEB60400"],
        "method": "WGBS",
        "estimated_samples": "246",
        "tissue_types": ["various"],
        "estimated_size_gb": "450-500",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJEB60400",
        "notes": "WGBS of oyster populations adapting to POMS (Crassostrea/Magallana gigas)"
    },
    "emseq_poms_gestinov": {
        "description": "Enzymatic methyl-seq of gill and mantle around POMS infection (GESTINOV 2021)",
        "bioprojects": ["PRJEB81880"],
        "method": "EM-seq (Bisulfite-Seq)",
        "estimated_samples": "40",
        "tissue_types": ["gill", "mantle"],
        "estimated_size_gb": "400-450",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJEB81880",
        "notes": "Enzymatic methyl sequencing before/after POMS infection"
    },
    "wgbs_aging_decicomp": {
        "description": "DNA methylation profiling across ages (DECICOMP)",
        "bioprojects": ["PRJEB105019"],
        "method": "WGBS",
        "estimated_samples": "60",
        "tissue_types": ["various"],
        "estimated_size_gb": "450-500",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJEB105019",
        "notes": "Epigenetic profiling of 4-, 16- and 28-month-old oysters under control conditions"
    },
    "wgbs_pesto": {
        "description": "Methylseq of Crassostrea gigas (PESTO project, 2022)",
        "bioprojects": ["PRJEB58545"],
        "method": "WGBS",
        "estimated_samples": "48",
        "tissue_types": ["various"],
        "estimated_size_gb": "350-400",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJEB58545",
        "notes": "Whole-genome methylation sequencing from the PESTO project"
    },
    "wgbs_transgen_infection": {
        "description": "Transgenerational infection-resistance plasticity after early microbial exposure",
        "bioprojects": ["PRJNA609264"],
        "method": "WGBS",
        "estimated_samples": "47",
        "tissue_types": ["various"],
        "estimated_size_gb": "600-650",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA609264",
        "notes": "Methylation underlying transgenerational adaptive phenotypic plasticity"
    },
    "wgbs_ph_ploidy": {
        "description": "WGBS of diploid and triploid ctenidia at different pH levels",
        "bioprojects": ["PRJNA682817"],
        "method": "WGBS",
        "estimated_samples": "24",
        "tissue_types": ["gill"],
        "estimated_size_gb": "100-110",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA682817",
        "notes": "DNA methylation response to pH in diploid vs triploid oysters"
    },
    "medip_development": {
        "description": "Developmental genome-wide methylome dynamics (MeDIP-seq)",
        "bioprojects": ["PRJNA324546"],
        "method": "MeDIP-seq",
        "estimated_samples": "21",
        "tissue_types": ["embryo", "larvae"],
        "estimated_size_gb": "5-10",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA324546",
        "notes": "MeDIP-seq across developmental stages"
    },
    "wgbs_epigenomics_series": {
        "description": "Magallana gigas epigenomics WGBS series",
        "bioprojects": ["PRJNA807732", "PRJNA562805", "PRJNA213124"],
        "method": "WGBS",
        "estimated_samples": "50",
        "tissue_types": ["various"],
        "estimated_size_gb": "300-330",
        "search_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA807732",
        "notes": "Three related 'Magallana gigas Epigenomics' bisulfite projects"
    }
}


def _eutils_params(extra: Dict[str, str]) -> Dict[str, str]:
    """Add shared identification params (api key, email, tool) to a request."""
    params = dict(extra)
    params["tool"] = "project-white-whale"
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    return params


def _http_get(endpoint: str, params: Dict[str, str]) -> Optional[str]:
    """GET a URL with retries and a polite delay. Returns text or None."""
    url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
    for attempt in range(1, HTTP_RETRIES + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "project-white-whale/1.0"})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            wait = 2 ** attempt
            logger.warning(f"NCBI request failed (attempt {attempt}/{HTTP_RETRIES}): {e}")
            if attempt < HTTP_RETRIES:
                time.sleep(wait)
    return None


def _parse_runs_from_runinfo(
    runinfo_csv: str,
    allowed_organisms: Optional[tuple] = ALLOWED_ORGANISMS,
    allowed_strategies: Optional[tuple] = ALLOWED_STRATEGIES,
) -> List[str]:
    """Extract run accessions from an SRA runinfo CSV.

    Only runs whose organism (ScientificName) and assay (LibraryStrategy) match
    the expected values are returned. This is the guard that stops the tool from
    downloading off-target data when a BioProject accession is wrong. Pass
    ``allowed_organisms=None`` and/or ``allowed_strategies=None`` to skip that
    check. Runs that are dropped are logged, grouped by organism/strategy, so a
    mismatched accession surfaces loudly instead of being fetched silently.
    """
    reader = csv.DictReader(io.StringIO(runinfo_csv))
    fields = reader.fieldnames or []
    if "Run" not in fields:
        return []

    check_org = bool(allowed_organisms) and "ScientificName" in fields
    check_strat = bool(allowed_strategies) and "LibraryStrategy" in fields
    org_set = tuple(o.lower() for o in (allowed_organisms or ()))
    strat_set = tuple(s.lower() for s in (allowed_strategies or ()))

    runs: List[str] = []
    skipped: Dict[tuple, int] = {}
    for row in reader:
        acc = (row.get("Run") or "").strip()
        if not acc.upper().startswith(RUN_ACCESSION_PREFIXES):
            continue
        organism = (row.get("ScientificName") or "").strip()
        strategy = (row.get("LibraryStrategy") or "").strip()
        org_ok = (not check_org) or organism.lower() in org_set
        strat_ok = (not check_strat) or strategy.lower() in strat_set
        if org_ok and strat_ok:
            runs.append(acc)
        else:
            key = (organism or "unknown", strategy or "unknown")
            skipped[key] = skipped.get(key, 0) + 1

    for (organism, strategy), n in sorted(skipped.items()):
        logger.warning(
            f"Skipped {n} off-target run(s) [organism='{organism}', "
            f"strategy='{strategy}']; expected an oyster methylation assay."
        )
    return runs


def fetch_runs_via_eutils(
    bioproject: str,
    allowed_organisms: Optional[tuple] = ALLOWED_ORGANISMS,
    allowed_strategies: Optional[tuple] = ALLOWED_STRATEGIES,
) -> Optional[List[str]]:
    """Discover SRA runs for a BioProject via NCBI E-utilities over HTTPS.

    Returns a list of accessions (possibly empty if the BioProject genuinely
    has none), or None if NCBI could not be reached at all.
    """
    logger.info(f"Querying NCBI for runs in {bioproject} (HTTPS E-utilities)")

    search = _http_get("esearch.fcgi", _eutils_params({
        "db": "sra",
        "term": f"{bioproject}[BioProject]",
        "retmax": "1",
        "usehistory": "y",
    }))
    if search is None:
        return None

    # Pull WebEnv / QueryKey / Count out of the XML without a heavy XML parser.
    def _between(text: str, tag: str) -> str:
        start = text.find(f"<{tag}>")
        end = text.find(f"</{tag}>")
        if start == -1 or end == -1:
            return ""
        return text[start + len(tag) + 2:end]

    count = _between(search, "Count")
    webenv = _between(search, "WebEnv")
    query_key = _between(search, "QueryKey")
    if not webenv or not query_key:
        logger.warning(f"No search history returned for {bioproject}; it may not exist.")
        return []
    if count == "0":
        logger.warning(f"BioProject {bioproject} has 0 SRA runs.")
        return []

    time.sleep(0.34 if not NCBI_API_KEY else 0.1)  # respect NCBI rate limits

    runinfo = _http_get("efetch.fcgi", _eutils_params({
        "db": "sra",
        "WebEnv": webenv,
        "query_key": query_key,
        "rettype": "runinfo",
        "retmode": "text",
    }))
    if runinfo is None:
        return None

    runs = _parse_runs_from_runinfo(runinfo, allowed_organisms, allowed_strategies)
    logger.info(f"Found {len(runs)} matching run(s) in {bioproject}")
    return runs


def fetch_runs_via_edirect(
    bioproject: str,
    allowed_organisms: Optional[tuple] = ALLOWED_ORGANISMS,
    allowed_strategies: Optional[tuple] = ALLOWED_STRATEGIES,
) -> Optional[List[str]]:
    """Discover runs using a local Entrez Direct install, if present."""
    if not (shutil.which("esearch") and shutil.which("efetch")):
        return None
    if not is_valid_bioproject(bioproject):
        logger.error(f"Refusing to query malformed BioProject accession: {bioproject!r}")
        return None
    logger.info(f"Querying {bioproject} via local Entrez Direct")
    # Run the two commands as a pipeline without a shell, so the accession is
    # passed as a plain argument and never interpreted by a shell.
    try:
        search = subprocess.run(
            ["esearch", "-db", "sra", "-query", f"{bioproject}[BioProject]"],
            capture_output=True, text=True, timeout=120,
        )
        if search.returncode != 0 or not search.stdout:
            return None
        fetch = subprocess.run(
            ["efetch", "-format", "runinfo"],
            input=search.stdout, capture_output=True, text=True, timeout=120,
        )
    except subprocess.TimeoutExpired:
        logger.warning(f"Entrez Direct timed out for {bioproject}")
        return None
    except OSError as e:
        logger.warning(f"Entrez Direct failed for {bioproject}: {e}")
        return None
    if fetch.returncode != 0 or not fetch.stdout:
        return None
    return _parse_runs_from_runinfo(fetch.stdout, allowed_organisms, allowed_strategies)


class MethylationDataDownloader:
    """Main class for downloading DNA methylation datasets."""

    def __init__(self, output_dir: str = "./methylation_data", max_parallel: int = 2,
                 validate: bool = True, run_timeout: Optional[int] = None):
        self.output_dir = Path(output_dir)
        self.max_parallel = max_parallel
        # Seconds allowed per run download; None means no limit. Large WGBS
        # runs on a slow link can take many hours, so there is no default cap.
        self.run_timeout = run_timeout if run_timeout and run_timeout > 0 else None

        # When validation is on (default), only oyster methylation runs are
        # accepted; setting these to None disables the respective check.
        self.allowed_organisms = ALLOWED_ORGANISMS if validate else None
        self.allowed_strategies = ALLOWED_STRATEGIES if validate else None
        if not validate:
            logger.warning(
                "Run validation is DISABLED (--skip-validation): runs will be "
                "downloaded regardless of organism or assay. Use with care."
            )

        # The SRA Toolkit is detected lazily, the first time a download is
        # attempted, so --list and --dry-run stay quiet and side-effect free.
        self._download_tool: Optional[str] = None
        self._tool_detected = False

    @property
    def download_tool(self) -> Optional[str]:
        """'fasterq-dump', 'fastq-dump', or None (detected on first use)."""
        if not self._tool_detected:
            self._download_tool = self._detect_download_tool()
            self._tool_detected = True
        return self._download_tool

    @download_tool.setter
    def download_tool(self, value: Optional[str]) -> None:
        self._download_tool = value
        self._tool_detected = True

    @staticmethod
    def _detect_download_tool() -> Optional[str]:
        """Return 'fasterq-dump', 'fastq-dump', or None."""
        for tool in ("fasterq-dump", "fastq-dump"):
            if shutil.which(tool):
                try:
                    result = subprocess.run([tool, '--version'],
                                            capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        logger.info(f"SRA Toolkit detected: {tool}")
                        return tool
                except (subprocess.TimeoutExpired, OSError):
                    continue
        logger.warning("SRA Toolkit not found. Install from: https://github.com/ncbi/sra-tools")
        return None

    def list_datasets(self) -> None:
        """List all available datasets."""
        print("\nAvailable DNA Methylation Datasets:")
        print("=" * 50)
        for dataset_id, info in METHYLATION_DATASETS.items():
            print(f"\nDataset ID: {dataset_id}")
            print(f"Description: {info['description']}")
            print(f"Method: {info['method']}")
            print(f"BioProjects: {', '.join(info['bioprojects'])}")
            print(f"Estimated Samples: {info['estimated_samples']}")
            print(f"Estimated Size: {info['estimated_size_gb']} GB")
            print(f"Tissue Types: {', '.join(info['tissue_types'])}")
            print(f"Notes: {info['notes']}")

    def get_bioproject_runs(self, bioproject: str) -> List[str]:
        """Get SRA run accessions for a BioProject from NCBI.

        Never returns fabricated data: if NCBI cannot be reached by any method,
        returns an empty list and logs the failure.
        """
        runs = fetch_runs_via_eutils(bioproject, self.allowed_organisms, self.allowed_strategies)
        if runs is None:
            logger.warning(f"HTTPS discovery unavailable for {bioproject}; trying Entrez Direct")
            runs = fetch_runs_via_edirect(bioproject, self.allowed_organisms, self.allowed_strategies)

        if runs is None:
            logger.error(
                f"Could not reach NCBI to look up runs for {bioproject}. "
                f"Check your internet connection and try again."
            )
            return []
        return runs

    @staticmethod
    def _estimated_upper_gb(dataset_info: Dict) -> Optional[int]:
        """Parse the upper bound of an 'X-Y' GB estimate, or None."""
        parts = str(dataset_info.get("estimated_size_gb", "")).split("-")
        try:
            return int(parts[-1])
        except (ValueError, IndexError):
            return None

    @classmethod
    def _estimated_needed_gb(cls, dataset_info: Dict, max_runs: Optional[int] = None,
                             n_bioprojects: Optional[int] = None) -> Optional[int]:
        """Estimate the disk space a download will need, in GB.

        Without a run limit this is the dataset's upper-bound estimate. With
        ``max_runs`` the estimate is scaled to the number of runs that will
        actually be fetched (per-run average x runs x BioProjects), so a small
        test download is not blocked by the size of the whole dataset.
        """
        upper = cls._estimated_upper_gb(dataset_info)
        if upper is None or not max_runs or max_runs <= 0:
            return upper
        try:
            total_runs = int(str(dataset_info.get("estimated_samples", "")).strip())
        except ValueError:
            return upper
        if total_runs <= 0:
            return upper
        n_bp = n_bioprojects or len(dataset_info.get("bioprojects", [])) or 1
        runs_to_fetch = min(max_runs * n_bp, total_runs)
        return min(upper, max(1, math.ceil(upper * runs_to_fetch / total_runs)))

    def _check_disk_space(self, needed_gb: Optional[int], force: bool) -> bool:
        """Warn (or block) if free disk space looks insufficient."""
        if needed_gb is None:
            return True
        self.output_dir.mkdir(parents=True, exist_ok=True)
        free_gb = shutil.disk_usage(self.output_dir).free / (1024 ** 3)
        logger.info(f"Free space at {self.output_dir}: {free_gb:.0f} GB (estimate needs ~{needed_gb} GB)")
        if free_gb < needed_gb:
            msg = (f"Only {free_gb:.0f} GB free but dataset may need up to {needed_gb} GB.")
            if force:
                logger.warning(msg + " Continuing because --force was given.")
                return True
            logger.error(msg + " Free up space or re-run with --force to override.")
            return False
        return True

    def download_dataset(self, dataset_id: str, bioproject: str = None,
                         max_runs: int = None, dry_run: bool = False,
                         force: bool = False) -> bool:
        """Download a specific dataset."""
        if dataset_id not in METHYLATION_DATASETS:
            logger.error(f"Unknown dataset ID: {dataset_id}. Use --list to see valid IDs.")
            return False

        dataset_info = METHYLATION_DATASETS[dataset_id]
        logger.info(f"Starting download for dataset: {dataset_info['description']}")

        if bioproject:
            bioproject = bioproject.strip()
            if not is_valid_bioproject(bioproject):
                logger.error(
                    f"'{bioproject}' is not a valid BioProject accession "
                    f"(expected e.g. PRJNA123456 or PRJEB12345)."
                )
                return False
            if bioproject not in dataset_info['bioprojects']:
                logger.error(
                    f"BioProject {bioproject} does not belong to dataset '{dataset_id}' "
                    f"(its BioProjects: {', '.join(dataset_info['bioprojects'])})."
                )
                return False
            bioprojects = [bioproject]
        else:
            bioprojects = list(dataset_info['bioprojects'])

        if not dry_run:
            needed = self._estimated_needed_gb(dataset_info, max_runs, len(bioprojects))
            if not self._check_disk_space(needed, force):
                return False

        if not dry_run and self.download_tool is None:
            logger.error(
                "SRA Toolkit is required to download data but was not found. "
                "Run ./install_dependencies.sh, or use --dry-run to preview."
            )
            return False

        dataset_dir = self.output_dir / dataset_id
        dataset_dir.mkdir(parents=True, exist_ok=True)

        metadata_file = dataset_dir / "dataset_info.json"
        with open(metadata_file, 'w') as f:
            json.dump(dataset_info, f, indent=2)

        success = True
        any_runs_found = False
        for bp in bioprojects:
            logger.info(f"Processing BioProject: {bp}")
            runs = self.get_bioproject_runs(bp)

            if not runs:
                logger.warning(f"No runs found for BioProject {bp}; skipping.")
                continue
            any_runs_found = True

            bp_dir = dataset_dir / bp
            bp_dir.mkdir(parents=True, exist_ok=True)
            # Record exactly which runs were discovered, for reproducibility.
            (bp_dir / "runs.txt").write_text("\n".join(runs) + "\n")

            if max_runs and len(runs) > max_runs:
                runs = runs[:max_runs]
                logger.info(f"Limiting to first {max_runs} run(s)")

            if dry_run:
                logger.info(f"DRY RUN: would download {len(runs)} run(s) from {bp}: {', '.join(runs)}")
                continue

            for run in runs:
                if not self._download_run(run, bp_dir):
                    logger.error(f"Failed to download run {run}")
                    success = False

        if not any_runs_found:
            logger.error("No runs were found for any BioProject in this dataset.")
            return False
        return success

    @staticmethod
    def _run_files(run_accession: str, output_dir: Path) -> List[Path]:
        """All FASTQ files (compressed or not) belonging to a run."""
        return sorted(
            list(output_dir.glob(f"{run_accession}*.fastq.gz"))
            + list(output_dir.glob(f"{run_accession}*.fastq"))
        )

    @staticmethod
    def _gzip_ok(path: Path) -> bool:
        """True if ``path`` is a complete, readable gzip file."""
        try:
            result = subprocess.run(['gzip', '-t', str(path)], capture_output=True, text=True)
            return result.returncode == 0
        except OSError:
            return False

    def _remove_run_files(self, run_accession: str, output_dir: Path) -> None:
        """Delete partial outputs of a run so a retry starts from scratch."""
        for p in self._run_files(run_accession, output_dir):
            try:
                p.unlink()
                logger.info(f"Removed partial file {p.name}")
            except OSError as e:
                logger.warning(f"Could not remove {p.name}: {e}")
        marker = output_dir / f"{run_accession}{DONE_SUFFIX}"
        if marker.exists():
            marker.unlink()

    def _is_run_complete(self, run_accession: str, output_dir: Path) -> bool:
        """Decide whether an earlier download of this run can be trusted.

        A run is complete when its ``.done`` marker exists. Downloads made by
        older versions of this tool have no marker: their compressed files are
        integrity-checked once with ``gzip -t`` and, if they pass, the marker
        is written so the check is not repeated. Files that fail the check are
        partial (an interrupted gzip, for example) and are removed so the run
        is fetched again.
        """
        marker = output_dir / f"{run_accession}{DONE_SUFFIX}"
        if marker.exists():
            return True
        files = self._run_files(run_accession, output_dir)
        if not files:
            return False
        gz = [p for p in files if p.suffix == ".gz"]
        plain = [p for p in files if p.suffix != ".gz"]
        if gz and not plain and all(p.stat().st_size > 0 and self._gzip_ok(p) for p in gz):
            logger.info(f"Run {run_accession}: existing files verified, marking complete")
            marker.touch()
            return True
        logger.warning(
            f"Run {run_accession}: found incomplete files from an earlier attempt; "
            f"removing them and downloading again"
        )
        self._remove_run_files(run_accession, output_dir)
        return False

    def _download_run(self, run_accession: str, output_dir: Path) -> bool:
        """Download and compress a single SRA run."""
        if self._is_run_complete(run_accession, output_dir):
            logger.info(f"Run {run_accession} already present, skipping")
            return True

        logger.info(f"Downloading run: {run_accession}")
        try:
            if self.download_tool == "fasterq-dump":
                cmd = ['fasterq-dump', '--split-files', '--outdir', str(output_dir),
                       '--progress', '--threads', '2', run_accession]
            else:
                cmd = ['fastq-dump', '--split-files', '--gzip',
                       '--outdir', str(output_dir), run_accession]

            logger.info(f"Executing: {' '.join(cmd)}")
            # Output is not captured: the SRA Toolkit's progress display goes
            # straight to the terminal so long downloads are visibly alive.
            result = subprocess.run(cmd, timeout=self.run_timeout)
            if result.returncode != 0:
                logger.error(
                    f"{self.download_tool} exited with code {result.returncode} for "
                    f"{run_accession} (see its output above)"
                )
                self._remove_run_files(run_accession, output_dir)
                return False

            # fasterq-dump writes uncompressed FASTQ; compress to save space.
            # gzip only removes the source once the .gz is complete, so an
            # interrupted compression leaves both files behind and the next
            # run detects that and redoes the run.
            if self.download_tool == "fasterq-dump":
                for fastq_file in output_dir.glob(f"{run_accession}*.fastq"):
                    logger.info(f"Compressing {fastq_file.name}")
                    try:
                        subprocess.run(['gzip', '-f', str(fastq_file)], check=True)
                    except (subprocess.CalledProcessError, OSError) as e:
                        logger.error(f"Could not gzip {fastq_file.name}: {e}")
                        self._remove_run_files(run_accession, output_dir)
                        return False

            if not self._validate_run(run_accession, output_dir):
                self._remove_run_files(run_accession, output_dir)
                return False
            (output_dir / f"{run_accession}{DONE_SUFFIX}").touch()
            return True

        except subprocess.TimeoutExpired:
            logger.error(
                f"Timeout downloading {run_accession} after {self.run_timeout} s; "
                f"removing partial files (raise or drop --run-timeout to allow longer)"
            )
            self._remove_run_files(run_accession, output_dir)
            return False
        except OSError as e:
            logger.error(f"Error downloading {run_accession}: {e}")
            self._remove_run_files(run_accession, output_dir)
            return False

    @classmethod
    def _validate_run(cls, run_accession: str, output_dir: Path) -> bool:
        """Confirm at least one non-empty, intact FASTQ file was produced."""
        produced = cls._run_files(run_accession, output_dir)
        non_empty = [p for p in produced if p.stat().st_size > 0]
        if not non_empty:
            logger.error(f"No non-empty FASTQ produced for {run_accession}")
            return False
        bad = [p for p in non_empty if p.suffix == ".gz" and not cls._gzip_ok(p)]
        if bad:
            logger.error(f"Corrupt gzip output for {run_accession}: {', '.join(p.name for p in bad)}")
            return False
        logger.info(f"Successfully downloaded {run_accession} "
                    f"({len(non_empty)} file(s), "
                    f"{sum(p.stat().st_size for p in non_empty) / 1e6:.1f} MB)")
        return True

    def create_download_script(self, dataset_ids: List[str],
                               output_file: str = "download_script.sh") -> None:
        """Create a standalone shell script that discovers and downloads runs.

        The generated script queries NCBI for the real runs in each BioProject
        at run time (using curl); it does not contain hardcoded accessions.
        """
        valid_ids = [d for d in dataset_ids if d in METHYLATION_DATASETS]
        unknown = [d for d in dataset_ids if d not in METHYLATION_DATASETS]
        for d in unknown:
            logger.warning(f"Skipping unknown dataset ID: {d}")
        if not valid_ids:
            logger.error("No valid dataset IDs provided; nothing to write.")
            return

        self.output_dir.mkdir(parents=True, exist_ok=True)
        script_path = self.output_dir / output_file
        lines: List[str] = []
        a = lines.append

        a("#!/bin/bash")
        a("# Auto-generated by download_methylation_data.py")
        a("# Discovers the real SRA runs for each BioProject from NCBI, then downloads them.")
        a("set -uo pipefail")
        a("")
        a('MAX_RUNS="${MAX_RUNS:-0}"   # 0 = all runs; set e.g. MAX_RUNS=3 to test')
        a('EUTILS="${EUTILS:-https://eutils.ncbi.nlm.nih.gov/entrez/eutils}"')
        a('NCBI_EMAIL="${NCBI_EMAIL:-}"; NCBI_API_KEY="${NCBI_API_KEY:-}"')
        a("")
        a("command_exists() { command -v \"$1\" >/dev/null 2>&1; }")
        a("")
        a("if ! command_exists curl; then echo \"Error: curl is required.\"; exit 1; fi")
        a("if ! command_exists python3; then echo \"Error: python3 is required (to parse NCBI run metadata).\"; exit 1; fi")
        a("if command_exists fasterq-dump; then")
        a('    DL="fasterq-dump"; DL_ARGS="--split-files --progress --threads 2"')
        a("elif command_exists fastq-dump; then")
        a('    DL="fastq-dump"; DL_ARGS="--split-files --gzip"')
        a("else")
        a('    echo "Error: SRA Toolkit not found. Install sra-tools first."; exit 1')
        a("fi")
        a('echo "Using $DL"; echo "Started: $(date)"')
        a("")
        a("# Extra identification parameters for NCBI (raise the rate limit if set).")
        a("ncbi_id_params() {")
        a('    local p=""')
        a('    [ -n "$NCBI_EMAIL" ] && p="${p}&email=${NCBI_EMAIL}"')
        a('    [ -n "$NCBI_API_KEY" ] && p="${p}&api_key=${NCBI_API_KEY}"')
        a('    echo "${p}&tool=project-white-whale"')
        a("}")
        a("")
        a("# Print the run accessions for a BioProject, one per line. Only oyster")
        a("# (Crassostrea/Magallana gigas) bisulfite/MeDIP runs are kept -- this")
        a("# mirrors the Python tool and stops off-target data being downloaded.")
        a("# The runinfo CSV is parsed with Python's csv module so quoted fields that")
        a("# contain commas do not shift the columns. Returns non-zero if NCBI could")
        a("# not be reached (as opposed to a project that simply has no runs).")
        a("# curl runs with -g so the [BioProject] field tag is not treated as a URL glob.")
        a("get_runs() {")
        a("    local bp=\"$1\"")
        a("    local hist webenv querykey")
        a('    if ! hist=$(curl -sfg "${EUTILS}/esearch.fcgi?db=sra&term=${bp}%5BBioProject%5D&usehistory=y&retmax=1$(ncbi_id_params)"); then')
        a('        echo "  Error: could not query NCBI for ${bp}" >&2; return 1')
        a("    fi")
        a("    webenv=$(echo \"$hist\" | sed -n 's:.*<WebEnv>\\(.*\\)</WebEnv>.*:\\1:p')")
        a("    querykey=$(echo \"$hist\" | sed -n 's:.*<QueryKey>\\(.*\\)</QueryKey>.*:\\1:p')")
        a('    [ -z "$webenv" ] && return 0')
        a('    sleep 0.4')
        a('    curl -sfg "${EUTILS}/efetch.fcgi?db=sra&WebEnv=${webenv}&query_key=${querykey}&rettype=runinfo&retmode=text$(ncbi_id_params)" \\')
        a("        | python3 -c '")
        a("import csv, sys")
        a("orgs = (\"crassostrea gigas\", \"magallana gigas\")")
        a("strats = (\"bisulfite-seq\", \"medip-seq\")")
        a("for row in csv.DictReader(sys.stdin):")
        a("    run = (row.get(\"Run\") or \"\").strip()")
        a("    org = (row.get(\"ScientificName\") or \"\").strip().lower()")
        a("    strat = (row.get(\"LibraryStrategy\") or \"\").strip().lower()")
        a("    if run[:3] in (\"SRR\", \"ERR\", \"DRR\") and org in orgs and strat in strats:")
        a("        print(run)")
        a("'")
        a("}")
        a("")
        a("# Download one run into a directory, compress it, and leave a .done marker")
        a("# so a re-run of this script skips it. Partial output is removed on failure.")
        a("download_run() {")
        a('    local run="$1" dir="$2"')
        a('    if [ -e "${dir}/${run}.done" ]; then echo "  ${run} already downloaded, skipping"; return 0; fi')
        a('    rm -f "${dir}/${run}"*.fastq "${dir}/${run}"*.fastq.gz')
        a('    echo "  Downloading ${run}..."')
        a('    if ! $DL $DL_ARGS --outdir "${dir}" "${run}"; then')
        a('        echo "  Error: ${DL} failed for ${run}" >&2')
        a('        rm -f "${dir}/${run}"*.fastq "${dir}/${run}"*.fastq.gz; return 1')
        a("    fi")
        a('    if [ "$DL" = "fasterq-dump" ]; then')
        a('        for f in "${dir}/${run}"*.fastq; do')
        a('            [ -e "$f" ] || continue')
        a('            if ! gzip -f "$f"; then echo "  Error: gzip failed for $f" >&2; rm -f "${dir}/${run}"*.fastq "${dir}/${run}"*.fastq.gz; return 1; fi')
        a("        done")
        a("    fi")
        a('    if ! ls "${dir}/${run}"*.fastq.gz >/dev/null 2>&1; then echo "  Error: no FASTQ produced for ${run}" >&2; return 1; fi')
        a('    touch "${dir}/${run}.done"')
        a("}")
        a("")
        a("FAILED=0")
        a("")

        total_size = 0
        for dataset_id in valid_ids:
            info = METHYLATION_DATASETS[dataset_id]
            upper = self._estimated_upper_gb(info)
            if upper:
                total_size += upper
            a(f"# Dataset: {info['description']} ({info['method']}, ~{info['estimated_size_gb']} GB)")
            a(f'echo "=== Dataset: {info["description"]} ==="')
            a(f'mkdir -p "{dataset_id}"')
            for bp in info['bioprojects']:
                a(f'echo "Discovering runs for {bp}..."')
                a(f'mkdir -p "{dataset_id}/{bp}"')
                a(f'if ! runs=$(get_runs "{bp}"); then')
                a('    FAILED=$((FAILED+1))')
                a(f'elif [ -z "$runs" ]; then echo "  No matching runs found for {bp}"; else')
                a(f'    printf "%s\\n" "$runs" > "{dataset_id}/{bp}/runs.txt"')
                a('    n=0')
                a('    for run in $runs; do')
                a('        n=$((n+1))')
                a('        if [ "$MAX_RUNS" -gt 0 ] && [ "$n" -gt "$MAX_RUNS" ]; then break; fi')
                a(f'        download_run "$run" "{dataset_id}/{bp}" || FAILED=$((FAILED+1))')
                a('    done')
                a('fi')
                a("")

        a('echo "Completed: $(date)"')
        a(f'echo "Estimated upper-bound size: ~{total_size} GB"')
        a('find . -name "*.fastq.gz" -type f | wc -l | xargs echo "Total compressed FASTQ files:"')
        a('du -sh */ 2>/dev/null | sort -h || true')
        a('if [ "$FAILED" -gt 0 ]; then echo "WARNING: $FAILED download(s)/lookup(s) failed; re-run this script to retry them." >&2; exit 1; fi')

        script_path.write_text("\n".join(lines) + "\n")
        os.chmod(script_path, 0o755)
        logger.info(f"Download script created: {script_path}")
        logger.info(f"Estimated upper-bound size: ~{total_size} GB")
        logger.info(f"Run it with: cd {self.output_dir} && ./{output_file}")
        logger.info("Tip: test first with  MAX_RUNS=2 ./" + output_file)


def main():
    parser = argparse.ArgumentParser(
        description="Download DNA methylation datasets for Crassostrea/Magallana gigas",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all available datasets
  python download_methylation_data.py --list

  # Preview what a dataset would download (writes only dataset_info.json
  # and runs.txt under the output directory; no FASTQ files are fetched)
  python download_methylation_data.py --dataset wgbs_ph_ploidy --dry-run

  # Download a specific dataset
  python download_methylation_data.py --dataset wgbs_ph_ploidy

  # Download a few runs of one BioProject to test your setup
  python download_methylation_data.py --dataset wgbs_ph_ploidy \\
      --bioproject PRJNA682817 --max-runs 3

  # Create a standalone shell script for several datasets
  python download_methylation_data.py --create-script \\
      --datasets wgbs_ph_ploidy medip_development

Tip: set NCBI_EMAIL (and optionally NCBI_API_KEY) in your environment to raise
NCBI rate limits.
        """
    )

    parser.add_argument('--list', action='store_true',
                        help='List all available datasets')
    parser.add_argument('--dataset', type=str,
                        help='Dataset ID to download')
    parser.add_argument('--bioproject', type=str,
                        help='Download only this BioProject of the dataset (optional; '
                             'must be one of the dataset\'s BioProjects)')
    parser.add_argument('--max-runs', type=int,
                        help='Maximum number of runs to download per BioProject')
    parser.add_argument('--output-dir', type=str, default='./methylation_data',
                        help='Output directory for downloaded data')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be downloaded without downloading')
    parser.add_argument('--force', action='store_true',
                        help='Proceed even if free disk space looks insufficient')
    parser.add_argument('--create-script', action='store_true',
                        help='Create a shell script for downloading datasets')
    parser.add_argument('--datasets', nargs='+',
                        help='List of dataset IDs for script generation')
    parser.add_argument('--max-parallel', type=int, default=2,
                        help='Reserved for future use; downloads currently run one at a time')
    parser.add_argument('--run-timeout', type=int, default=0,
                        help='Seconds allowed per run download before it is aborted '
                             '(0 = no limit, the default)')
    parser.add_argument('--skip-validation', action='store_true',
                        help='Do NOT filter discovered runs by organism/assay '
                             '(unsafe: may download non-oyster or non-methylation data)')

    args = parser.parse_args()

    downloader = MethylationDataDownloader(
        output_dir=args.output_dir,
        max_parallel=args.max_parallel,
        validate=not args.skip_validation,
        run_timeout=args.run_timeout,
    )

    if args.list:
        downloader.list_datasets()
        return

    if args.create_script:
        if not args.datasets:
            logger.error("--datasets is required when using --create-script")
            sys.exit(1)
        downloader.create_download_script(args.datasets)
        return

    if args.dataset:
        success = downloader.download_dataset(
            args.dataset,
            bioproject=args.bioproject,
            max_runs=args.max_runs,
            dry_run=args.dry_run,
            force=args.force,
        )
        if not success:
            sys.exit(1)
        return

    parser.print_help()


if __name__ == "__main__":
    main()
