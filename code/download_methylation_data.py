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
import sys
import csv
import io
import json
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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_methylation_data.log'),
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
    logger.info(f"Querying {bioproject} via local Entrez Direct")
    try:
        cmd = f'esearch -db sra -query "{bioproject}[BioProject]" | efetch -format runinfo'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        logger.warning(f"Entrez Direct timed out for {bioproject}")
        return None
    if result.returncode != 0 or not result.stdout:
        return None
    return _parse_runs_from_runinfo(result.stdout, allowed_organisms, allowed_strategies)


class MethylationDataDownloader:
    """Main class for downloading DNA methylation datasets."""

    def __init__(self, output_dir: str = "./methylation_data", max_parallel: int = 2,
                 validate: bool = True):
        self.output_dir = Path(output_dir)
        self.max_parallel = max_parallel
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # When validation is on (default), only oyster methylation runs are
        # accepted; setting these to None disables the respective check.
        self.allowed_organisms = ALLOWED_ORGANISMS if validate else None
        self.allowed_strategies = ALLOWED_STRATEGIES if validate else None
        if not validate:
            logger.warning(
                "Run validation is DISABLED (--skip-validation): runs will be "
                "downloaded regardless of organism or assay. Use with care."
            )

        # Detect the SRA Toolkit downloader once, up front.
        self.download_tool = self._detect_download_tool()

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

    def _check_disk_space(self, needed_gb: Optional[int], force: bool) -> bool:
        """Warn (or block) if free disk space looks insufficient."""
        if needed_gb is None:
            return True
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

        if not dry_run and not self._check_disk_space(self._estimated_upper_gb(dataset_info), force):
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

        bioprojects = [bioproject] if bioproject else dataset_info['bioprojects']

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

    def _download_run(self, run_accession: str, output_dir: Path) -> bool:
        """Download and compress a single SRA run."""
        existing = list(output_dir.glob(f"{run_accession}*.fastq.gz"))
        if existing:
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
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            if result.returncode != 0:
                logger.error(f"{self.download_tool} failed for {run_accession}: {result.stderr.strip()}")
                return False

            # fasterq-dump writes uncompressed FASTQ; compress to save space.
            if self.download_tool == "fasterq-dump":
                for fastq_file in output_dir.glob(f"{run_accession}*.fastq"):
                    logger.info(f"Compressing {fastq_file.name}")
                    try:
                        subprocess.run(['gzip', '-f', str(fastq_file)], check=True, timeout=1800)
                    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                        logger.warning(f"Could not gzip {fastq_file.name}: {e}")

            return self._validate_run(run_accession, output_dir)

        except subprocess.TimeoutExpired:
            logger.error(f"Timeout downloading {run_accession}")
            return False
        except OSError as e:
            logger.error(f"Error downloading {run_accession}: {e}")
            return False

    @staticmethod
    def _validate_run(run_accession: str, output_dir: Path) -> bool:
        """Confirm at least one non-empty FASTQ file was produced."""
        produced = list(output_dir.glob(f"{run_accession}*.fastq.gz")) + \
            list(output_dir.glob(f"{run_accession}*.fastq"))
        non_empty = [p for p in produced if p.stat().st_size > 0]
        if not non_empty:
            logger.error(f"No non-empty FASTQ produced for {run_accession}")
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

        script_path = self.output_dir / output_file
        lines: List[str] = []
        a = lines.append

        a("#!/bin/bash")
        a("# Auto-generated by download_methylation_data.py")
        a("# Discovers the real SRA runs for each BioProject from NCBI, then downloads them.")
        a("set -euo pipefail")
        a("")
        a('MAX_RUNS="${MAX_RUNS:-0}"   # 0 = all runs; set e.g. MAX_RUNS=3 to test')
        a("")
        a("command_exists() { command -v \"$1\" >/dev/null 2>&1; }")
        a("")
        a("if ! command_exists curl; then echo \"Error: curl is required.\"; exit 1; fi")
        a("if command_exists fasterq-dump; then")
        a('    DL="fasterq-dump"; DL_ARGS="--split-files --progress --threads 2"')
        a("elif command_exists fastq-dump; then")
        a('    DL="fastq-dump"; DL_ARGS="--split-files --gzip"')
        a("else")
        a('    echo "Error: SRA Toolkit not found. Install sra-tools first."; exit 1')
        a("fi")
        a('echo "Using $DL"; echo "Started: $(date)"')
        a("")
        a("# Print the run accessions for a BioProject, one per line. Only oyster")
        a("# (Crassostrea/Magallana gigas) bisulfite/MeDIP runs are kept -- this")
        a("# mirrors the Python tool and stops off-target data being downloaded.")
        a("get_runs() {")
        a("    local bp=\"$1\"")
        a('    local base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"')
        a('    local hist; hist=$(curl -s "${base}/esearch.fcgi?db=sra&term=${bp}[BioProject]&usehistory=y&retmax=1")')
        a("    local webenv querykey")
        a("    webenv=$(echo \"$hist\" | sed -n 's:.*<WebEnv>\\(.*\\)</WebEnv>.*:\\1:p')")
        a("    querykey=$(echo \"$hist\" | sed -n 's:.*<QueryKey>\\(.*\\)</QueryKey>.*:\\1:p')")
        a('    [ -z "$webenv" ] && return 0')
        a('    curl -s "${base}/efetch.fcgi?db=sra&WebEnv=${webenv}&query_key=${querykey}&rettype=runinfo&retmode=text" \\')
        a("        | awk -F, 'NR==1{for(i=1;i<=NF;i++){if($i==\"Run\")r=i; if($i==\"ScientificName\")s=i; if($i==\"LibraryStrategy\")l=i}; next}")
        a("                   $r ~ /^[SED]RR/ && tolower($s) ~ /gigas/ && tolower($l) ~ /bisulfite|medip/ {print $r}'")
        a("}")
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
                a(f'runs=$(get_runs "{bp}")')
                a(f'if [ -z "$runs" ]; then echo "  No runs found for {bp}"; else')
                a('    n=0')
                a('    for run in $runs; do')
                a('        n=$((n+1))')
                a('        if [ "$MAX_RUNS" -gt 0 ] && [ "$n" -gt "$MAX_RUNS" ]; then break; fi')
                a('        echo "  Downloading $run..."')
                a(f'        $DL $DL_ARGS --outdir "{dataset_id}/{bp}" "$run"')
                a(f'        [ "$DL" = "fasterq-dump" ] && gzip -f "{dataset_id}/{bp}/$run"*.fastq 2>/dev/null || true')
                a('    done')
                a('fi')
                a("")

        a('echo "Completed: $(date)"')
        a(f'echo "Estimated upper-bound size: ~{total_size} GB"')
        a('find . -name "*.fastq.gz" -type f | wc -l | xargs echo "Total compressed FASTQ files:"')
        a('du -sh */ 2>/dev/null | sort -h || true')

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

  # Preview what a dataset would download (no files written)
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
                        help='Specific BioProject to download (optional)')
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
                        help='Maximum number of parallel downloads (reserved)')
    parser.add_argument('--skip-validation', action='store_true',
                        help='Do NOT filter discovered runs by organism/assay '
                             '(unsafe: may download non-oyster or non-methylation data)')

    args = parser.parse_args()

    downloader = MethylationDataDownloader(
        output_dir=args.output_dir,
        max_parallel=args.max_parallel,
        validate=not args.skip_validation,
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
