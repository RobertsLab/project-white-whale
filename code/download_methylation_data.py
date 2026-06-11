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

# Dataset information extracted from repository documentation.
#
# NOTE: the BioProject accessions below were compiled from literature and
# should be treated as a starting point. The script verifies each one against
# NCBI at run time -- if an accession is wrong or has no runs, the script tells
# you rather than guessing.
METHYLATION_DATASETS = {
    "wgbs_roberts": {
        "description": "Roberts Lab WGBS Studies",
        "bioprojects": ["PRJNA316216", "PRJNA394801"],
        "method": "WGBS",
        "estimated_samples": "30-50",
        "tissue_types": ["gonad", "gill", "mantle", "digestive_gland"],
        "estimated_size_gb": "200-400",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=roberts+crassostrea+bisulfite",
        "notes": "High-quality WGBS from University of Washington Roberts Lab"
    },
    "wgbs_ocean_acidification": {
        "description": "Ocean Acidification Methylation Study",
        "bioprojects": ["PRJNA394801", "PRJNA316216"],
        "method": "WGBS",
        "estimated_samples": "20-30",
        "tissue_types": ["gill", "mantle"],
        "estimated_size_gb": "150-250",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=crassostrea+pH+methylation",
        "notes": "DNA methylation response to ocean acidification"
    },
    "rrbs_developmental": {
        "description": "Developmental Methylation Studies",
        "bioprojects": ["PRJNA486983", "PRJNA273482"],
        "method": "RRBS",
        "estimated_samples": "25-35",
        "tissue_types": ["gonad", "larvae", "spat"],
        "estimated_size_gb": "50-100",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=crassostrea+RRBS",
        "notes": "RRBS during oyster development and reproduction"
    },
    "rrbs_environmental_stress": {
        "description": "Environmental Stress RRBS",
        "bioprojects": ["PRJNA506631", "PRJNA413624"],
        "method": "RRBS",
        "estimated_samples": "20-30",
        "tissue_types": ["various_adult_tissues"],
        "estimated_size_gb": "40-80",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=crassostrea+stress+methylation",
        "notes": "Methylation changes under environmental stress"
    },
    "medip_seq": {
        "description": "Genome-wide Methylation Profiling",
        "bioprojects": ["PRJNA348937", "PRJNA394425"],
        "method": "MeDIP-seq",
        "estimated_samples": "15-25",
        "tissue_types": ["adult_tissues"],
        "estimated_size_gb": "30-60",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=crassostrea+MeDIP",
        "notes": "MeDIP-seq for genome-wide methylation patterns"
    },
    "targeted_bisulfite": {
        "description": "Gene-specific Methylation Studies",
        "bioprojects": ["PRJNA311096", "PRJNA381456"],
        "method": "Targeted Bisulfite",
        "estimated_samples": "20-40",
        "tissue_types": ["multiple_tissue_types"],
        "estimated_size_gb": "10-30",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=crassostrea+targeted+bisulfite",
        "notes": "Targeted analysis of specific gene regions"
    },
    "magallana_recent": {
        "description": "Recent Magallana gigas Studies",
        "bioprojects": ["PRJNA725689", "PRJNA688412"],
        "method": "Mixed methods",
        "estimated_samples": "15-25",
        "tissue_types": ["gonad", "gill", "mantle"],
        "estimated_size_gb": "100-200",
        "search_url": "https://www.ncbi.nlm.nih.gov/sra/?term=magallana+methylation",
        "notes": "Studies using updated Magallana gigas nomenclature"
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


def _parse_runs_from_runinfo(runinfo_csv: str) -> List[str]:
    """Extract run accessions from an SRA runinfo CSV using the 'Run' column."""
    runs: List[str] = []
    reader = csv.DictReader(io.StringIO(runinfo_csv))
    if not reader.fieldnames or "Run" not in reader.fieldnames:
        return runs
    for row in reader:
        acc = (row.get("Run") or "").strip()
        if acc.upper().startswith(RUN_ACCESSION_PREFIXES):
            runs.append(acc)
    return runs


def fetch_runs_via_eutils(bioproject: str) -> Optional[List[str]]:
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

    runs = _parse_runs_from_runinfo(runinfo)
    logger.info(f"Found {len(runs)} run(s) in {bioproject}")
    return runs


def fetch_runs_via_edirect(bioproject: str) -> Optional[List[str]]:
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
    return _parse_runs_from_runinfo(result.stdout)


class MethylationDataDownloader:
    """Main class for downloading DNA methylation datasets."""

    def __init__(self, output_dir: str = "./methylation_data", max_parallel: int = 2):
        self.output_dir = Path(output_dir)
        self.max_parallel = max_parallel
        self.output_dir.mkdir(parents=True, exist_ok=True)

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
        runs = fetch_runs_via_eutils(bioproject)
        if runs is None:
            logger.warning(f"HTTPS discovery unavailable for {bioproject}; trying Entrez Direct")
            runs = fetch_runs_via_edirect(bioproject)

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
        a("# Print the SRR/ERR/DRR run accessions for a BioProject, one per line.")
        a("get_runs() {")
        a("    local bp=\"$1\"")
        a('    local base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"')
        a('    local hist; hist=$(curl -s "${base}/esearch.fcgi?db=sra&term=${bp}[BioProject]&usehistory=y&retmax=1")')
        a("    local webenv querykey")
        a("    webenv=$(echo \"$hist\" | sed -n 's:.*<WebEnv>\\(.*\\)</WebEnv>.*:\\1:p')")
        a("    querykey=$(echo \"$hist\" | sed -n 's:.*<QueryKey>\\(.*\\)</QueryKey>.*:\\1:p')")
        a('    [ -z "$webenv" ] && return 0')
        a('    curl -s "${base}/efetch.fcgi?db=sra&WebEnv=${webenv}&query_key=${querykey}&rettype=runinfo&retmode=text" \\')
        a("        | awk -F, 'NR>1 && $1 ~ /^[SED]RR/ {print $1}'")
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
  python download_methylation_data.py --dataset wgbs_roberts --dry-run

  # Download a specific dataset
  python download_methylation_data.py --dataset wgbs_roberts

  # Download a few runs of one BioProject to test your setup
  python download_methylation_data.py --dataset rrbs_developmental \\
      --bioproject PRJNA486983 --max-runs 3

  # Create a standalone shell script for several datasets
  python download_methylation_data.py --create-script \\
      --datasets wgbs_roberts rrbs_developmental

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

    args = parser.parse_args()

    downloader = MethylationDataDownloader(
        output_dir=args.output_dir,
        max_parallel=args.max_parallel
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
