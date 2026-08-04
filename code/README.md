# Code Directory

This directory contains scripts for downloading and processing the DNA methylation datasets identified in this repository.

## Scripts

- `download_methylation_data.py` - Main script for downloading DNA methylation datasets from NCBI SRA
- `install_dependencies.sh` - Installation script for required dependencies (SRA Toolkit, etc.)
- `requirements.txt` - List of dependencies and tools needed

## Quick Start

1. **Install dependencies**:
   ```bash
   ./install_dependencies.sh
   ```

2. **List available datasets**:
   ```bash
   python download_methylation_data.py --list
   ```

3. **Download a small dataset for testing**:
   ```bash
   python download_methylation_data.py --dataset medip_development --max-runs 3 --dry-run
   ```

4. **Generate download scripts**:
   ```bash
   python download_methylation_data.py --create-script --datasets medip_development wgbs_ph_ploidy
   ```

## Available Datasets

The script downloads from 8 verified DNA methylation datasets (10 BioProjects,
all confirmed to contain *Crassostrea/Magallana gigas* Bisulfite-Seq/MeDIP-Seq
runs against NCBI SRA). Highlights:

- **WGBS studies** (POMS adaptation, aging/DECICOMP, PESTO, transgenerational
  infection, pH/ploidy, epigenomics series): ~100-650 GB each
- **EM-seq** (GESTINOV POMS gill/mantle): ~400-450 GB
- **MeDIP-seq** (developmental methylome dynamics): ~5-10 GB

See `USAGE.md` for the full per-dataset table. Total estimated size for all
datasets: **~2.6-3.0 TB** (~536 runs).

## Requirements

- Python 3.6+
- SRA Toolkit (fastq-dump/fasterq-dump)
- NCBI Entrez Direct tools (optional, for enhanced functionality)
- Sufficient storage space (see individual dataset sizes)
- Stable internet connection

## Documentation

- `HOW-TO.md` - Beginner-friendly, step-by-step walkthrough (start here)
- `USAGE.md` - Comprehensive usage guide with examples
- Repository documentation in `../ncbi-datasets/dna-methylation-datasets.md`

## How run discovery works

The downloader never hardcodes sequence accessions. For each BioProject it asks
NCBI (via E-utilities over HTTPS — no extra software required) which sequencing
runs belong to it. Each run is then validated against its NCBI metadata and
only *Crassostrea/Magallana gigas* bisulfite/MeDIP runs are kept; off-target
runs are logged and skipped so a wrong accession is reported rather than
downloaded. If NCBI can't be reached, it reports the problem and stops rather
than downloading placeholder data.

## Support

For detailed usage instructions and troubleshooting, see `USAGE.md`.

## Last Updated
December 2024