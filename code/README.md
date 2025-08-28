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
   python download_methylation_data.py --dataset targeted_bisulfite --max-runs 3 --dry-run
   ```

4. **Generate download scripts**:
   ```bash
   python download_methylation_data.py --create-script --datasets medip_seq targeted_bisulfite
   ```

## Available Datasets

The script can download data from 7 major DNA methylation study categories:

- **WGBS Studies**: Whole genome bisulfite sequencing (200-400 GB)
- **RRBS Studies**: Reduced representation bisulfite sequencing (40-100 GB)
- **MeDIP-seq Studies**: Methylated DNA immunoprecipitation sequencing (30-60 GB)
- **Targeted Bisulfite**: Gene-specific methylation analysis (10-30 GB)
- **Magallana Studies**: Recent studies with updated nomenclature (100-200 GB)

Total estimated size for all datasets: **580-1,220 GB**

## Requirements

- Python 3.6+
- SRA Toolkit (fastq-dump/fasterq-dump)
- NCBI Entrez Direct tools (optional, for enhanced functionality)
- Sufficient storage space (see individual dataset sizes)
- Stable internet connection

## Documentation

- `USAGE.md` - Comprehensive usage guide with examples
- Repository documentation in `../ncbi-datasets/dna-methylation-datasets.md`

## Support

For detailed usage instructions and troubleshooting, see `USAGE.md`.

## Last Updated
December 2024