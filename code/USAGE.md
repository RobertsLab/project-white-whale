# DNA Methylation Data Downloader - Usage Guide

## Overview

This script downloads DNA methylation datasets for *Crassostrea gigas* and *Magallana gigas* from NCBI SRA based on the datasets documented in this repository.

## Installation

### Prerequisites

1. Python 3.6 or higher
2. SRA Toolkit
3. NCBI Entrez Direct tools (optional, for enhanced functionality)

### Quick Installation

Run the installation script:
```bash
./install_dependencies.sh
```

### Manual Installation

#### SRA Toolkit
- **Conda**: `conda install -c bioconda sra-tools`
- **Ubuntu/Debian**: `sudo apt-get install sra-toolkit`
- **Manual**: Download from https://github.com/ncbi/sra-tools

#### NCBI Entrez Direct (optional)
```bash
curl -s https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh | bash
export PATH=${PATH}:${HOME}/edirect
```

## Usage Examples

### List Available Datasets
```bash
python download_methylation_data.py --list
```

### Download a Complete Dataset
```bash
# Download the POMS adaptation WGBS study
python download_methylation_data.py --dataset wgbs_poms_adaptation

# Download with custom output directory
python download_methylation_data.py --dataset wgbs_pesto --output-dir /data/methylation
```

### Download Specific BioProject
```bash
# Download only PRJNA682817 from the pH/ploidy WGBS dataset
python download_methylation_data.py --dataset wgbs_ph_ploidy --bioproject PRJNA682817
```

### Limited Downloads for Testing
```bash
# Download only first 5 runs for testing
python download_methylation_data.py --dataset medip_development --max-runs 5 --dry-run

# Remove --dry-run to actually download
python download_methylation_data.py --dataset medip_development --max-runs 5
```

### Generate Download Scripts
```bash
# Create shell script for multiple datasets
python download_methylation_data.py --create-script --datasets wgbs_ph_ploidy medip_development wgbs_pesto

# Execute the generated script
cd methylation_data
chmod +x download_script.sh
./download_script.sh
```

## Available Datasets

All BioProjects below were verified live against NCBI SRA (2026-08): each
contains *Crassostrea/Magallana gigas* Bisulfite-Seq or MeDIP-Seq runs. Every
BioProject appears in exactly one dataset (no duplicates). "Samples" is the
real run count and "Size" is the summed SRA download size; both may grow if
submitters add runs.

| Dataset ID | Method | BioProjects | Size (GB) | Samples |
|------------|--------|-------------|-----------|---------|
| `wgbs_poms_adaptation` | WGBS | PRJEB60400 | 450-500 | 246 |
| `emseq_poms_gestinov` | EM-seq (Bisulfite-Seq) | PRJEB81880 | 400-450 | 40 |
| `wgbs_aging_decicomp` | WGBS | PRJEB105019 | 450-500 | 60 |
| `wgbs_pesto` | WGBS | PRJEB58545 | 350-400 | 48 |
| `wgbs_transgen_infection` | WGBS | PRJNA609264 | 600-650 | 47 |
| `wgbs_ph_ploidy` | WGBS | PRJNA682817 | 100-110 | 24 |
| `medip_development` | MeDIP-seq | PRJNA324546 | 5-10 | 21 |
| `wgbs_epigenomics_series` | WGBS | PRJNA807732, PRJNA562805, PRJNA213124 | 300-330 | 50 |

**Totals:** 8 datasets, 10 BioProjects, ~536 runs, ~2.6-3.0 TB.

## Output Structure

```
methylation_data/
├── download_script.sh          # Generated download script (if --create-script used)
├── <dataset_id>/
│   ├── dataset_info.json       # Dataset metadata
│   ├── <bioproject>/
│   │   ├── runs.txt            # Run accessions discovered from NCBI
│   │   ├── <SRR_accession>_1.fastq.gz
│   │   ├── <SRR_accession>_2.fastq.gz
│   │   └── ...
│   └── ...
└── download_methylation_data.log  # Download progress log
```

## NCBI Rate Limits (optional but recommended)

Run discovery queries NCBI. To raise the request limit, set these before running:

```bash
export NCBI_EMAIL="you@example.com"     # identifies you to NCBI
export NCBI_API_KEY="your_ncbi_api_key" # optional, from your NCBI account
```

## Advanced Options

### Parallel Downloads
```bash
# Use 4 parallel download processes
python download_methylation_data.py --dataset wgbs_poms_adaptation --max-parallel 4
```

### Custom Selection
```bash
# Download specific runs (requires modification of script)
# See source code for adding custom run lists
```

### Run Validation (on by default)

Discovered runs are filtered against their NCBI metadata: only
*Crassostrea/Magallana gigas* runs from a methylation assay (Bisulfite-Seq or
MeDIP-Seq) are kept. Off-target runs are logged as `Skipped N off-target
run(s)` and never downloaded. This is what prevents a wrong or stale BioProject
accession from pulling non-oyster or non-methylation data.

```bash
# Disable the filter (NOT recommended): download whatever a BioProject contains
python download_methylation_data.py --dataset wgbs_ph_ploidy --skip-validation
```

## Storage Requirements

- **Minimum**: <10 GB for the smallest dataset (`medip_development`)
- **Typical**: 100-500 GB per WGBS dataset
- **Largest single dataset**: ~650 GB (`wgbs_transgen_infection`)
- **Full collection (all 8 datasets)**: ~2.6-3.0 TB estimated

## Download Speed Estimates

Based on typical academic internet connections:
- **100 Mbps**: ~45 GB/hour → 4-5 hours for 200 GB dataset
- **1 Gbps**: ~450 GB/hour → ~30 minutes for 200 GB dataset
- **Residential**: Highly variable, plan for overnight downloads

## Troubleshooting

### SRA Toolkit Issues
```bash
# Configure SRA toolkit
vdb-config --interactive

# Test SRA toolkit
fastq-dump --version
fasterq-dump --version
```

### Network Issues
- Large datasets may take hours to days to download
- Use `--max-runs` to test with smaller subsets
- Consider institutional high-speed networks
- Resume interrupted downloads (SRA toolkit handles this automatically)

### Storage Issues
- Monitor disk space during downloads
- Consider downloading to external drives for large datasets
- Use compression tools if needed

### Authentication Issues
- Most datasets are public and require no authentication
- Some may require NCBI account for controlled access

## Quality Control

After download:
1. Check file integrity with MD5 sums (if provided)
2. Verify file sizes match expectations
3. Test random samples with quality assessment tools
4. Check for complete paired-end files

## Data Citation

When using these datasets, cite:
1. The original publications (see literature-review/ directory)
2. NCBI SRA: "Data were obtained from the NCBI Sequence Read Archive (SRA)"
3. Specific BioProject accessions used

## Support

For issues with:
- **This script**: Create issue in project repository
- **SRA Toolkit**: https://github.com/ncbi/sra-tools/issues
- **Data access**: Contact original data submitters or NCBI help

## License

This script is provided as-is for research purposes. Respect the original data licenses and publication requirements.