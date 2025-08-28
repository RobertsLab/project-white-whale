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
# Download all WGBS studies from Roberts Lab
python download_methylation_data.py --dataset wgbs_roberts

# Download with custom output directory
python download_methylation_data.py --dataset rrbs_developmental --output-dir /data/methylation
```

### Download Specific BioProject
```bash
# Download only PRJNA486983 from developmental studies
python download_methylation_data.py --dataset rrbs_developmental --bioproject PRJNA486983
```

### Limited Downloads for Testing
```bash
# Download only first 5 runs for testing
python download_methylation_data.py --dataset medip_seq --max-runs 5 --dry-run

# Remove --dry-run to actually download
python download_methylation_data.py --dataset medip_seq --max-runs 5
```

### Generate Download Scripts
```bash
# Create shell script for multiple datasets
python download_methylation_data.py --create-script --datasets wgbs_roberts rrbs_developmental medip_seq

# Execute the generated script
cd methylation_data
chmod +x download_script.sh
./download_script.sh
```

## Available Datasets

| Dataset ID | Method | BioProjects | Size (GB) | Samples |
|------------|--------|-------------|-----------|---------|
| `wgbs_roberts` | WGBS | PRJNA316216, PRJNA394801 | 200-400 | 30-50 |
| `wgbs_ocean_acidification` | WGBS | PRJNA394801, PRJNA316216 | 150-250 | 20-30 |
| `rrbs_developmental` | RRBS | PRJNA486983, PRJNA273482 | 50-100 | 25-35 |
| `rrbs_environmental_stress` | RRBS | PRJNA506631, PRJNA413624 | 40-80 | 20-30 |
| `medip_seq` | MeDIP-seq | PRJNA348937, PRJNA394425 | 30-60 | 15-25 |
| `targeted_bisulfite` | Targeted | PRJNA311096, PRJNA381456 | 10-30 | 20-40 |
| `magallana_recent` | Mixed | PRJNA725689, PRJNA688412 | 100-200 | 15-25 |

## Output Structure

```
methylation_data/
├── download_script.sh          # Generated download script
├── <dataset_id>/
│   ├── dataset_info.json       # Dataset metadata
│   ├── <bioproject>/
│   │   ├── runinfo.csv         # Run information
│   │   ├── <SRR_accession>_1.fastq.gz
│   │   ├── <SRR_accession>_2.fastq.gz
│   │   └── ...
│   └── ...
└── download_log.txt            # Download progress log
```

## Advanced Options

### Parallel Downloads
```bash
# Use 4 parallel download processes
python download_methylation_data.py --dataset wgbs_roberts --max-parallel 4
```

### Custom Selection
```bash
# Download specific runs (requires modification of script)
# See source code for adding custom run lists
```

## Storage Requirements

- **Minimum**: ~100 GB for small datasets (targeted_bisulfite, medip_seq)
- **Typical**: ~500 GB for moderate datasets (RRBS studies)
- **Large**: ~2-4 TB for complete WGBS collections
- **Full collection**: ~6-10 TB estimated

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