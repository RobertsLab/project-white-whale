# White Whale - Crassostrea gigas Public Datasets

This directory contains research and documentation of public datasets with RNA-seq and DNA methylation data for *Crassostrea gigas* (Pacific oyster) and *Magallana gigas* (updated genus name).

## Directory Structure

- `code/` - Scripts for downloading the identified datasets (see [`code/HOW-TO.md`](code/HOW-TO.md))
- `ncbi-datasets/` - Documentation of datasets found in NCBI databases
- `literature-review/` - Information from published studies
- `dataset-summaries/` - Consolidated summaries of identified datasets
- `file-size-estimates/` - Estimates of cumulative file sizes for datasets

## Research Scope

This research focuses on identifying publicly available datasets containing:
- RNA-seq data
- DNA methylation data (including bisulfite sequencing, MeDIP-seq, etc.)

For each dataset, we aim to collect:
- Number of samples
- Tissue types analyzed
- Environmental conditions
- Estimated cumulative file size
- Access information (SRA accessions, GEO series, etc.)

## Search Terms Used

### Species Names
- *Crassostrea gigas*
- *Magallana gigas* (updated genus classification)

### Data Types
- RNA-seq
- Transcriptome
- Bisulfite sequencing
- DNA methylation
- MeDIP-seq
- WGBS (Whole Genome Bisulfite Sequencing)
- RRBS (Reduced Representation Bisulfite Sequencing)

## Downloading Data

The `code/` directory contains a downloader that fetches the identified DNA
methylation datasets from NCBI SRA. Sequencing runs are discovered live from
NCBI (no accession numbers are hardcoded), so you always get the real files.

```bash
cd code
./install_dependencies.sh                                   # one-time: install SRA Toolkit
python3 download_methylation_data.py --list                 # see available datasets
python3 download_methylation_data.py --dataset wgbs_roberts --dry-run  # preview safely
python3 download_methylation_data.py --dataset wgbs_roberts # download
```

New to this? Start with the step-by-step [`code/HOW-TO.md`](code/HOW-TO.md).
For all options and troubleshooting, see [`code/USAGE.md`](code/USAGE.md).

## Last Updated
- Created: December 2024
- Updated: June 2026 — added data downloader (`code/`)