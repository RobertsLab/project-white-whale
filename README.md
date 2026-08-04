# White Whale - Crassostrea gigas Public Datasets

This repository documents public RNA-seq and DNA methylation datasets for
*Crassostrea gigas* (Pacific oyster) and *Magallana gigas* (updated genus name),
and includes a downloader for the curated methylation BioProjects.

## At a Glance

| | Value |
|---|---|
| Verified methylation runs | ~536 (across 10 BioProjects) |
| Verified methylation volume | ~2.6–3.0 TB |
| Downloadable methylation datasets | 8 verified sets |
| RNA-seq accessions | Pending verification |

See [`dataset-summaries/consolidated-summary.md`](dataset-summaries/consolidated-summary.md)
for the full catalog.

## Directory Structure

- [`action-plan.md`](action-plan.md) — Systematic search and cataloging plan
- [`code/`](code/) — Scripts for downloading identified methylation datasets ([`HOW-TO.md`](code/HOW-TO.md))
- [`ncbi-datasets/`](ncbi-datasets/) — Datasets found in NCBI (SRA, GEO, BioProject)
- [`literature-review/`](literature-review/) — Key publications with associated data
- [`dataset-summaries/`](dataset-summaries/) — Consolidated inventory and metrics
- [`file-size-estimates/`](file-size-estimates/) — Size estimation methodology

## Research Scope

Focus: publicly available datasets with:

- RNA-seq / transcriptome data
- DNA methylation data (WGBS, RRBS, MeDIP-seq, targeted bisulfite, etc.)

For each dataset, the documentation aims to capture sample counts, tissue types,
environmental conditions, estimated file size, and access information (SRA,
GEO, BioProject).

## Search Terms Used

### Species Names
- *Crassostrea gigas*
- *Magallana gigas* (updated genus classification)

### Data Types
- RNA-seq / Transcriptome
- Bisulfite sequencing / DNA methylation
- MeDIP-seq
- WGBS (Whole Genome Bisulfite Sequencing)
- RRBS (Reduced Representation Bisulfite Sequencing)

## Downloading Data

The `code/` directory downloads DNA methylation datasets from NCBI SRA.
Sequencing runs are discovered live from NCBI (no accessions hardcoded) and each
run is validated as a *Crassostrea/Magallana gigas* bisulfite/MeDIP run before
download, so you only ever get on-target files.

```bash
cd code
./install_dependencies.sh                                       # one-time: install SRA Toolkit
python3 download_methylation_data.py --list                     # see available datasets
python3 download_methylation_data.py --dataset wgbs_ph_ploidy --dry-run  # preview safely
python3 download_methylation_data.py --dataset wgbs_ph_ploidy   # download
```

**Dataset IDs:** `wgbs_poms_adaptation`, `emseq_poms_gestinov`, `wgbs_aging_decicomp`,
`wgbs_pesto`, `wgbs_transgen_infection`, `wgbs_ph_ploidy`, `medip_development`,
`wgbs_epigenomics_series`

Start with the step-by-step [`code/HOW-TO.md`](code/HOW-TO.md).
For options and troubleshooting, see [`code/USAGE.md`](code/USAGE.md).

## Last Updated
- Created: December 2024
- Updated: August 2026 — methylation catalog rebuilt from verified NCBI SRA
  queries; downloader now validates organism/assay per run
