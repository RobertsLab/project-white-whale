# DNA Methylation Datasets for Crassostrea / Magallana gigas

All BioProjects listed here were **verified live against NCBI SRA (2026-08)**:
each one actually contains *Crassostrea gigas* / *Magallana gigas* runs from a
methylation assay (Bisulfite-Seq or MeDIP-Seq). NCBI treats the two genus names
as a single taxon, so run metadata may show either name. Sample counts are the
real number of SRA runs and sizes are the summed SRA download sizes at the time
of verification; both can grow if submitters add runs.

These datasets are the ones the downloader in [`../code/`](../code/) uses. Each
BioProject appears in exactly one dataset (no duplicates), so the totals below
do not double-count.

## Verified Methylation Studies

### 1. POMS adaptation (WGBS)
- **Dataset ID**: `wgbs_poms_adaptation`
- **BioProject**: PRJEB60400
- **Title**: Genetic and epigenetic variations underpin rapid adaptation of oyster populations to Pacific Oyster Mortality Syndrome (POMS)
- **Method**: WGBS (Bisulfite-Seq)
- **Runs**: 246
- **Size**: ~450-500 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJEB60400

### 2. GESTINOV POMS gill/mantle (EM-seq)
- **Dataset ID**: `emseq_poms_gestinov`
- **BioProject**: PRJEB81880
- **Title**: Enzymatic Methyl sequencing of gills and mantle before and after POMS infection (GESTINOV 2021)
- **Method**: Enzymatic methyl-seq (labelled Bisulfite-Seq in SRA)
- **Tissues**: gill, mantle
- **Runs**: 40
- **Size**: ~400-450 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJEB81880

### 3. Aging / DECICOMP (WGBS)
- **Dataset ID**: `wgbs_aging_decicomp`
- **BioProject**: PRJEB105019
- **Title**: Epigenetic profiling (DNA methylation) of Pacific oysters under control conditions across three ages (4, 16, 28 months), DECICOMP project
- **Method**: WGBS (Bisulfite-Seq)
- **Runs**: 60
- **Size**: ~450-500 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJEB105019

### 4. PESTO methylseq (WGBS)
- **Dataset ID**: `wgbs_pesto`
- **BioProject**: PRJEB58545
- **Title**: Methylseq of Crassostrea gigas for PESTO project, 2022
- **Method**: WGBS (Bisulfite-Seq)
- **Runs**: 48
- **Size**: ~350-400 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJEB58545

### 5. Transgenerational infection resistance (WGBS)
- **Dataset ID**: `wgbs_transgen_infection`
- **BioProject**: PRJNA609264
- **Title**: Transgenerational adaptive phenotypic plasticity in infection resistance upon early microbial exposure
- **Method**: WGBS (Bisulfite-Seq)
- **Runs**: 47
- **Size**: ~600-650 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA609264

### 6. pH / ploidy WGBS
- **Dataset ID**: `wgbs_ph_ploidy`
- **BioProject**: PRJNA682817
- **Title**: Whole genome bisulfite sequencing (WGBS) of diploid and triploid ctenidia exposed to different pH levels
- **Method**: WGBS (Bisulfite-Seq)
- **Tissues**: gill (ctenidia)
- **Runs**: 24
- **Size**: ~100-110 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA682817

### 7. Developmental methylome (MeDIP-seq)
- **Dataset ID**: `medip_development`
- **BioProject**: PRJNA324546
- **Title**: Crassostrea gigas developmental genome-wide methylome dynamics
- **Method**: MeDIP-seq
- **Tissues**: embryo / larval developmental stages
- **Runs**: 21
- **Size**: ~5-10 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA324546

### 8. Epigenomics WGBS series
- **Dataset ID**: `wgbs_epigenomics_series`
- **BioProjects**: PRJNA807732, PRJNA562805, PRJNA213124
- **Title**: "Magallana gigas Epigenomics" (three related bisulfite projects)
- **Method**: WGBS (Bisulfite-Seq)
- **Runs**: 50 (24 + 12 + 14)
- **Size**: ~300-330 GB
- **URL**: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA807732

## Summary

| Metric | Value |
|--------|-------|
| Datasets | 8 |
| BioProjects | 10 |
| Total runs (samples) | ~536 |
| Total size | ~2.6-3.0 TB |
| Methods | WGBS, EM-seq, MeDIP-seq |

## Key Research Themes
- **Disease adaptation**: POMS resistance and infection response (PRJEB60400, PRJEB81880, PRJNA609264)
- **Environmental epigenetics**: pH / acidification response (PRJNA682817)
- **Developmental epigenetics**: methylome dynamics across development (PRJNA324546)
- **Aging**: methylation across age classes (PRJEB105019)

## Data Access Notes
- All datasets are public in NCBI SRA and require no authentication.
- The downloader re-checks every accession against NCBI at run time and keeps
  only oyster bisulfite/MeDIP runs, so a stale accession is reported rather than
  silently downloaded.
- WGBS/EM-seq datasets are large (hundreds of GB); MeDIP-seq is comparatively
  small.

## How These Were Verified
Discovered with NCBI E-utilities using:
```
"Crassostrea gigas"[Organism] AND "Bisulfite-Seq"[Strategy]
"Crassostrea gigas"[Organism] AND "MeDIP-Seq"[Strategy]
```
Runs were grouped by BioProject and each project's organism, library strategy,
run count and download size were read from the SRA `runinfo` records.

## Last Updated
August 2026 — rebuilt from verified NCBI SRA queries (replaces earlier
unverified accession list).
