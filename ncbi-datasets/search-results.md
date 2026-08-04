# NCBI Database Search Results

## Search Strategy

Searched NCBI databases using the following approach:
1. **NCBI SRA** (Sequence Read Archive) - for raw sequencing data - https://www.ncbi.nlm.nih.gov/sra/
2. **NCBI GEO** (Gene Expression Omnibus) - for processed datasets - https://www.ncbi.nlm.nih.gov/geo/
3. **NCBI BioProject** - for project-level information - https://www.ncbi.nlm.nih.gov/bioproject/

## Search Terms Used

### Primary searches:
- "Crassostrea gigas" AND "RNA-seq" - https://www.ncbi.nlm.nih.gov/sra/?term=Crassostrea+gigas+RNA-seq
- "Crassostrea gigas" AND "transcriptome" - https://www.ncbi.nlm.nih.gov/sra/?term=Crassostrea+gigas+transcriptome
- "Crassostrea gigas" AND "bisulfite" - https://www.ncbi.nlm.nih.gov/sra/?term=Crassostrea+gigas+bisulfite
- "Crassostrea gigas" AND "methylation" - https://www.ncbi.nlm.nih.gov/sra/?term=Crassostrea+gigas+methylation
- "Magallana gigas" AND "RNA-seq" - https://www.ncbi.nlm.nih.gov/sra/?term=Magallana+gigas+RNA-seq
- "Magallana gigas" AND "transcriptome" - https://www.ncbi.nlm.nih.gov/sra/?term=Magallana+gigas+transcriptome
- "Magallana gigas" AND "bisulfite" - https://www.ncbi.nlm.nih.gov/sra/?term=Magallana+gigas+bisulfite
- "Magallana gigas" AND "methylation" - https://www.ncbi.nlm.nih.gov/sra/?term=Magallana+gigas+methylation

## Datasets Identified

*Note: This section will be populated with specific dataset information*

### RNA-seq Datasets

> ⚠️ **Unverified.** The RNA-seq BioProject accessions in
> [`rna-seq-datasets.md`](rna-seq-datasets.md) have **not** yet been confirmed
> against NCBI. A spot check found several that are empty or point to the wrong
> organism, so they must be re-verified before use. The downloader in
> [`../code/`](../code/) does not rely on them.

### DNA Methylation Datasets (verified 2026-08)

**BioProjects confirmed to contain *Crassostrea/Magallana gigas* bisulfite/MeDIP runs:**
- PRJEB60400: POMS adaptation WGBS (246 runs)
- PRJEB81880: GESTINOV POMS gill/mantle EM-seq (40 runs)
- PRJEB105019: Aging / DECICOMP WGBS (60 runs)
- PRJEB58545: PESTO methylseq WGBS (48 runs)
- PRJNA609264: Transgenerational infection resistance WGBS (47 runs)
- PRJNA682817: pH / ploidy WGBS (24 runs)
- PRJNA324546: Developmental methylome MeDIP-seq (21 runs)
- PRJNA807732, PRJNA562805, PRJNA213124: "Magallana gigas Epigenomics" WGBS series (50 runs)

See [`dna-methylation-datasets.md`](dna-methylation-datasets.md) for full details.

## Summary Statistics

- **DNA methylation BioProjects verified**: 10 (across 8 datasets)
- **DNA methylation samples (runs)**: ~536
- **DNA methylation total size**: ~2.6-3.0 TB
- **RNA-seq BioProjects**: pending re-verification (see warning above)

## Last Updated
December 2024