# PubMed Literature Mining Results (Phase 2)

This document records the **executed** Phase 2 literature-mining searches from
[`action-plan.md`](../action-plan.md). All searches were run live against the
NCBI PubMed database using the E-utilities API. Result counts and query
translations are reproducible with the query URLs provided.

- **Search date**: 2026-08-04
- **Tool**: NCBI E-utilities (`esearch`, `esummary`, `elink`)
- **Databases queried**: PubMed, BioProject, SRA, GEO DataSets

> Note on syntax: the draft queries in `action-plan.md` used `[Organism]`, which
> is **not** a valid PubMed field (it is an SRA/nucleotide field). PubMed treats
> species names as text, so the executed queries below use free-text / `[Author]`
> fields instead. The `[Organism]` form is still correct for the Phase 1 SRA
> searches.

## 2.1 PubMed Search Strategy — queries and result counts

Species clause used below: `("Crassostrea gigas" OR "Magallana gigas")`

| # | Query | Hits |
|---|-------|------|
| 1 | species AND ("RNA-seq" OR "transcriptome" OR "transcriptomic") | 248 |
| 2 | species AND ("methylation" OR "bisulfite" OR "epigenetic") | 54 |
| 3 | species AND ("genome" OR "genomic") AND ("sequencing") | 77 |
| 4 | species AND "RNA-seq" | 50 |
| 5 | species AND ("bisulfite" OR "WGBS" OR "RRBS" OR "MeDIP") | 14 |
| 6 | "Magallana gigas" AND ("transcriptome" OR "RNA-seq" OR "methylation" OR "bisulfite") | 11 |
| 7 | species AND ("heat stress" OR "temperature" OR "acidification" OR "salinity" OR "hypoxia") | 417 |

Reproduce any count by pasting the query into the PubMed search box, or via the
E-utilities endpoint, e.g. the methylation query:

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=("Crassostrea gigas" OR "Magallana gigas") AND ("methylation" OR "bisulfite" OR "epigenetic")
```

The methylation set (54 records) was retrieved in full; the RNA-seq/transcriptome
set was retrieved as the 30 most recent records (of 248). See
[`key-publications.md`](key-publications.md) for the curated, verified entries.

## 2.2 Journal-Specific Distribution

Journal counts across the retrieved methylation, recent RNA-seq, and Roberts-lab
record sets (deduplicated within each set). This shows where oyster genomics/
epigenetics data-generating papers concentrate — useful for targeted browsing.

| Count | Journal |
|-------|---------|
| 13 | Comp Biochem Physiol Part D Genomics Proteomics |
| 9 | BMC Genomics |
| 8 | Mar Biotechnol (NY) |
| 4 | Genes (Basel) |
| 4 | Fish Shellfish Immunol |
| 4 | Front Physiol |
| 3 | Front Immunol |
| 3 | Glob Chang Biol |
| 3 | PeerJ |
| 3 | Mar Environ Res |
| 2 | Genomics / Sci Total Environ / Int J Biol Macromol / Gene / PLoS One / Environ Epigenet |

Takeaway: *BMC Genomics*, *Marine Biotechnology*, and *Comp Biochem Physiol Part
D (Genomics & Proteomics)* are the highest-yield journals for datasets; immunity
data cluster in *Fish & Shellfish Immunology* / *Frontiers in Immunology*, and
environmental-epigenetics data in *Sci Total Environ* / *Environmental
Epigenetics* / *Marine Environmental Research*.

## 2.3 Author-Based Searches

Query form: `<Author>[Author] AND (Crassostrea OR Magallana OR oyster)`

| Hits | Author | Group / focus |
|------|--------|----------------|
| 24 | Roberts SB | Univ. Washington — oyster DNA methylation / functional genomics |
| 50 | Favrel P | Univ. Caen (France) — oyster development / epigenetics |
| 22 | Riviere G | Univ. Caen (France) — DNA methylation & development |
| 159 | Li L | IOCAS (China) — oyster genome / population genomics (broad name, needs disambiguation) |
| 7 | Grunau C | Univ. Perpignan (France) — environmental epigenetics |

The Roberts SB set (24 records) was retrieved in full and cross-checked; note it
includes closely related species studied by the same lab (*Crassostrea
virginica*, *Ostrea lurida*, *Panopea generosa*) that share methods and are
relevant methodological references. "Li L" is a common name and returns
unrelated authors — disambiguate by affiliation (IOCAS) before use.

## Dataset Linkage (PubMed → sequence archives)

Each retrieved PMID was cross-linked to BioProject, SRA, and GEO via `elink`.

**Key finding:** NCBI maintains PubMed↔archive links for only a small fraction of
these papers. Of the featured methylation studies, only **Wang et al. 2014**
(PMID 25514978) carries live links. For most papers, accession numbers are stated
only in the article's *Data Availability* section and must be verified
per-paper (a task folded into Phase 4 cataloging/verification).

**Verified oyster BioProjects (resolved and confirmed as _Magallana gigas_):**

| Accession | Title | SRA experiments | Source |
|-----------|-------|-----------------|--------|
| PRJNA173440 | Genome-wide and single-base resolution DNA methylomes of the Pacific oyster | 4 | linked from PMID 25514978 |
| PRJNA146329 | RNA-seq of Pacific oyster under different developmental stages | 112 | linked from PMID 25514978 |

### Data-integrity correction

The BioProject accessions previously listed in
[`../dataset-summaries/consolidated-summary.md`](../dataset-summaries/consolidated-summary.md)
and [`../ncbi-datasets/search-results.md`](../ncbi-datasets/search-results.md)
were checked against NCBI and **do not correspond to oyster data**. They resolve
to unrelated organisms or do not exist:

| Previously claimed | Actually resolves to |
|--------------------|----------------------|
| PRJNA85067 ("developmental transcriptomes") | *Saccharomyces cerevisiae* (yeast) |
| PRJNA316216 ("Roberts WGBS") | *Panopea generosa* (geoduck) |
| PRJNA394801 ("OA methylation") | *Sorghum bicolor* |
| PRJNA506631 ("salinity stress") | *Erwinia amylovora* |
| PRJNA248740 ("temperature stress") | *Arabidopsis thaliana* |
| PRJNA422851 ("pathogen response") | *Triuncina daii* |
| PRJNA486983 ("reproductive dev.") | not found |
| PRJNA725689 ("recent Magallana") | not found |

These should be treated as **unverified/incorrect** and replaced with
per-paper–verified accessions during Phase 4. Do not cite them as oyster
datasets.

## Last Updated
2026-08-04 — Phase 2 executed (searches run live against NCBI PubMed/BioProject/SRA/GEO)
