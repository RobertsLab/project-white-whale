# Enrichment tier: MeDIP and MBD-seq

MeDIP-seq (an anti-5mC antibody) and MBD-seq (a methyl-CpG binding protein)
pull down methylated DNA fragments and sequence them. They tell you **where
methylated fragments are enriched**, not what fraction of a CpG's copies is
methylated. The assay_chemistry vocabulary gives both `tier: enrichment`, and
PLAN.md §4.1 says such data is processed and shown but never enters a
base-resolution comparison. This page describes how White Whale processes
them (`workflows/enrichment`), what the numbers mean, and where they show up.

Studies in this tier now: RIVIERE2017_DEVELOPMENT_MEDIP (PRJNA324546, 21
runs) and the Gavery 2013 gill MBD-seq set (probably inside PRJNA213124, see
`registry/backlog.tsv`).

## The workflow

```
FASTQC -> TRIM_GALORE -> BOWTIE2_ALIGN -> MARKDUP -> ENRICHMENT_COUNT -> ENRICHMENT_QC_GATE -> MANIFEST
```

| Step | What it does | Container |
|---|---|---|
| `FASTQC`, `TRIM_GALORE` | the methylation workflow's modules; plain adapter/quality trimming (no bisulfite preset) | as there |
| `BOWTIE2_BUILD` | bowtie2 index of the reference bundle's `genomic.fna`, kept as `bowtie2_genome/` in the bundle (or `$WHALE_BOWTIE2_INDEX`), with the FASTA's SHA-256 so a changed genome is noticed | `mulled-v2-ac74a7…:f70b31a2…-0` (bowtie2 2.5.2, samtools 1.18) |
| `BOWTIE2_ALIGN` | end-to-end alignment (`--sensitive`, plus `params.bowtie2_args`), collated, mate-fixed and sorted in one stream | same |
| `MARKDUP` | `samtools markdup -r` (duplicates removed, as `deduplicate: true` asks); then the fragment table: primary, mapped, non-duplicate records with MAPQ ≥ `params.min_mapq` (10), proper pairs only for paired libraries | same |
| `ENRICHMENT_COUNT` | `workflows/enrichment/bin/enrichment_counts.py count`: tiles, genes, QC (below) | `python:3.12.12` (standard library only) |
| `ENRICHMENT_QC_GATE` | `enrichment_counts.py qc-table` → `qc_summary.tsv` | same |
| `MANIFEST` | the methylation workflow's module: parameters (including `tier: enrichment` and `enrichment_method`), tool versions, input checksums | same |

Profiles, `WHALE_*` variables and study configs are the same as for the
methylation workflow (`config/base.config`; docs/compute_profiles.md). A
study config sets `chemistry = 'medip'` or `'mbd'`, and `whale fetch-fastq
--launch` then starts `workflows/enrichment` instead of `workflows/methylation`:

```bash
nextflow run workflows/enrichment -c config/RIVIERE2017_DEVELOPMENT_MEDIP.config -profile desktop,docker
```

Test it on the fixtures (two simulated MeDIP libraries on the 30 kb TEST
genome, `tests/fixtures/make_enrichment_fixtures.py`):

```bash
make enrichment-stub         # the DAG with -stub
make enrichment-real-run     # every tool in Docker, then tests/check_enrichment_run.py
```

A `-stub` run leaves an empty bowtie2 index in `results/test/bowtie2_index`
(storeDir); a real run refuses it, so delete it (the make target does) before
a real run.

## How enrichment is measured (and why this way)

### Current practice

- **MEDIPS** (Lienhard et al. 2014, Bioconductor) counts reads in genome-wide
  windows, reports counts / RPKM, and offers `rms`: signal divided by a CpG
  coupling factor from a linear fit of signal against local CpG density. Its
  QC — CpG enrichment (relH, GoGe), coverage of the genome's CpGs, and
  saturation from two random halves of a library — is the de facto standard.
- **QSEA** (Lienhard et al. 2017) replaces the linear fit with a sigmoidal
  CpG-density curve in a Bayesian model and estimates an absolute,
  bisulfite-like methylation level per window, using TMM between libraries.
- BayMeth, MeSiC and MeDEStrand also model the CpG-density bias to estimate
  absolute levels; all need a calibration (a fully methylated control, an
  input library, or a genome-wide curve shape assumed in advance).

### What White Whale does

**CpG-density-normalised window counts** (`enrichment_method.name:
cpg_density_normalised_window_counts`, version 1). For each 1 kb tile (the
tiles of PLAN.md §6.2), each gene body (`features/genes.bed`) and each
promoter (`features/promoters.bed`, intervals merged per gene):

| Column | Definition |
|---|---|
| `count` | deduplicated fragments whose midpoint lies in the region. Paired: the template; single-end: the read extended to `params.fragment_length` (200 bp, as MEDIPS' `extend`) |
| `cpm` | `count × 10⁶ / N`, N = fragments counted in the library (library-size normalisation) |
| `cpgs` | CpGs in the region (strand-merged, from the bundle's genome) |
| `obs_exp` | `(count / N) / (cpgs / C)`, C = CpGs in the genome: the region's reads per CpG over the genome's reads per CpG. 1 = as enriched as the average CpG; empty below `enrichment_min_cpgs` (5) CpGs |

Why not MEDIPS' `rms` or QSEA's absolute levels:

1. **Their calibrations assume a mammalian methylome.** Both fit the
   signal-vs-CpG-density curve over the whole genome, which works where most
   CpGs are methylated and CpG-poor DNA is the background. The oyster
   methylome is mosaic: roughly 15% of CpGs methylated, concentrated in gene bodies,
   with most intergenic DNA and many CpG-rich promoters unmethylated. A
   curve fitted on that confounds "few CpGs" with "unmethylated", so the
   absolute levels would look precise and be wrong. With no fully methylated
   control or input libraries in the public oyster sets, nothing here could
   check them.
2. **Both need a BSgenome package** of the assembly, which would have to be
   forged and pinned for GCF_963853765.1; a standard-library counter over the
   reference bundle's own FASTA and BEDs keeps one reference by construction
   (PLAN.md §2.1).
3. **`obs_exp` is the transparent part of the same idea.** It removes the two
   effects every method agrees on (library size; more CpGs, more capture)
   without asserting a curve shape, and keeps raw counts in the tables so a
   calibrated model (QSEA, or edgeR/DESeq2 on counts with CpG count as an
   offset) can be run downstream if a study ever has the controls for it.

What it gives up: `obs_exp` is not linear in methylation (capture saturates
in CpG-dense methylated regions and is near zero in CpG-poor ones), and a
region with few CpGs is noisy. Use it for ranks and fold differences, not
as a percentage.

### Per-library QC (`qc_summary.tsv`)

| Column | Meaning | Gate |
|---|---|---|
| `mapping_pct` | bowtie2 overall alignment rate | ≥ `enrichment_min_mapping` (50) |
| `duplication_pct` | `samtools markdup` duplicates / examined reads | reported |
| `fragments` | fragments counted (N) | ≥ `enrichment_min_fragments` (1 M) |
| `relh`, `goge` | MEDIPS.CpGenrich: CpG frequency (relH) and observed/expected CpG ratio (GoGe) in the counted fragments over the genome's. > 1 means the pull-down captured CpG-rich DNA | `relh` ≥ `enrichment_min_relh` (1.0) |
| `cpgs_1x_pct`, `cpgs_5x_pct` | the genome's CpGs covered by ≥ 1 / ≥ 5 fragments | reported |
| `saturation_r`, `saturation_est` | MEDIPS.saturation's idea: Pearson r of log tile counts between two random halves of the library, and the Spearman-Brown estimate for the whole library (`2r / (1 + r)`); near 1 means more sequencing would not change the profile | reported |
| `tiles_with_reads_pct` | tiles with ≥ 1 CpG that have ≥ 1 fragment | reported |
| `gene_body_obs_exp`, `intergenic_obs_exp` | `obs_exp` of all gene bodies and of all intergenic DNA together. In a mosaic methylome gene bodies should be well above intergenic DNA; a ratio near 1 is a failed or non-specific pull-down | reported |

The thresholds are first guesses; set them per study in its config when the
first real study has run, and record why.

## Recording a run: `whale record-enrichment`

On a study branch, a successful run calls `whale record-enrichment` (the
counterpart of `whale record-run`; `record-run` refuses an enrichment run and
vice versa). It commits, under `data/studies/<STUDY_ID>/`:

| File | What |
|---|---|
| `qc_summary.tsv`, `run_manifest.json` | as produced by the run |
| `enrichment/{gene_body,promoter}.{count,cpm,obs_exp}.tsv.gz` | gene × library matrices, every gene of the bundle, libraries that passed QC |
| `enrichment/cpgs.tsv.gz` | each gene's CpGs in gene body and promoter |
| `methylome_manifest.tsv` | URL, size and SHA-256 of every per-sample file on the data host |
| `study_summary.yaml` | `tier: enrichment`, a `comparability` statement, the method, QC medians, gene bodies vs intergenic, and how many genes are enriched |
| `status.yaml` | aligned → called → qc_passed → summarized |

The per-sample tile tables, gene tables and QC reports go to the data host
(`methylomes/<STUDY_ID>/<sample>/`), as methylomes do. Matrices go in
`enrichment/`, never `genes/`, so nothing that reads base-resolution gene
tables can pick them up by accident. No CpG methylation percentage is
written for these studies, anywhere.

## What an enrichment result can and cannot be compared with

| Comparison | Allowed? |
|---|---|
| Libraries of the same study: ranks of genes / tiles, fold differences between groups (on `cpm` or `obs_exp`, or counts in a count model) | yes |
| Two enrichment studies: rank agreement (Spearman) of gene-body `obs_exp`, and whether the same genes are top-enriched | yes, as ranks; MeDIP and MBD capture differently, so never as values |
| An enrichment study vs a base-resolution study of the same tissue: rank concordance of gene-body enrichment with gene-body methylation % | yes, as a **validation** of the assay, reported separately |
| Pooling enrichment values with % methylation in the atlas tables, variance partitioning, or the core/variable methylome | **no** (PLAN.md §4.1) |
| Reading `obs_exp` or `cpm` as "x % methylated", or averaging it with WGBS / EM-seq values | **no** |
| Expression analyses (PLAN.md §7.2) | only within the study, with methylation entered as an enrichment rank, and never pooled across tiers |

## Where enrichment studies are shown

- **Atlas**: `whale build-atlas` leaves every enrichment-tier study out of the
  pooled base-resolution tables (by the summary's `tier`, or by its chemistry's
  tier in the vocabulary), and lists it under *excluded* with the reason
  `enrichment tier (MeDIP/MBD)…`; it does not vote on the atlas assembly
  either. `tests/test_atlas.py::test_enrichment_studies_never_enter_the_pooled_tables`
  holds this.
- **Dashboard**: the study's finding is written from its summary as an
  enrichment-tier finding (QC, gene bodies vs intergenic, enriched genes)
  and says it is relative enrichment, not % methylation.
- **Next** (not built yet): an enrichment section on the study page (gene ×
  library `obs_exp` with the gene names, like the base-resolution gene
  table) and, once there are two enrichment studies or one with a
  base-resolution partner of the same tissue, an atlas panel of rank
  concordance only, labelled by tier.
