# Gene-level methylation

Every library is summarized over the same genes, so studies can be compared
gene by gene (PLAN.md §5 SUMMARIES, §6.2 `gene_features`).

## How it is made

`FEATURE_SUMMARY` (workflows/methylation/bin/methylome_features.py) reads a
library's strand-merged CpG calls and the reference bundle's
`features/genes.bed`, `promoters.bed`, `exons.bed` and `introns.bed`
([reference bundle](reference_bundle.md)). Promoters are 1 kb upstream to
200 bp downstream of the TSS. Every gene in `genes.bed` (33,331 for
GCF_963853765.1, annotation RS_2024_06) gets a row in every class, covered or
not, so every sample of every study has the same rows.

## Where it is

| File | Where | What |
|---|---|---|
| `<sample>.features.tsv.gz` | gannet, `methylomes/<STUDY_ID>/<sample>/`, listed in `methylome_manifest.tsv` | one row per gene and class: `cpgs`, `cpgs_5x`, `meth`, `unmeth`, `meth_pct`, `mean_meth_5x_pct` (~1.7 MB) |
| `genes/{gene_body,promoter}.meth_pct.tsv.gz` | committed, `data/studies/<STUDY_ID>/` | genes x libraries: methylated / all calls, % (empty: no call) |
| `genes/{gene_body,promoter}.cpgs.tsv.gz` | committed | CpGs with at least one call |
| `genes/{gene_body,promoter}.calls.tsv.gz` | committed | methylated + unmethylated calls |

The matrices hold the libraries that passed the QC gate, in `qc_summary.tsv`
order. `whale record-run` builds them when the run ends and adds a `genes`
block to `study_summary.yaml` (genes with >= 5 covered CpGs in every library,
their median methylation and how many are >= 50%), which the dashboard's
finding quotes. No coverage cut-off is applied to the matrices: choose one
from `cpgs` and `calls` for the question at hand.

## Before comparing studies

- **No C/T SNP mask yet** (PLAN.md §6.1). A C→T SNP reads as an unmethylated
  C, so populations with different SNPs differ in apparent methylation.
- **Coverage.** Only CpGs with a call are counted; a gene's `cpgs` is not the
  number of CpGs it holds. ROBERTS2022 averages ~5.5 calls per CpG.
- **Design.** Tissue, chemistry, pooling and ploidy differ between studies
  and must be modelled. A triploid's value averages three alleles.
- **Mosaic methylation.** Oyster genes tend to be either methylated or not,
  so a methylated / unmethylated call may compare across studies more
  safely than small percentage differences.
