# Replication reports

For every study, White Whale asks whether the paper's findings hold when its
data are re-derived from raw reads with the shared pipeline (PLAN.md §7.1,
item 5). Each answer is a per-claim verdict, from rules simple enough to
disagree with.

```
data/studies/<STUDY_ID>/claims.yaml       the paper's claims (written by hand or by an agent)
data/studies/<STUDY_ID>/replication.md    the report, written by `whale replicate`
data/studies/<STUDY_ID>/replication.json  the same, for the dashboard and the atlas
```

## Claims

`claims.yaml` lists the paper's sources, optional `caveats` shown at the top
of the report, and the claims. Each claim paraphrases one finding and gives
it a `kind` and a `check` that a program can test:
[replication_claims_brief.md](replication_claims_brief.md) defines them. The
claims are written **from the paper alone**, without sight of our results, so
the list is not shaped by what replicated. Claims our data cannot reach stay in
the list as `not_testable` with a reason, so a report shows what replication
could not cover as well as what it could.

## Verdicts

`whale replicate [STUDY_ID ...]` tests every claim against the study's
committed results: per-library CpG methylation (`qc_summary.tsv`), the per-gene
tables (`genes/`) and, for differential claims, the shared DSS model
(`differential/`, [differential.md](differential.md)).

| kind | test | replicated when |
|---|---|---|
| `global_level` | mean CpG methylation per sample (or, with `region`, the mean per-gene level) | within 2 points (partly within 5) |
| `group_difference` | difference of the sample means, permutation test | same direction and p < 0.05, or complete separation when the group sizes cannot reach 0.05; a claim of *no* difference is inconclusive at such sizes |
| `region_contrast` | CpG calls pooled over all genes, gene body vs promoter vs whole library | same direction |
| `methylated_fraction` | share of genes at or above the paper's cut-off | within 5 points (partly within 15); inconclusive when the paper has no cut-off; not testable per CpG |
| `dm_count` | genes at FDR < 0.05 in the matching job | genes vs genes: within 3x; other units: presence only, at most partly |
| `dm_direction` | more significant genes up or down | the same majority, from at least 5 genes |
| `dm_genes` | each named gene, matched by symbol or product description | FDR < 0.05 in the claimed direction, corrected for the number of matching paralogs; the claim: most genes |
| `clustering` | adjusted R² of the factor on PC1/PC2 of the largest PCA in which it varies | >= 0.5 (partly >= 0.2); for a claim of no separation, < 0.2 |

Two conditions cap a verdict at **partly**:

- **The paper's contrast is narrower than any job.** For example, the paper
  tests time within one phenotype while the shared model pools phenotypes.
  The nearest broader job is used and the verdict is marked approximate.
- **The job's test statistics are inflated** (genomic λ > 1.3), so its
  significant-gene count is not trusted.

Other verdicts are `inconclusive` (tested, too little signal), `not_testable`
and `pending` (waits for a missing input, usually the differential model).

## Reading the numbers

Our results are per gene on GCF_963853765.1. Most papers count CpGs or regions
on older assemblies, so counts are compared by size only when both count genes.
Gene names are matched by symbol or product description, not by coordinates.
Oyster gene families are large (85 annotated genes match "heat shock 70 kDa
protein 12A"), so a named gene is only replicated when it survives correction
for the number of genes it could be. A re-run of a study changes its verdicts:
run `whale replicate` again after `whale differential`.
