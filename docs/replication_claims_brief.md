# Claims extraction — brief

White Whale re-processes public oyster (*Magallana gigas*) DNA-methylation
data from raw reads with one pipeline and one reference genome, then asks of
each study: **do the paper's findings replicate when the data are re-derived?**
Your job is the first half: write down what the paper claims, precisely enough
that a program can test each claim against our re-derived results. You do not
see those results, and must not look for them: a claim list written with the
answers in view would be biased.

## What you may read

- The paper(s) named in your task, their methods, results, figures, tables and
  supplements (PubMed / PMC, publisher pages, HAL, preprints, theses, the
  authors' public repositories). Use local copies where your task gives them.
- The study's registry entry, for its design factors only:
  `registry/studies/<STUDY_ID>.yaml`
  (the `factors:` block: `factor_id`, `levels`, `reference_level`).
- The study's sample sheet header and a few rows:
  `data/studies/<STUDY_ID>/samplesheet.tsv`.

## What you must not read

Anything else in the repository: not `study_summary.yaml`, `qc_summary.tsv`,
`genes/`, `differential/`, `analysis/`, `docs/`, other studies, git history,
or other folders under the scratchpad. Do not modify the repository.

## What our re-derived results can test

Know this so you can classify claims, not so you can look at them.

- Per library: overall CpG methylation (%, all CpG calls pooled), non-CpG
  methylation, CpGs covered.
- Per gene (NCBI annotation of GCF_963853765.1, not the paper's genome):
  methylation of the **gene body** and the **promoter** (1 kb upstream to
  200 bp downstream of the TSS), as summed calls per gene.
- Per registered contrast: one shared beta-binomial model (DSS) **per gene**,
  gene body and promoter separately, FDR < 0.05; counts of significant genes,
  how many are higher vs lower in the non-reference level, each gene's
  difference in percentage points, and its symbol and description.
- PCA of gene-body methylation across samples.

Not available: per-CpG or per-region (DMR) calls, exons/introns/intergenic/
repeats as separate matrices, the paper's genome coordinates, RNA-seq.

## Output

Write `claims.yaml` to your output folder:

```yaml
study_id: WANG2020_TIDAL_HEAT_WGBS
sources:
  - id: P1
    citation: "Authors (year). Title. Journal vol:pages."
    doi: "..."
    read: "full text and supplements"     # what you actually read
claims:
  - id: C1
    source: P1
    where: "Results 3.1; Fig. 2a"          # where in the paper
    statement: "Your paraphrase of the claim, one sentence, with its numbers."
    kind: global_level
    check: {level_pct: 12.1, groups: all, measure: "mean CpG methylation level across samples"}
```

**Paraphrase; never paste the paper's sentences.** Numbers, gene names and
table references are fine.

### Kinds and their `check` fields

Use the registry's `factor_id` and level strings exactly (for example
`factor: ph`, `level: "7.7"`), and say which level is higher with `level` =
the non-reference level, as the registry defines the reference.

| kind | the claim is about | `check` fields |
|---|---|---|
| `global_level` | the overall CpG methylation level | `level_pct` (number), `groups`: `all` or `{factor: level}`, `measure` (the paper's definition) |
| `group_difference` | one group's overall methylation is higher than another's | `factor`, `higher` (level), `lower` (level), optional `difference_pct` |
| `region_contrast` | gene bodies vs promoters (or vs genome) | `higher`: `gene_body`/`promoter`, `lower`: `gene_body`/`promoter`/`genome`, optional values |
| `methylated_fraction` | the share of genes (or CpGs) that are methylated | `unit`: `genes`/`cpgs`, `region`, `threshold_pct` (the paper's cut-off, or null), `fraction_pct`, `definition` |
| `dm_count` | how many differentially methylated CpGs/regions/genes a contrast gave | `factor`, `level`, `stratum` (`{factor: level}` or null), `n`, `unit`: `DMCs`/`DMRs`/`DMGs`, `method` (tool, thresholds) |
| `dm_direction` | most differences go one way | `factor`, `level`, `stratum`, `direction`: `hyper`/`hypo` (in `level` vs the reference), optional `share_pct` |
| `dm_genes` | named genes are differentially methylated | `factor`, `level`, `stratum`, `genes`: list of `{name: symbol or product description, published_id, direction: hyper/hypo/null, region: gene_body/promoter/null}` (at most 25, the ones the paper highlights) |
| `clustering` | samples group by a factor in a PCA or clustering | `factor`, `component` (e.g. `PC1`), `statement` |
| `not_testable` | anything else | `reason`: `needs_cpg_level` / `needs_other_regions` / `needs_rnaseq` / `needs_coordinates` / `other_data` / `not_quantitative`, `note` |

Rules:

1. **Cover the paper's main findings**: everything in the abstract, plus the
   central quantitative results. 8-25 claims per paper is typical. Include
   claims you expect cannot be tested, as `not_testable`, so the report shows
   what replication could not reach.
2. **A claim must be falsifiable by its `check`.** If the numbers you would
   need are not in the paper, say so in `statement` and keep the fields you
   can fill.
3. **One claim, one kind.** Split compound sentences.
4. **The paper's contrast may not match a registered factor.** Then use
   `not_testable` with `reason: other_data` and explain, rather than forcing it.
5. When you can read only the abstract, extract from the abstract and say so
   in `sources[].read`.
6. Finish with a short reply: the output file, the papers read and how, the
   number of claims by kind, and anything that made extraction hard.
