# Re-derived differential methylation

`whale differential` runs each study's registered contrasts through one shared
model ([PLAN.md §6.3](../PLAN.md#63-re-derived-differential-results)), so "does
this replicate" asks the same question of every study. It runs after a study
reaches `summarized`, either on the study branch or on main once the study has
merged, and writes only to `data/studies/<STUDY_ID>/differential/`. Its PR
therefore passes the scope check like any other study PR.

```bash
whale differential --all --plan                 # what would be fitted, and what is skipped and why
whale differential --study WANG2020_TIDAL_HEAT_WGBS            # fit it (Docker; ~1-5 min per study)
whale differential --study <ID> --engine apptainer             # on klone
```

## What is fitted

- **Units, not libraries.** Libraries of one BioSample are summed into one unit
  first, as in the atlas, unless the study is flagged `shared_biosample`.
- **Sites are genes.** In each region (gene body; promoter, 1 kb upstream to
  200 bp downstream of the TSS), a gene's summed methylated and total calls form
  one site. A gene is tested only if every unit in the job has ≥ 5 covered CpGs
  and ≥ 10 calls.
- **Test.** DSS 2.58 beta-binomial regression (`DMLfit.multiFactor`, Wald test
  per coefficient), in the pinned container
  `bioconductor-dss:2.58.0--r45h01b2380_0`. Benjamini–Hochberg is applied
  within each job, region and coefficient, with significance at FDR < 0.05.
  `diff_*` columns give the difference in raw group means in percentage points
  (level minus reference).
- **Sensitivity.** Each job is fitted again without units lying more than
  3 SD from the mean on gene-body PC1 or PC2. The counts are recorded in
  `tests.tsv` and `summary.yaml`, and the genes in `significant.tsv`.

### Default models

With no `differential` block in the study YAML:

1. One additive model over every **testable** factor. A testable factor has a
   `reference_level`, and both the reference and at least one other level have
   ≥ 2 units. Each such factor is tested.
2. If exactly two factors are testable and every combination of their levels has
   ≥ 2 units, the full model `~ a * b` is also fitted, and its interaction is
   tested.

A factor with `reference_level: null` (for example, population as a blocking
variable) is left out of the default model. Everything left out is listed under
`skipped` along with the reason.

### Named models: repeated measures and nested designs

The default model treats every unit as a separate animal. When animals were
sampled more than once, or the design is nested, name the models in the study
YAML instead. VALDIVIESO2025 is an example: each oyster was biopsied in gill
and mantle at T0 and T1, and resistant and susceptible animals are paired by
population.

```yaml
differential:
  individual: [population, phenotype_class, code]   # sample-sheet columns that together name one animal
  models:
    - id: phenotype                                 # between animals, so within each tissue x time
      formula: "~ population + phenotype_class"
      within: [tissue_type, time]
    - id: time                                      # within animals: a fixed effect per animal
      formula: "~ individual + time"
      within: [tissue_type]
      test: [time]
```

`within` fits the model separately for each combination of those factors'
levels; each combination becomes a job named `<model>.<level>…`. `test`
defaults to every factor in the formula that has a reference level, plus the
formula's interactions. `whale validate-study` checks that every term is a
factor id (or `individual`).

## Output (`data/studies/<ID>/differential/`)

| File | Contents |
|---|---|
| `design.tsv`, `jobs.tsv` | The plan: one row per unit (its libraries and levels); one row per job (formula, tests, level order, units) |
| `<job>.<region>.tsv.gz` | Every tested gene: symbol, means overall and per level, then `diff_`, `stat_`, `p_`, `fdr_` per coefficient |
| `tests.tsv` | Per job, analysis (`all_units`, `without_pca_outliers`), region and coefficient: genes tested, significant genes (up, down, and those with \|diff\| ≥ 10 pp), λ |
| `significant.tsv` | Every gene at FDR < 0.05, from both analyses |
| `pca.tsv` | Gene-body PC1 and PC2 per unit and job, with the outliers flagged |
| `summary.yaml` | All of the above in one place, with the method, the SNP mask, the run behind the gene tables and the study's quality flags |

## Reading the results

- **Check λ before the counts.** λ near 1 means the test is calibrated.
  Promoters sit around 0.8–0.9 because most of them are near 0% methylated,
  which leaves the beta-binomial little to work with. A gene-body λ well
  above 1 (for example 1.4–1.9 for the tidal-origin and within-animal contrasts)
  means either a broad true difference or unmodelled structure: pooling,
  family, private C/T SNPs, or batch. Look at `pca.tsv` and the study's
  `confounded_design` notes before reading those counts as biology.
- **Quality flags carry over.** A significant origin effect in a study flagged
  `confounded_design`, with one population per class, is a population
  difference and not an origin effect.
- **Fewer than two units per level means no test.** Neither WANG2014 (one
  animal per group) nor PESTO2022 (no factors registered yet) has a job.
