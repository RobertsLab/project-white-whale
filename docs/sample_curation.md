# Curated sample metadata

The atlas compares methylation across tissues, life stages, sexes, ploidies
and pooling designs, so every library needs those attributes in the shared
vocabularies. A study's sample sheet cannot supply them: it is the deposit as
ENA and BioSample have it, its columns differ from study to study, and it is
regenerated from ENA whenever the study is fetched again. The curated
attributes live in their own file beside it:

```
data/studies/<STUDY_ID>/sample_curation.tsv
```

One row per library of `samplesheet.tsv`, **in the same order, every library**
(including ones that fail QC), tab-separated:

| column | content |
|---|---|
| `sample_id` | as in the sample sheet |
| `animal_id` | the animal (or pool) the library came from; libraries of one animal (two tissues, two time points, technical replicates) share it. Empty when not known |
| `tissue` | a `tissue.yaml` term, or `unknown` |
| `life_stage` | a `life_stage.yaml` term, or `unknown` |
| `age` | as stated (`10 months`, `48 hpf`), else empty |
| `sex` | `female`, `male`, `hermaphrodite` or `unknown` |
| `ploidy` | `diploid`, `triploid`, `tetraploid` or `unknown` |
| `sample_unit` | a `sample_unit.yaml` term, or `unknown` |
| `pool_size` | animals per library when `sample_unit` is `pool` (at least 2), else empty |
| `population` | short label of the source population, line or family |
| `treatment` | the library's experimental group in the paper's terms |
| `confidence` | `high` (stated per sample), `medium` (stated for the study, applied to each library), `low` (inferred) |
| `sources` | where each value came from: `field: source` pairs separated by `; ` |
| `notes` | anything a reader must know about the row |

`whale validate-curation` (also run by `make validate`) checks every table:
the same libraries in the same order as the sample sheet, vocabulary terms,
pool sizes, and a source on every row.

## How the atlas uses it

`whale build-atlas` takes animal, tissue, life stage, sample unit, pool size,
ploidy, sex and population from the curated table for every library that has
a row, and records `metadata_source = curated` in `sample_metadata.tsv`
(otherwise `deposit`, with the registry's study-level terms or the deposit's
own columns, as before). `unknown` becomes an empty cell. Units (libraries
summed before comparison) still come from BioSamples: `animal_id` is
reported, not used to merge libraries, because two tissues or two time points
of one animal are different samples.

## How a table is made

Two curators read the primary sources independently: the ENA run records,
every BioSample record, and the paper with its methods and supplements. Neither
reads the study YAML or the other's work, so where they agree it is because the
sources say so, not because one copied the other. A third person compares the
two field by field and settles each disagreement from the sources, and the
settled table is committed with the curation notes in the study YAML. The
instructions each curator gets are in [sample_curation_brief.md](sample_curation_brief.md)
(its output is the table above plus a `study.json` of study-level facts:
publication, library kit, published DML/DMR tables, paired RNA-seq,
discrepancies). Values neither source states are `unknown`, never a guess;
an inference is marked `low` and explained in `notes`.
