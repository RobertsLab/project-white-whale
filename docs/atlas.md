# The atlas: pooled tables, landscape and concordance

The first cross-study views (PLAN.md §6.2, §7.1 items 1 and 4), built by
`whale build-atlas` from the studies merged to `main` and published with the
dashboard at `atlas/landscape.html` and `atlas/concordance.html`.

## What is pooled

Every study on `main` with a `study_summary.yaml` and gene tables
([gene-level methylation](gene_methylation.md)), except:

- gene tables not SNP-masked ([the C/T SNP mask](snp_mask.md)): re-run the
  study with `WHALE_SNP_MASK` and `-resume`, which redoes only the gene step;
- gene tables on another assembly than most studies.

Left-out studies are listed on both pages with the reason. Study branches are
not pooled: a study joins the atlas when its PR merges, and the Pages workflow
rebuilds the atlas on every push to `main`.

**Samples, not libraries.** Libraries of one BioSample are technical replicates
of one animal or pool; they are merged (methylated and total calls summed)
before anything is compared, so they do not count as independent samples. A
study whose libraries share one BioSample for convenience carries the
`shared_biosample` quality flag, and its libraries stay separate.

## Outputs (`site/atlas/data/`, built in CI, not committed)

| File | What |
|---|---|
| `tables/{gene_body,promoter}.{meth,calls,cpgs}.tsv.gz` | genes x libraries, every pooled study side by side, columns in `sample_metadata.tsv` order |
| `tables/{gene_body,promoter}.meth_pct.tsv.gz` | methylated / all calls, % (empty: no call) |
| `tables/sample_metadata.tsv` | one row per library: unit (sample) id, study, BioProject, chemistry, tissue, life stage, sample unit, pool size, ploidy, sex, condition, factor levels, assembly, SNP mask, run |
| `tables/sample_qc.tsv` | every library's `qc_summary.tsv` row, with its study |
| `atlas.json` | what the pages draw |

The tables are derived entirely from the committed per-study tables, so they
are rebuilt rather than committed (a study PR may only change its own files).
Build them locally with:

```bash
whale build-atlas && python3 -m http.server --directory site 8765
# then open http://localhost:8765/atlas/landscape.html
```

## What the pages compute

A gene counts in a sample when it has >= 5 CpGs with a call there (the
dashboard's rule). Thresholds: "unmethylated" < 10%, "methylated" >= 50%.

- **Landscape**: per sample, the distribution of gene-body and promoter
  methylation over genes (2% bins), the median, and the fractions under 10% and
  at or over 50%; per study, the medians of those over its samples. One panel
  per study, with the median of all pooled samples as a reference line.
- **Concordance**: Spearman's rho of gene-body methylation between every two
  samples (ties averaged, over the genes measured in both), summarised within
  each study and between each pair of studies. Per-sample rho is held down by
  sampling noise, so a shallow study agrees less with itself; read a
  between-study value against the two studies' within-study values. Study
  against study: each study's mean per gene (kept where at least half its
  samples measure it), their rho, a density grid, and how many genes both call
  unmethylated, both methylated, or call oppositely.
- **Design caveats**, generated from the pooled studies: any tissue, life
  stage, chemistry or sample unit that comes from a single study is
  inseparable from that study, and a factor every study shares says nothing.

Every statement carries the number of studies and samples it rests on.

## Not yet

Variance partitioning and the core / variable classification (§7.1 items 2-3)
need more than one study per tissue and life stage; 1 kb tiles need the
`cpg_calls` table on the data host.
