# White Whale — plan

**A harmonized public DNA-methylation resource for the Pacific oyster
(*Crassostrea / Magallana gigas*).**

Every public bisulfite / EM-seq dataset for the species, re-processed from raw
reads with one pipeline against one reference, into one comparable set of
per-sample methylomes — and, where paired RNA-seq exists, into a direct test of
what methylation does to gene expression.

This plan takes its structure from [AREE](https://github.com/sr320/AREE)
(registry → workflow → harmonization → synthesis → site, with the branch-per-
study claim model for running on several machines) and narrows it to one
species, one assay family, and two questions. It replaces the previous White
Whale repository, which was a dataset *catalogue* with a downloader; that
catalogue is the seed of the new registry (see §4).

---

## Status (2026-11-25)

Sections 1, 2 and 4–8 are the design. This section, §9 and §10 say how far it
has been built; the live numbers are on the
[dashboard](https://robertslab.github.io/project-white-whale/).

- **Built and in use:** registry, claims, compute profiles (desktop, raven,
  klone), the methylation workflow, reference bundle, C/T SNP mask v1,
  gene-level tables, the shared DSS model, replication reports, sample
  curation, the dashboard and three atlas pages, `whale mcp`, and the
  autopilot's deterministic tick.
- **Backlog:** 8 of 21 studies merged as `summarized`, 3 in flight, 10
  unclaimed (§4.2).
- **Goal 1:** first cross-study views are live (landscape, concordance, gene
  methylation classes), from per-gene tables only. Per-CpG and 1 kb tile
  tables, variance partitioning and genome tracks are not started.
- **Goal 2:** the RNA-seq workflow is merged and tested on fixtures. No
  study's RNA-seq has been run and no methylation × expression analysis
  exists yet.

**Needs attention, in order:**

1. **Re-run two merged studies whose committed results come from flawed
   runs.** ROBERTS2022_PH_PLOIDY_WGBS was aligned in directional mode (run
   `curious_mestorf`, 30% mapping, about half depth); its config now has
   `--non_directional`. WANG2023_TIDAL_TRANSGEN_WGBS was not clipped on read 2
   (run `elated_pesquet`, CpG methylation ~1.5 points high, conversion ~97%);
   its config now clips. Both feed the atlas, the DSS results and the
   replication verdicts, which must be re-derived afterwards.
2. **Decide on SNP mask v2.** Mask v1 was called from ROBERTS2022 alone, at
   half depth. §6.1 asks for the deepest studies.
3. **Resolve the low QC pass rates** (§10, open decision 1):
   GAWRA2023 132 of 246, PESTO2022 13 of 48, ROBERTSLAB2013 1 of 14.
4. **Finish the in-flight studies:** ROBERTSLAB2013 (at `qc_passed` on Hete,
   which has no gannet key), FALLET2022 (staged on klone), OUC2024_SEX_REVERSAL
   (registered on raven).
5. **Start Goal 2** with WANG2020_TIDAL_HEAT_WGBS, the only merged study with
   same-animal RNA-seq, and claim the unclaimed studies that have RNA-seq
   partners next.

---

## 1. Goals

### Goal 1 — Comparative characterization of the methylome across studies

Describe DNA methylation across the *C. gigas* genome using every public
dataset at once, on a common footing, so that statements are about the species
and not about one study:

- the global landscape: CpG methylation level distribution, bimodality,
  fraction of methylated CpGs, per-chromosome and per-feature profiles
  (promoter / 5'UTR / exon / intron / 3'UTR / intergenic / repeat / TE);
- what is **invariant** across tissues, life stages, populations, ploidy and
  conditions (the "core methylome"), and what is **variable**, and at which
  resolution (CpG, 1 kb tile, gene);
- partitioning of methylation variance into study, tissue, life stage,
  genotype/population, treatment and residual;
- how the published studies' own claims look when re-derived with a shared
  method (do the DML/DMR sets replicate within the pooled data?).

### Goal 2 — The role of DNA methylation in regulating gene expression

For every dataset with matched (same animals, or same design) RNA-seq:

- gene-body methylation vs expression level (the classic monotonic / bell-
  shaped relationship) and vs expression *variability* across samples
  (transcriptional noise);
- promoter methylation vs expression, by promoter CpG density class;
- exon-level methylation vs exon inclusion / alternative splicing;
- within-study changes: do DMLs / DMRs under a treatment co-locate with DEGs,
  in which direction, and is that consistent across studies;
- plasticity: is the magnitude of expression change under stress predicted by
  baseline gene-body methylation (the stable-genes / plastic-genes hypothesis);
- the same relationships across the meta-set, with study as a random effect,
  so the answer is "across N studies" and not "in one".

### What this project is not

- Not a resilience-biomarker engine; AREE already does that for phenotype
  contrasts, and the two methylation studies it has registered will be shared,
  not duplicated.
- Not a browser of other people's tables. Nothing enters the comparative
  analysis unless it was re-processed from raw reads here. Published
  processed results are kept only as a replication target.

---

## 2. Design principles (carried over from AREE)

1. **One reference.** Every sample maps to `GCF_963853765.1` (xbMagGiga1.1,
   annotation `RS_2024_06`). Coordinates and gene IDs are comparable by
   construction. Lift-over is never used.
2. **Raw reads only, one pipeline, pinned versions.** Container digests and
   parameters are part of every output manifest.
3. **Registry first, data second.** A study is registered (YAML, validated
   against a schema and controlled vocabularies) before a byte is downloaded.
   The sample sheet is read from ENA/SRA metadata, not transcribed from the
   paper.
4. **Explicit missingness and confounding.** Pooled libraries, shared
   BioSamples, one-population-per-class designs, EM-seq vs WGBS vs RRBS
   chemistry, missing spike-ins — all recorded as flags, never silently
   dropped or silently merged.
5. **Provenance everywhere.** Source accession, checksums, versions, machine,
   date, who.
6. **Generated outputs are not committed.** Small summaries (gene × sample
   matrices, QC tables, dashboard JSON) are committed; large per-sample files
   live on a data host and are referenced by manifest + checksum.

---

## 3. Repository layout

As built (2026-11-25):

```
project-white-whale/
├── README.md, PLAN.md, LICENSE, VENDORED_FROM_AREE.md
├── pyproject.toml, Makefile
├── schemas/                      study.schema.json, status.schema.json
├── registry/
│   ├── backlog.tsv               every public study, claimed or not
│   ├── studies/<STUDY_ID>.yaml   one per BioProject (or per coherent sub-project)
│   └── controlled_vocabularies/  assay_chemistry, tissue, life_stage, stressor, sample_unit,
│                                 quality_flags, pipeline_stage, compute_profile, rnaseq_relationship
├── workflows/
│   ├── methylation/              Nextflow DSL2: FASTQ → per-sample CpG calls → QC → gene tables
│   ├── snp_mask/                 C/T SNP mask from pooled WGBS libraries
│   ├── rnaseq/                   paired-expression studies (fastp → salmon → gene matrices)
│   ├── enrichment/               MeDIP / MBD tier (bowtie2 → counts)
│   └── differential/             dss.R, the shared model run by `whale differential`
├── modules/                      Nextflow processes (vendored from AREE where identical)
├── config/
│   ├── base.config               profiles: desktop | raven | klone
│   ├── autopilot.yaml            per-machine limits and environments for `whale autopilot`
│   └── <STUDY_ID>.config
├── data/
│   ├── data_host.yaml            gannet
│   ├── reference/                genome_assemblies.yaml, snp_masks.yaml; the bundle itself is not in git
│   ├── studies/<STUDY_ID>/       status.yaml, samplesheet.tsv, fastq_manifest.tsv, fastq_staged.tsv,
│   │                             ena_provenance.json, run_manifest.json, qc_summary.tsv,
│   │                             study_summary.yaml, methylome_manifest.tsv, genes/,
│   │                             sample_curation.tsv, claims.yaml, replication.{md,json},
│   │                             differential/ (committed, small)
│   └── raw/, methylomes/         → this machine's drives (not in git)
├── src/whale/                    the `whale` CLI and library, one module per concern:
│   │                             cli, registry, validate, vocab, claims, scope, ena, fastq, mirror,
│   │                             reference, snp_mask, results, curation, differential, replication,
│   │                             rnaseq, enrichment, atlas, dashboard
│   ├── autopilot/                tick, decide, playbook, slurm, machines, digest
│   └── query/                    the MCP server (`whale mcp`)
├── analysis/<STUDY_ID>/          per-study analyses (so far ROBERTS2022_PH_PLOIDY_WGBS)
├── docs/                         how-tos and design notes, Markdown
├── site/                         the published site: index.html (dashboard), atlas/, genes/;
│                                 data files are generated by `whale build-dashboard` / `build-atlas`
├── tests/
└── .github/workflows/            ci.yml (tests, schema validation, scope check, workflow runs on
                                  fixtures), pages.yml (site, on push and every three hours)
```

Not built, and still wanted: `src/whale` code for Goal 2 (§7.2), analysis
notebooks for the pooled results, `CITATION.cff`, `CONTRIBUTING.md`.

---

## 4. Layer 1 — Study registry

### 4.1 Study YAML

Adopt AREE's study schema and trim it to this project. Every study records:
`study_id`, accessions (BioProject / SRA / ENA / GEO / DOI / PMID), citation,
species name as deposited, population / line / ploidy, assay chemistry
(`wgbs`, `emseq`, `rrbs`, `bs_capture`, `medip`, `mbd`), library layout, platform,
`paired_rnaseq` (BioProject + relationship: `same_animals` /
`same_design_different_animals` / `none`), per-sample attributes
(tissue, life stage, sex, ploidy, treatment, timepoint, pooling), data
availability, QC and processing status, quality flags, limitations, provenance.

Two vocabularies matter more here than in AREE and are added:

- `assay_chemistry` with a **quantitative-comparability tier**: base-resolution
  (WGBS, EM-seq, RRBS and bisulfite capture, each with a coverage caveat) vs enrichment (MeDIP, MBD).
  Enrichment data is processed and shown but never enters a base-resolution
  comparison.
- `sample_unit`: individual / pool-of-N / tissue-pooled, because pooling
  changes what variance means.

### 4.2 The backlog

The seed catalogue from the previous repository, the studies AREE had
registered, and two SRA/ENA sweeps (2026-08, 2026-10) are now
[`registry/backlog.tsv`](registry/backlog.tsv): **21 studies, about 4.5 TB of
FASTQ**. `whale backlog` prints the current state of each; the table below is
a snapshot from 2026-11-25. Where AREE has a YAML the `study_id` is kept so
the two projects can cross-reference.

**Merged (`summarized`), 8 studies**

| Study | BioProject | Chemistry | QC pass | Machine | Paired RNA-seq | Note |
|---|---|---|---|---|---|---|
| `WANG2014_FIRST_METHYLOMES_WGBS` | PRJNA173440 | WGBS | 4 / 4 | raven | related (PRJNA146329), not one-to-one | |
| `ROBERTS2022_PH_PLOIDY_WGBS` | PRJNA682817 | WGBS | 24 / 24 | klone | none | **re-run needed** (non-directional libraries) |
| `WANG2020_TIDAL_HEAT_WGBS` | PRJNA562805 | WGBS | 12 / 12 | klone | same animals (PRJNA407831) | first Goal 2 study |
| `WANG2023_TIDAL_TRANSGEN_WGBS` | PRJNA807732 | WGBS | 24 / 24 | raven | unknown | **re-run needed** (read 2 clipping) |
| `PESTO2022_METHYLSEQ_WGBS` | PRJEB58545 | WGBS | 13 / 48 | klone | PRJEB58194, pairing unknown | conversion 90–97% in most libraries; factors not registered |
| `VALDIVIESO2025_POMS_RS_EMSEQ` | PRJEB81880 | EM-seq | 40 / 40 | klone | none | |
| `DECICOMP2025_POMS_AGE_EMSEQ` | PRJEB105019 | EM-seq | 57 / 60 | klone | none (metabolomics only) | |
| `GAWRA2023_POMS_ADAPTATION_BSCAPTURE` | PRJEB60400 | bisulfite capture | 132 / 246 | klone | none | 114 libraries fail on mapping, most near 0.1% |

**In flight, 3 studies**

| Study | BioProject | Stage | Machine | Note |
|---|---|---|---|---|
| `ROBERTSLAB2013_EPIGENOMICS_WGBS` | PRJNA213124 | `qc_passed` | Hete | 1 of 14 libraries passes (the deep male gonad library of Olson & Roberts 2014); kept as a single-sample study. The gill library is MBD-enriched. |
| `FALLET2022_MICROBIAL_TRANSGEN_WGBS` | PRJNA609264 | `fastq_staged` | klone | 47 runs, 916 GB |
| `OUC2024_SEX_REVERSAL_WGBS` | PRJNA917710 | `registered` | raven | RNA-seq partner PRJNA917707 |

**Unclaimed, 10 studies**

| Suggested ID | BioProject | Chemistry | Runs | GB | Paired RNA-seq |
|---|---|---|---|---|---|
| `OUC2023_TRIPLOID_GONAD_WGBS` | PRJNA684407 | WGBS | 18 | 94 | PRJNA690125 |
| `OUC2025_VALGINOLYTICUS_GILL_WGBS` | PRJNA1113357 | WGBS | 9 | 110 | PRJNA830931 |
| `OUC2024_GROWTH_STRAINS_WGBS` | PRJNA833956 | WGBS | 6 | 67 | PRJNA1113800 |
| `PERPIGNAN2016_DIURON_PARENTAL_WGBS` | PRJNA291364 | WGBS | 6 | 89 | same project |
| `WHOI2022_GONAD_LOWPH_WGBS` | PRJNA806944 | WGBS | 8 | 139 | check |
| `UW2020_PLOIDY_DESICCATION_HEAT_WGBS` | PRJNA678408 | WGBS | 10 | 91 | check |
| `IOCAS2021_WILD_TIDAL_WGBS` | PRJNA689936 | WGBS | 6 | 52 | none |
| `CIBNOR2021_HEAT_RRSS` | PRJNA690951 | RRBS | — | — | none; import AREE's registration |
| `CEMM2023_580SPECIES_RRBS` | PRJNA802599 | RRBS | 13 | 6 | none |
| `RIVIERE2017_DEVELOPMENT_MEDIP` | PRJNA324546 | MeDIP | 21 | 8 | published; enrichment tier |

The first four unclaimed studies have RNA-seq partners and are the priority
for Goal 2. `backlog.tsv` itself has rows that predate what the registrations
found (PRJEB60400 is still listed as `POMS2023_ADAPTATION_WGBS` / wgbs,
PRJEB105019 as wgbs, several RNA-seq cells say `check`); bring it in line
with the study YAMLs.

### 4.3 Claims — running studies on several machines

Adopt AREE's model unchanged in spirit, because it needs no server and git
push is atomic:

```
whale start-study --accession PRJNA682817 --study-id ROBERTS2022_PH_PLOIDY_WGBS --machine raven
```

1. fetches `origin/main`, refuses if the study or accession is already
   registered or claimed by a pushed `study/` branch;
2. creates `study/<ACCESSION>/<STUDY_ID>` in its own worktree;
3. writes the stub YAML and `config/<STUDY_ID>.config` with the machine's
   profile pre-selected;
4. pushes the branch (**the push is the claim**) and opens a draft PR whose
   body names the machine and carries the study checklist.

Additions for this project:

- `data/studies/<STUDY_ID>/status.yaml` on the branch records stage
  (`registered → fastq_staged → aligned → called → qc_passed → summarized`),
  machine, dates, and the Nextflow run name, so the dashboard shows where
  every study is and on which computer without reading PR bodies.
- `whale claims` prints the ledger from `git ls-remote` + the status files.
- `whale release-study` drops a claim (deletes the remote branch) when a
  machine gives up, so another can take it.
- A study PR may touch only its own paths (`registry/studies/<ID>.yaml`,
  `config/<ID>.config`, `data/studies/<ID>/`, `tests/test_<id>.py`); CI
  enforces it, exactly as AREE's `check-study-scope`. Studies therefore merge
  in any order.
- Big studies can be claimed as **parts** (`--part gill`, `--part 1of4`) that
  share a study YAML; each part is its own branch and status file. Not needed
  so far: PRJEB60400 (246 runs) ran as one study.
- A merged study's remote branch should be deleted, or `whale claims` and
  `whale backlog` keep reporting the branch's last stage instead of
  `summarized`.

### 4.4 Compute profiles

| Profile | Machine | Executor | Scratch | Typical job |
|---|---|---|---|---|
| `desktop` | lab Macs / workstations | local | `$WHALE_WORKDIR` on a fast internal/SSD drive | small studies (≤ 150 GB), QC, summaries, Quarto |
| `raven` | Roberts Lab server | local, high `maxForks` | raven's local scratch | medium studies (150–500 GB), whole-study Bismark without chunking |
| `klone` | UW Hyak | Slurm `ckpt-all`, 4.5 h tasks, chunked alignment, Apptainer | `/scr` per node; work on `/gscratch/srlab/$USER/whale`, staged FASTQ on `/gscratch/scrubbed/$USER/whale` | the deep studies (PRJEB60400, PRJEB81880, PRJEB105019, PRJNA609264) |

Every profile reads the same four environment variables
(`WHALE_WORKDIR`, `WHALE_RAW`, `WHALE_REFERENCE`, `WHALE_METHYLOMES`) so a
config file never contains a machine-specific path. AREE's Klone lessons
(`NXF_TEMP` not `TMPDIR`, link directories not brace patterns, chunked
alignment at ~20 M pairs per task, ENA at ~12 MB/s through the proxy) are
carried over verbatim into `docs/running_on_klone.md`.

---

## 5. Layer 2 — The methylation workflow

One Nextflow pipeline, `workflows/methylation`, FASTQ in, **per-sample
methylome out**. The deliverable unit is the sample, not the contrast — that
is the main departure from AREE, whose methylation workflow ends in a DMR
table for one comparison.

```
FASTQC → TRIM_GALORE (chemistry-aware: EM-seq vs WGBS vs RRBS presets)
  → BISMARK_ALIGN (--score_min L,0,-0.6; chunked on klone)
  → DEDUPLICATE (skipped for RRBS)
  → METHYLATION_EXTRACTOR (M-bias → per-study trim decision recorded)
  → CONVERSION_EFFICIENCY (lambda / pUC19 spike-in if present; else CHH-based estimate, flagged)
  → SNP_MASK (remove CpGs overlapping C/T SNPs from the project-wide mask, §6.1)
  → PER_SAMPLE_CALLS: CpG_report.txt.gz (strand-merged), bedGraph, bigWig
  → SUMMARIES per sample: 1 kb tiles, gene features (gene body / promoter / exon / intron / UTRs), repeats
  → QC_GATE (coverage, mapping rate, conversion, duplication, M-bias, sex check) → qc_summary.tsv
  → MANIFEST (versions, params, input MD5s, machine) + study HTML report
```

Vendor `bismark_align.nf`, `bismark_chunked.nf`,
`bismark_methylation_extractor.nf`, `trim_galore.nf`, `fastqc.nf` and
`collect_versions.py` from AREE with attribution; write the rest here.
Differential analysis (methylKit / DSS) is **not** in the pipeline; it is
downstream, run from the harmonized calls, so it can be rerun with a shared
model for every study (§6.3).

Per-sample outputs (~0.5–2 GB each) go to `WHALE_METHYLOMES` on the machine
that ran them, then are synced to the lab data host (gannet, `data/data_host.yaml`) under
`methylomes/<STUDY_ID>/<sample>/`; `methylome_manifest.tsv` (committed)
holds the URL + SHA-256 + size of every file. `whale fetch-methylomes
--study X` restores them on any machine.

---

## 6. Layer 3 — Harmonization and the shared matrices

### 6.1 Reference bundle and derived tracks (built once, checksummed)

- genome FASTA + GTF from NCBI, verified against NCBI MD5s (`whale
  fetch-reference`, from AREE);
- CpG index (every CpG, strand-merged, with local CpG o/e and GC);
- feature BEDs: gene body, promoter (−1 kb/+200 bp of TSS, plus a CpG-density
  class), exons, introns, UTRs, intergenic, repeats / TEs (RepeatMasker of the
  assembly, or NCBI's masked track);
- **C/T SNP mask**: the genome is highly heterozygous and a C→T SNP is
  indistinguishable from unmethylated C. Build the mask from the WGBS data
  itself (Bis-SNP / BS-SNPer on the deepest studies) and union it with any
  published SNP panel on this assembly. Record its version in every manifest.

### 6.2 Matrices (the atlas tables)

For every base-resolution sample that passed QC:

| Table | Rows × cols | Size | Where |
|---|---|---|---|
| `cpg_calls` | ~10–15 M CpGs × samples (count M, count U) | tens of GB | data host, Parquet partitioned by study |
| `tiles_1kb` | ~550 k tiles × samples (% + coverage) | ~1–2 GB | data host, Parquet |
| `gene_features` | ~30 k genes × feature class × samples | tens of MB | built in CI from the committed per-study tables, published with the site |
| `sample_qc` | samples × QC metrics + flags | small | built in CI, published |
| `sample_metadata` | samples × harmonized attributes (tissue, stage, sex, ploidy, population, treatment, chemistry, pool size, study) | small | built in CI, published |

Every row carries `study_id`, `chemistry`, and the manifest hash, so a table
can always be filtered to "WGBS only" or "individual animals only".

**Status:** `gene_features`, `sample_qc` and `sample_metadata` are built by
`whale build-atlas` and published. `cpg_calls` and `tiles_1kb` are not built;
per-sample calls are on gannet but not pooled. Their absence is why many
replication claims (per-CpG counts, windowed density, DMR overlap) are
`not_testable`.

### 6.3 Re-derived differential results

For every study with a registered within-study contrast, run one shared
model (DSS or methylKit with the same cutoffs, same tile definition, same
minimum coverage) from the harmonized calls, and store the result next to the
authors' published DML/DMR list when one exists. This is what makes "did it
replicate" answerable, and it feeds Goal 2's DML-vs-DEG tests.

**Status:** `whale differential` runs one DSS model per gene (gene body and
promoter). Run for ROBERTS2022, WANG2020, WANG2023 and VALDIVIESO2025; not
yet for WANG2014 (1 vs 1), PESTO2022 (no factors registered), DECICOMP2025 or
GAWRA2023. `whale replicate` compares each paper's claims with the result
(docs/replication.md); six studies have a report. The model works on genes,
not CpGs or tiles, and WANG2023's test statistics are inflated (lambda 1.9).

---

## 7. Layer 4 — Analyses

### 7.1 Goal 1 — Comparative atlas (`src/whale/atlas`, `analysis/atlas_*.qmd`)

1. **Landscape per sample, then per group**: methylation-level histograms,
   bimodality index, % methylated CpGs at several thresholds, feature-class
   profiles, metagene plots (TSS → TES), chromosome-scale profiles.
2. **Core vs variable methylome**: for each gene and 1 kb tile, mean and
   variance across all samples, and across samples within tissue / stage /
   study; classify tiles as constitutively methylated, constitutively
   unmethylated, variable. Report how much of the variable fraction is
   explained by study vs biology.
3. **Variance partitioning**: linear mixed model per gene-body (or on PCs of
   the tile matrix) with study, chemistry, tissue, life stage, sex, ploidy,
   population, treatment as terms; report variance components and the PCs'
   correlation with each covariate. This is also the batch-effect diagnostic:
   if study explains most variance after accounting for tissue, say so on the
   dashboard.
4. **Cross-study concordance**: rank correlation of gene-body methylation
   between every pair of samples / studies; compare with the known result
   that gene-body methylation rank is stable across tissues in this species.
5. **Replication of published findings**: for each study with a published
   DML/DMR set, overlap with the re-derived set (§6.3) and with every other
   study's variable tiles.
6. **Genome tracks**: per-group mean bigWigs and the core/variable
   classification as a track hub for IGV / JBrowse.

**Status:** live at gene resolution: 1 (landscape page), 2 (gene methylation
classes page, genes only, no tiles) and 4 (concordance page). 5 exists per
study as replication reports, without the cross-study overlap. Not started:
3 (variance partitioning) and 6 (tracks).

### 7.2 Goal 2 — Methylation and expression (`src/whale/regulation`, `analysis/regulation_*.qmd`)

Inputs: the `gene_features` matrix, the RNA-seq quantifications from
`workflows/rnaseq` (salmon → gene and transcript TPM, same GTF), the
re-derived DML/DMR sets, and DESeq2 results for the same contrasts.

1. **Level**: gene-body methylation vs log TPM, per sample and pooled, with
   the bell-shaped fit tested against monotonic; promoter methylation vs TPM
   stratified by promoter CpG o/e.
2. **Variability**: gene-body methylation vs coefficient of variation of
   expression across replicates (transcriptional noise), within each study,
   then meta-analysed across studies.
3. **Splicing**: exon methylation vs exon inclusion (PSI from the same
   RNA-seq), per study.
4. **Co-change**: within each contrast, enrichment of DEGs among genes with
   DMLs/DMRs in body / promoter, with direction tables (hyper + up, hyper +
   down, …), and a random-effects pooled odds ratio across studies with I².
5. **Plasticity**: baseline gene-body methylation (control samples) vs |log2FC|
   under treatment in the same study; across studies with study as random
   effect.
6. Every result is reported with the number of studies it rests on and is
   never collapsed across chemistry tiers.

Analyses read only the committed small tables and the data-host files, so
anyone can re-run them without the FASTQ.

**Status:** none of 1–5 is started. `workflows/rnaseq` and `whale
fetch-rnaseq-samplesheet` / `record-rnaseq` exist (docs/rnaseq.md) and pass on
fixtures. Of the eight merged studies only WANG2020 has same-animal RNA-seq;
PESTO2022's pairing is unknown and the rest have none, so Goal 2 depends on
the unclaimed studies with RNA-seq partners (§4.2). A pooled estimate "across
N studies" needs N of at least three or four.

---

## 8. Layer 5 — Outputs

- **Site** (GitHub Pages, static pages in `site/` over generated JSON),
  regenerated on every push to `main` and every three hours. The dashboard,
  per-study findings and gene sections, and three atlas pages are live;
  regulation pages and the how-tos (which are Markdown in `docs/`) are not on
  the site. As designed:
  - dashboard: studies registered / claimed / processed, by machine; samples
    passing QC by tissue × stage × chemistry; TB processed; what the pooled
    data currently supports and does not;
  - one page per study (metadata, QC, replication of its published findings,
    link to its PR);
  - the atlas pages (Goal 1) and the regulation pages (Goal 2);
  - how-tos: adding a study, claiming on desktop / raven / klone, fetching
    methylomes, interpreting the QC gate, interpreting variance components.
- **Downloads**: the committed matrices, the track hub, the per-sample
  manifests.
- **MCP server** (`whale mcp`, built; docs/mcp.md): read-only tools over the
  summarized tables, as AREE's, so Claude can answer "which genes are
  hypermethylated in gill across every study" with study IDs cited.
- **Autopilot** (`whale autopilot`, partly built; docs/autopilot.md,
  docs/autopilot_plan.md): a scheduled tick that moves each study one stage,
  with a playbook for known failures. The deterministic tick, klone queue
  visibility and launch paths are merged; the model-driven register, diagnose
  and review steps are not.

---

## 9. Phases

The project started 2026-10-23. State as of 2026-11-25:

| Phase | Deliverable | Done when | State |
|---|---|---|---|
| **0. Scaffold** | Repo layout, `whale` CLI, schemas, vocabularies, `start-study` / claims, CI with scope check, profiles | `whale validate-study` and `whale start-study --no-push` work on a stub; CI green | **Done** |
| **1. Registry** | SRA/ENA sweep; backlog with sizes and RNA-seq partners; AREE YAMLs imported | Dashboard shows the full backlog | **Done** (21 studies; sweeps 2026-08 and 2026-10). Repeat the sweep before write-up. |
| **2. Reference + pipeline** | Reference bundle, feature BEDs; pipeline end to end on PRJNA173440 and PRJNA682817; first SNP mask | Two studies `qc_passed`, methylomes on the data host | **Done**, with two debts: mask v1 comes from one half-depth study, and no CpG index is published |
| **3. Klone** | Klone profile proven on a deep study (PRJEB81880) | Study merged as `summarized` | **Done** (six studies run on klone) |
| **4. Backlog** | Every base-resolution study processed; enrichment tier through the MeDIP/MBD workflow | All studies `summarized` or `excluded` with a reason | **In progress**: 8 merged, 3 in flight, 10 unclaimed; 2 merged studies need re-runs; no enrichment study run yet |
| **5. Goal 1 analyses** | Atlas matrices, variance partitioning, core/variable classification, track hub, atlas pages | Pages render from committed tables; findings carry study counts | **Started**: gene-level landscape, concordance and classes are live. Tiles, per-CpG tables, variance partitioning and tracks not started |
| **6. Goal 2 analyses** | RNA-seq on paired studies; the five regulation analyses; pooled estimates | Regulation pages render; each claim lists its supporting studies | **Workflow only**: no RNA-seq run, no analysis code |
| **7. Write-up** | Methods page, implementation-status table, manuscript draft from the site | — | Not started |

Phases 4–6 overlap by design: analyses are rebuilt as each study merges.

### Next steps

In this order, since each later step reuses the earlier one's output:

1. Re-run ROBERTS2022 (`--non_directional`; needs a non-preemptible klone
   partition or raven) and WANG2023 (clipped); `record-run`; re-run
   `whale differential` and `whale replicate` for both.
2. Build SNP mask v2 from several deep, well-mapped studies if v1 proves
   thin, then re-mask the gene tables of every study (§10, open decision 2).
3. Check the GAWRA2023 libraries at ~0.1% mapping against the deposit, and
   settle the gate policy for PESTO2022 (§10, open decision 1).
4. Finish ROBERTSLAB2013, FALLET2022 and OUC2024_SEX_REVERSAL.
5. Run RNA-seq for WANG2020 (PRJNA407831); then claim OUC2023_TRIPLOID_GONAD,
   OUC2025_VALGINOLYTICUS_GILL, OUC2024_GROWTH_STRAINS and
   PERPIGNAN2016_DIURON_PARENTAL together with their RNA-seq.
6. Register PESTO2022's stage, generation and exposure factors so its 37
   claims become testable; run `whale differential` for DECICOMP2025 and
   GAWRA2023.
7. Pool 1 kb tiles and per-CpG calls on gannet (§6.2), then variance
   partitioning (§7.1 item 3).
8. First Goal 2 analysis (§7.2 item 1) on WANG2020, written so further
   studies drop in.
9. Housekeeping: delete merged study and feature branches, prune worktrees,
   bring `backlog.tsv` in line with the registrations, remove the old
   downloader's `code/` leftovers.

---

## 10. Decisions

Decided 2026-10-23:

1. **Reuse of AREE code: vendor.** AREE code is copied in and adapted, with
   the source file and AREE commit recorded in
   [VENDORED_FROM_AREE.md](VENDORED_FROM_AREE.md). Revisit if both projects
   keep changing the same modules.
2. **Data host for per-sample methylomes: gannet**, at
   `/volume1/v1_web/project-white-whale`, served as
   `https://gannet.fish.washington.edu/v1_web/project-white-whale/`
   (`data/data_host.yaml`). The volume has 171 TB free.
3. **Raven** was inspected over SSH and the `raven` profile sized from it:
   24 cores / 48 threads, 247 GB RAM, Singularity CE 3.10 as the container
   engine, NVMe RAID0 work directory. See docs/compute_profiles.md.

Decided 2026-10-25:

4. **Pooled atlas tables are derived, not committed.** `gene_features`,
   `sample_qc` and `sample_metadata` are rebuilt by `whale build-atlas` in the
   Pages workflow from the per-study tables on `main`, since a study PR may
   only change its own files. See docs/atlas.md.

Decided since, in the course of the work:

5. **PRJEB60400 is one study**, `GAWRA2023_POMS_ADAPTATION_BSCAPTURE`, not
   parts. It is exome-capture bisulfite, so `bs_capture` was added to the
   base-resolution tier with a coverage caveat (§4.1).
6. **The site is static pages over generated JSON**, not Quarto; docs stay
   Markdown in `docs/`.
7. **Differential methylation is per gene** (gene body and promoter) with one
   DSS model, `whale differential`; each paper's claims are tested against it
   by `whale replicate`.
8. **Enrichment-tier data has its own workflow** (`workflows/enrichment`) and
   is excluded from the atlas.
9. **The backlog is driven by an autopilot** where rules suffice
   (docs/autopilot_plan.md); judgment steps stay with a person for now.

Open:

1. **QC gate policy for whole-study failures.** PESTO2022 fails the QC
   gate in 35 of 48 libraries (34 on the 98% conversion gate, one on mapping),
   with unstated chemistry; GAWRA2023
   passes 132 of 246. Options: keep the single gate and pool only what
   passes (current behaviour), or allow a recorded per-study gate with a
   quality flag. Affects how many samples the atlas rests on.
2. **SNP mask v2.** Which studies to call it from, and whether every merged
   study's gene tables are re-masked when the mask changes.
3. **Resolution of the pooled tables.** Whether to build `tiles_1kb` and
   `cpg_calls` (§6.2) now, or stay at gene level until the backlog is done.
   Variance partitioning, DMR replication and tracks all wait on this.
4. **Goal 2 scope.** One same-animal study is merged. Decide whether
   same-design, different-animal pairs (and PESTO2022's unknown pairing)
   count for the pooled estimates, or only for per-study results.
5. **Enrichment tier.** Whether RIVIERE2017 (MeDIP) and the MBD library in
   PRJNA213124 are worth running now that the workflow exists.
6. **Name.** The repository stays `project-white-whale`; the resource and CLI
   are called `whale` here. The site is live as "Project White Whale"; a
   descriptive public name (e.g. "Pacific oyster methylome atlas") is still
   open.
