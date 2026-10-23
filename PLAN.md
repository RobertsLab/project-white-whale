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

```
project-white-whale/
├── README.md, PLAN.md, LICENSE, CITATION.cff, CONTRIBUTING.md
├── pyproject.toml, Makefile, environment.yml
├── schemas/                      study.schema.json, sample.schema.json, methylome.schema.json
├── registry/
│   ├── studies/<STUDY_ID>.yaml   one per BioProject (or per coherent sub-project)
│   ├── controlled_vocabularies/  assay_chemistry, tissue, life_stage, stressor, ploidy, sex, …
│   ├── claims/                   (generated) who holds which study, from git branches
│   └── study_registry.csv        (generated) index
├── workflows/
│   ├── methylation/              Nextflow DSL2: FASTQ → per-sample CpG calls → QC → summaries
│   └── rnaseq/                   only for paired-expression studies (salmon → gene/transcript TPM)
├── modules/                      Nextflow processes (vendored from AREE where identical)
├── config/
│   ├── base.config               profiles: desktop | raven | klone (+ docker/apptainer/conda)
│   └── <STUDY_ID>.config
├── data/
│   ├── reference/GCF_963853765.1/   genome, GTF, feature BEDs, CpG index, C/T SNP mask (not in git)
│   ├── raw/                      → symlink to this machine's FASTQ drive
│   ├── studies/<STUDY_ID>/       samplesheet.tsv, fastq_manifest.tsv, ena_provenance.json,
│   │                             qc_summary.tsv, methylome_manifest.tsv (committed, small)
│   └── methylomes/               → symlink; per-sample calls land here (not in git)
├── src/whale/                    the `whale` CLI and library
│   ├── intake/                   schema validation, ENA sample sheets, claims
│   ├── reference/                bundle fetch, feature annotation, CpG index, SNP mask
│   ├── qc/                       pipeline QC gates
│   ├── atlas/                    Goal 1: matrices, variance partitioning, core/variable methylome
│   ├── regulation/               Goal 2: methylation × expression models
│   ├── reporting/                study pages, dashboard JSON
│   └── query/                    (later) MCP server over the summarized tables
├── analysis/                     Quarto notebooks for Goals 1 and 2 (rendered into docs/)
├── docs/                         Quarto site: dashboard, study pages, atlas, regulation results, how-tos
├── tests/
└── .github/workflows/            ci.yml (tests, schema validation, scope check on study PRs), pages.yml
```

---

## 4. Layer 1 — Study registry

### 4.1 Study YAML

Adopt AREE's study schema and trim it to this project. Every study records:
`study_id`, accessions (BioProject / SRA / ENA / GEO / DOI / PMID), citation,
species name as deposited, population / line / ploidy, assay chemistry
(`wgbs`, `emseq`, `rrbs`, `medip`, `mbd`), library layout, platform,
`paired_rnaseq` (BioProject + relationship: `same_animals` /
`same_design_different_animals` / `none`), per-sample attributes
(tissue, life stage, sex, ploidy, treatment, timepoint, pooling), data
availability, QC and processing status, quality flags, limitations, provenance.

Two vocabularies matter more here than in AREE and are added:

- `assay_chemistry` with a **quantitative-comparability tier**: base-resolution
  (WGBS, EM-seq, RRBS-with-coverage-caveat) vs enrichment (MeDIP, MBD).
  Enrichment data is processed and shown but never enters a base-resolution
  comparison.
- `sample_unit`: individual / pool-of-N / tissue-pooled, because pooling
  changes what variance means.

### 4.2 Seed list

The verified catalogue from the previous repository plus the studies AREE has
already registered become the first registrations. Where AREE has a YAML
(`WANG2020_TIDAL_HEAT_WGBS`, `VALDIVIESO2025_POMS_RS_EMSEQ`,
`DECICOMP2025_POMS_AGE_EMSEQ`, `CIBNOR2021_HEAT_RRSS`), copy it and keep its
`study_id` so the two projects can cross-reference.

| Study (working ID) | BioProject | Chemistry | Runs | ~Size | Paired RNA-seq | Goal |
|---|---|---|---|---|---|---|
| POMS adaptation, Ifremer | PRJEB60400 | WGBS | 246 | 450–500 GB | check | 1 |
| GESTINOV POMS gill/mantle (`VALDIVIESO2025_POMS_RS_EMSEQ`) | PRJEB81880 | EM-seq | 40 | 400–450 GB | same design (AREE `DECICOMP2025_POMS_AGE` RNA-seq is a sister) | 1, 2 |
| DECICOMP aging (`DECICOMP2025_POMS_AGE_EMSEQ`) | PRJEB105019 | WGBS | 60 | 450–500 GB | yes, same project | 1, 2 |
| PESTO methylseq | PRJEB58545 | WGBS | 48 | 350–400 GB | check | 1 |
| Transgenerational infection resistance | PRJNA609264 | WGBS | 47 | 600–650 GB | check | 1 |
| pH × ploidy ctenidia, Roberts Lab | PRJNA682817 | WGBS | 24 | 100–110 GB | yes, same animals | 1, 2 |
| Tidal origin × heat (`WANG2020_TIDAL_HEAT_WGBS`) | PRJNA562805 | WGBS | 12 | 141 GB | yes, same animals (PRJNA407831) | 1, 2 |
| Epigenomics series | PRJNA807732, PRJNA213124 | WGBS | 38 | ~190 GB | check | 1 |
| Heat RRBS (`CIBNOR2021_HEAT_RRSS`) | PRJNA690951 | RRBS | — | — | — | 1 (RRBS tier) |
| First oyster methylomes, Wang 2014 | PRJNA173440 | WGBS | 4 | small | developmental RNA-seq PRJNA146329 | 1, 2 |
| Developmental methylome, Riviere 2017 | PRJNA324546 | MeDIP | 21 | 5–10 GB | yes (published) | enrichment tier |
| Sperm methylome + RNA-seq, Olson 2014 | to verify | WGBS | — | — | same animals | 2 |
| Gill MBD-seq + RNA-seq, Gavery 2013 | to verify | MBD | — | — | yes | enrichment tier |

Plus candidates to verify against SRA before registering: Li 2024 *Vibrio*
WGBS, Sun 2022 gonad sex, Tan 2022 growth/shell colour, Wang 2023
transgenerational tidal, Song 2017 exon methylation. The first task of Phase 1
is a fresh `"Crassostrea gigas"[Organism] AND ("Bisulfite-Seq" OR "MeDIP-Seq"
OR "MBD-Seq")[Strategy]` sweep to make this list complete and to find which
BioProjects also carry RNA-seq.

Roughly **3 TB of FASTQ, ~550 runs**, of which the base-resolution set is
~500 runs.

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
- Big studies (PRJEB60400, 246 runs) can be claimed as **parts**
  (`--part gill`, `--part 1of4`) that share a study YAML; each part is its own
  branch and status file.

### 4.4 Compute profiles

| Profile | Machine | Executor | Scratch | Typical job |
|---|---|---|---|---|
| `desktop` | lab Macs / workstations | local | `$WHALE_WORKDIR` on a fast internal/SSD drive | small studies (≤ 150 GB), QC, summaries, Quarto |
| `raven` | Roberts Lab server | local, high `maxForks` | raven's local scratch | medium studies (150–500 GB), whole-study Bismark without chunking |
| `klone` | UW Hyak | Slurm `ckpt-all`, 4 h 45 m tasks, chunked alignment, Apptainer | `/scr` per node; `/gscratch/scrubbed/$USER/whale` | the deep studies (PRJEB60400, PRJEB81880, PRJEB105019, PRJNA609264) |

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
that ran them, then are synced to the lab data host (gannet / owl) under
`whale/methylomes/<STUDY_ID>/<sample>/`; `methylome_manifest.tsv` (committed)
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
| `gene_features` | ~30 k genes × feature class × samples | tens of MB | **committed** |
| `sample_qc` | samples × QC metrics + flags | small | **committed** |
| `sample_metadata` | samples × harmonized attributes (tissue, stage, sex, ploidy, population, treatment, chemistry, pool size, study) | small | **committed** |

Every row carries `study_id`, `chemistry`, and the manifest hash, so a table
can always be filtered to "WGBS only" or "individual animals only".

### 6.3 Re-derived differential results

For every study with a registered within-study contrast, run one shared
model (DSS or methylKit with the same cutoffs, same tile definition, same
minimum coverage) from the harmonized calls, and store the result next to the
authors' published DML/DMR list when one exists. This is what makes "did it
replicate" answerable, and it feeds Goal 2's DML-vs-DEG tests.

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

Analyses are Quarto documents that read only the committed small tables and
the data-host Parquet, so anyone can re-render them without the FASTQ.

---

## 8. Layer 5 — Outputs

- **Quarto site** (GitHub Pages), regenerated on every push to `main`:
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
- **MCP server** (`whale mcp`, later): read-only tools over the summarized
  tables, as AREE's, so Claude can answer "which genes are hypermethylated in
  gill across every study" with study IDs cited.

---

## 9. Phases

| Phase | Deliverable | Done when |
|---|---|---|
| **0. Scaffold** (weeks 1–2) | Repo layout, `whale` CLI skeleton, schemas, vocabularies, `start-study` / claims, CI with scope check, profiles, docs stubs | `whale validate-study` and `whale start-study --no-push` work on a stub; CI green |
| **1. Registry** (weeks 2–3) | Fresh SRA/ENA sweep; every base-resolution BioProject registered with ENA-derived sample sheets and FASTQ manifests; paired RNA-seq identified; AREE YAMLs imported | Dashboard shows the full backlog with sizes and which have RNA-seq |
| **2. Reference + pipeline** (weeks 3–6) | Reference bundle, feature BEDs, CpG index; pipeline runs end to end on the smallest WGBS set (PRJNA173440, 4 runs) on a desktop, then PRJNA682817 on raven; first SNP mask | Two studies `qc_passed`, methylomes on the data host, manifests committed |
| **3. Klone** (weeks 5–8) | Klone profile proven on one deep study (PRJEB81880, 40 EM-seq runs); chunked alignment; download pattern | Study merged with status `summarized` |
| **4. Backlog** (weeks 6–16, parallel across machines) | Every registered base-resolution study claimed and processed; enrichment-tier studies processed with a MeDIP/MBD branch | All studies `summarized` or explicitly `excluded` with a reason |
| **5. Goal 1 analyses** (from week 8, as studies land) | Atlas matrices, variance partitioning, core/variable classification, track hub, atlas pages | Pages render from committed tables; findings carry study counts |
| **6. Goal 2 analyses** (from week 10) | RNA-seq workflow on paired studies; the five regulation analyses; pooled estimates | Regulation pages render; each claim lists its supporting studies |
| **7. Write-up** | Methods page, implementation-status table, manuscript draft from the site | — |

Phases 4–6 overlap by design: analyses are re-rendered as each study merges.

---

## 10. Decisions to make before Phase 0

1. **Reuse of AREE code: vendor or depend?** Vendoring the five Nextflow
   modules and the ENA sample-sheet code is simplest and keeps this repo
   self-contained; depending on AREE as a git submodule keeps fixes in one
   place. Recommendation: vendor, with a `VENDORED_FROM_AREE.md` noting the
   commit, and revisit if both projects keep changing the same modules.
2. **Data host for per-sample methylomes.** gannet / owl with HTTPS URLs in
   the manifests is the lab's established pattern; confirm space for ~1 TB of
   calls + bigWigs.
3. **Raven specifics** (cores, RAM, scratch path, container engine) to size
   the `raven` profile.
4. **Study granularity for PRJEB60400** (246 runs): one study claimed in
   parts, or split by sub-experiment once its sample metadata is read.
5. **Name.** The repository stays `project-white-whale`; the resource and CLI
   are called `whale` here. A public-facing name (e.g. "Pacific oyster
   methylome atlas") can be chosen when the site goes live.
