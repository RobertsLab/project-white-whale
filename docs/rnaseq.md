# RNA-seq: quantifying a study's paired expression data

Goal 2 (PLAN.md §7.2) relates each study's methylation to gene expression.
For every study whose registration names RNA-seq under
`paired_rnaseq.bioproject`, `workflows/rnaseq` quantifies those libraries
against the same reference bundle, with the same gene IDs, as the methylomes.

```
FASTQC
FASTP ─► SALMON_QUANT ─► GENE_MATRICES ─► RNASEQ_QC ─► MULTIQC ─► MANIFEST
            ▲
SALMON_INDEX (rna.fna + genomic.fna as decoys; built once per assembly, cached)
```

## Running a study

1. **Reference.** The bundle needs the transcript files, once per machine:
   `whale fetch-reference --rnaseq` adds `rna.fna` (NCBI's transcripts) and
   `tx2gene.tsv` (transcript → NCBI GeneID, from the GTF).
2. **Sample sheet.** From ENA, never from the paper, like the methylation sheet:

   ```bash
   whale fetch-rnaseq-samplesheet --study WANG2020_TIDAL_HEAT_WGBS --show-attributes
   ```

   ```bash
   whale fetch-rnaseq-samplesheet --study WANG2020_TIDAL_HEAT_WGBS --condition-attribute treatment
   ```

   This writes `data/studies/<STUDY_ID>/rnaseq/samplesheet.tsv`,
   `fastq_manifest.tsv` and `ena_provenance.json`, keeping the BioProject's
   `RNA-Seq` runs (`--library-strategy` to change). The sheet has the
   methylation sheet's columns (`sample_id` = run accession, `run_accession`,
   `fastq_1`, `fastq_2`, `library_id`, `condition`, the attributes). It warns
   when `paired_rnaseq.relationship` is not `same_animals`: such libraries can
   only be compared with the methylomes by group, not one to one.
3. **FASTQ.** Stage them under `$WHALE_RAW/<STUDY_ID>/rnaseq/`:

   ```bash
   whale fetch-fastq --study WANG2020_TIDAL_HEAT_WGBS --manifest data/studies/WANG2020_TIDAL_HEAT_WGBS/rnaseq/fastq_manifest.tsv --dest $WHALE_RAW/WANG2020_TIDAL_HEAT_WGBS/rnaseq --accession PRJNA407831
   ```

   Leave out `--mark-staged` and `--launch`: those move the methylation stage.
4. **Run**, with the study's config and the machine's profile:

   ```bash
   nextflow run workflows/rnaseq -c config/WANG2020_TIDAL_HEAT_WGBS.config -profile raven
   ```

   RNA-seq settings (`tag_seq`, `se_adapter`, the QC thresholds) go in the
   same `config/<STUDY_ID>.config`, the one config a study PR may change; the
   methylation workflow ignores them.

When the run succeeds it calls `whale record-rnaseq` (below). Turn that off
with `--record_status false`.

### Which branch

The run records itself on the study's own branch
(`study/<ACCESSION>/<STUDY_ID>`) while the study is claimed, so the RNA-seq
lands in the same PR. For a study whose methylation has already merged, start
a branch for it:

```bash
git switch -c rnaseq/WANG2020_TIDAL_HEAT_WGBS origin/main
```

`whale record-rnaseq` refuses to commit on any other branch, and never moves
`status.yaml`: its stages are the methylation workflow's.

## Design

**Quantifier: Salmon 1.10.3, selective alignment.** Salmon estimates
transcript abundance directly, models fragment GC bias and handles
multi-mapping reads by EM, and its quantifications are what tximport and
DESeq2 expect. A genome aligner (STAR, HISAT2) would add splice junctions and
intron coverage, which Goal 2's splicing analysis will need (below), at
several times the compute; it can be added as a second branch later without
changing these outputs. Salmon 1.10.3 rather than 2.x: it is the version the
literature and nf-core use, and its index and output formats are stable.

**Decoys: the whole genome.** The index is built from `rna.fna` plus
`genomic.fna` marked as decoy sequence (Srivastava et al. 2020, *Genome
Biology* 21:239). A fragment that matches the genome better than any
transcript, typically from an intron, an unannotated locus or genomic DNA, is
discarded instead of being forced onto a similar transcript. This matters in
oyster: the genome is highly heterozygous, the annotation is incomplete, and
many libraries carry intronic reads. The fixture run shows what decoys do
and do not do: unspliced fragments across an exon-intron boundary are partly
taken by the decoys (about a third paired-end, a sixth single-end, on 75 bp
reads); the rest still pass Salmon's alignment score with their exonic part
and count towards the gene. `--salmon_decoys false` builds a transcript-only
index (`salmon_index_nodecoy/`). The index of the ~650 Mb genome needs
roughly 20-30 GB of memory while it builds (the raven and klone profiles give
`SALMON_INDEX` 64 GB) and is cached in the reference bundle (`storeDir`), so
it is built once per machine. Its `whale_index.json` records the SHA-256 of
`rna.fna` and `genomic.fna`, k and the Salmon version; `main.nf` stops if the
bundle's files or this run's k / index arguments no longer match.

**Trimming: fastp**, vendored from AREE (VENDORED_FROM_AREE.md) with AREE's
lesson: single-end adapters are given explicitly (`se_adapter`, Illumina
TruSeq read 1 by default), because fastp's single-end auto-detection once
picked an abundant oyster transcript as the "adapter". Paired-end detection
uses read overlap and is kept. `--tag_seq` (3' tag libraries such as QuantSeq)
adds poly-X trimming and turns off Salmon's length correction and GC bias.
Trimmed reads are not published.

**Library type** is detected per library (`-l A`). Salmon's auto-detection
falls back to an unstranded type on small libraries, so `qc_summary.tsv`
reports `strandedness` (reverse / forward / unstranded) and `sense_pct` from
Salmon's per-orientation fragment counts as well; `rnaseq_summary.yaml` flags
a study whose libraries differ.

**Gene aggregation: tximport's rule, in Python.** `GENE_MATRICES`
(`workflows/rnaseq/bin/rnaseq_quant.py`) sums each gene's transcripts' TPM
and estimated counts, and takes the TPM-weighted mean of their effective
lengths as the gene length (the plain mean where the gene is not expressed in
that sample), which is what `tximport(..., countsFromAbundance = "no")`
returns. Downstream, DESeq2 can use the counts with the length matrix as an
offset (`DESeqDataSetFromTximport`), or a study can be re-imported with
tximport from the published `quant.sf` files. Gene IDs are `tx2gene.tsv`'s,
NCBI GeneIDs, the same as `features/genes.bed` and the methylation gene
matrices (`data/studies/<ID>/genes/`), so the two join on `gene_id`.
Transcripts missing from `tx2gene.tsv` and genes with no transcript in the
index are counted in `matrices/matrices_report.json` and the run manifest.

## Outputs

Under `$WHALE_METHYLOMES/<STUDY_ID>/rnaseq/` on the machine that ran it:

| Path | What |
|---|---|
| `salmon/<sample>/` | `quant.sf`, `aux_info/meta_info.json`, `lib_format_counts.json`, `cmd_info.json`, `logs/` |
| `matrices/gene_{tpm,counts,length}.tsv.gz` | gene × sample, every library in sample-sheet order |
| `matrices/transcript_{tpm,counts}.tsv.gz` | transcript × sample |
| `qc/` | FastQC and fastp reports |
| `qc_summary.tsv` | per library: reads, trimming, Salmon mapping and decoy rates, library type, strandedness, genes with TPM ≥ 1, `qc_pass`, `qc_reasons` |
| `multiqc_report.html` | for people |
| `run_manifest.json` | parameters, tool versions, the Salmon index description, input checksums |

### Recording (`whale record-rnaseq`)

Committed to `data/studies/<STUDY_ID>/rnaseq/`: the three gene matrices,
`qc_summary.tsv`, `run_manifest.json`, `rnaseq_summary.yaml` (QC medians,
strandedness, genes detected in every passing library) and
`rnaseq_manifest.tsv`. Uploaded to the data host (`data/data_host.yaml`) under
`rnaseq/<STUDY_ID>/<sample_id>/` (each library's Salmon files and reports)
and `rnaseq/<STUDY_ID>/_study/` (transcript matrices, MultiQC report);
`rnaseq_manifest.tsv` lists every one with its size, SHA-256 and URL. The
matrices keep every library; filter on `qc_summary.tsv`'s `qc_pass`. Running
it again is safe and commits nothing when nothing changed; a failed upload is
finished by running it again.

### QC gate

| Parameter | Default | Test profile |
|---|---|---|
| `rnaseq_min_mapping` | 30 % of fragments assigned to transcripts | 30 |
| `rnaseq_min_mapped` | 1,000,000 fragments | 1,000 |
| `rnaseq_min_genes` | 1,000 genes with TPM ≥ 1 | 3 |

First guesses: with the genome as decoy, intronic and intergenic reads do not
count as mapped. Tune them once the first studies have run.

## Testing

```bash
make rnaseq-stub        # the DAG, -stub
make rnaseq-real-run    # every tool in Docker, then tests/check_rnaseq_run.py
```

The fixtures (`tests/fixtures/make_rnaseq_fixtures.py`) build five
transcripts of the four fixture genes from the fixture genome's exons (G1
with two isoforms), and a stranded (dUTP) library per sample, R1 paired-end
and R2 single-end, with known fragments per transcript plus intronic,
intergenic and unspliced exon-intron fragments. The check requires every
gene's estimate between its simulated count and that plus its own unspliced
fragments, no expression for the gene with no reads, isoforms summing to their
gene, fragments discarded by the decoys, reverse strandedness and both
libraries through the gate. CI runs both (`.github/workflows/ci.yml`, job
`rnaseq`).

## Not yet provided (for the Goal 2 analyses)

- **Differential expression.** No DESeq2 step; PLAN.md §6.3 keeps contrasts
  downstream so every study shares one model. The counts and length matrices
  are its input.
- **Splicing.** Exon methylation vs exon inclusion (PLAN.md §7.2.3) needs PSI
  from genome alignments (STAR + rMATS / SUPPA from Salmon transcript TPM is a
  weaker proxy). Not here.
- **Pairing table.** Which RNA-seq library belongs to which methylation
  library (same animal or pool) is not recorded anywhere yet; one-to-one
  analyses need a curated `rnaseq/pairs.tsv` per `same_animals` study.
- **Runs per library.** Libraries sequenced on several runs are quantified per
  run, as in the methylation workflow; summing them (or passing all runs to
  one `salmon quant`) is left for when a study needs it.
- **Pooled matrices and dashboard.** `whale build-atlas` and the dashboard do
  not read `rnaseq_summary.yaml` or the gene matrices yet.
