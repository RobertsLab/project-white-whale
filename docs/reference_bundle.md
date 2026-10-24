# The reference bundle

Every sample is processed against one assembly, the `default: true` entry in
`data/reference/genome_assemblies.yaml` (xbMagGiga1.1, `GCF_963853765.1`,
annotation `RS_2024_06`). Each machine fetches its bundle once:

```bash
whale fetch-reference --plan      # what would be fetched, with NCBI's MD5s
whale fetch-reference             # genome + GTF + derived files (~180 MB download, ~2 min)
whale fetch-reference --rnaseq    # also transcripts and tx2gene.tsv, for Goal 2
whale fetch-reference --check     # verify an existing bundle against its provenance
```

The bundle goes to `$WHALE_REFERENCE/GCF_963853765.1/` (default
`data/reference/GCF_963853765.1/`, not committed), which is where the
workflow looks for it.

## Contents

| File | Origin |
|---|---|
| `genomic.fna` | NCBI `*_genomic.fna.gz`, decompressed (Bismark) |
| `genomic.gtf` | NCBI `*_genomic.gtf.gz`, decompressed |
| `rna.fna`, `tx2gene.tsv` | with `--rnaseq`: NCBI `*_rna.fna.gz`; transcript → GeneID from the GTF |
| `chrom.sizes` | sequence lengths, from the genome |
| `features/genes.bed` | GTF `gene` features; name = NCBI GeneID, score column = biotype |
| `features/exons.bed` | every transcript's exons, de-duplicated |
| `features/introns.bed` | gaps between consecutive exons of each transcript |
| `features/promoters.bed` | 1 kb upstream to 200 bp downstream of each gene's TSS, strand-aware, clipped |
| `features/intergenic.bed` | sequence covered by no gene |
| `md5checksums.txt`, `*.gz` | NCBI's checksum list and the originals, kept so re-runs download nothing |
| `reference_provenance.json` | URLs, NCBI MD5 and SHA-256 of every file, GTF header, feature rules and counts |

BED files are 0-based, half-open, BED6, sorted in genome order.

On 2026-10-23 the bundle had 29 sequences (564.0 Mb), 33,331 genes, 284,442
exons, 249,701 introns and 27,908 intergenic intervals.

## Checks

- Every download is compared with NCBI's published MD5 before use; a mismatch
  deletes the file and stops.
- The GTF's `annotation-source` must name the registered annotation release.
  NCBI can re-annotate an assembly in place; if it has, the command exits
  non-zero, and `genome_assemblies.yaml` must be updated (in its own PR)
  before the bundle is used.
- `--check` recomputes every file's SHA-256 against the provenance file.

## The Bismark index

The workflow builds the Bismark index the first time it runs on a machine and
keeps it in the bundle as `bismark_genome/` (with its own copy of the genome,
as `genome.fa`, the name Bismark accepts); every later run, for any study,
reuses it rather than spending hours rebuilding it. Set `WHALE_BISMARK_INDEX`
to keep it somewhere else (it then lives in `$WHALE_BISMARK_INDEX/bismark_genome/`).
On klone it must be on storage the compute nodes see, like the bundle itself.

`bismark_genome/genome.sha256` records which genome was indexed. If the
bundle's `genomic.fna` no longer matches it (the bundle was fetched again for
a new release), the workflow stops and asks for `bismark_genome/` to be
deleted, so nothing is aligned to an index of a different genome. Start the
first study on a new machine alone: two runs building the index at the same
moment would both write it.

`bismark_genome/` is not listed in `reference_provenance.json`, so
`whale fetch-reference --check` ignores it.

Not in the bundle yet: the CpG index and the C/T SNP mask (later Phase 1,
PLAN.md §6.1).
