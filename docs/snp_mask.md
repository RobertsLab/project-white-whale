# The C/T SNP mask

In bisulfite data a C→T SNP is indistinguishable from an unmethylated C, and
the oyster genome is highly heterozygous. Before methylation is compared,
CpGs on SNPs are masked (PLAN.md §6.1). The SNPs are called from the WGBS
libraries themselves, pooled, with `workflows/snp_mask`.

## How it calls SNPs

```
SORT_BAM (each library) -> PILEUP (BISCUIT, every library, per 5 Mb chunk) -> BUILD_MASK
```

- **Input**: each library's deduplicated Bismark BAM. The methylation workflow
  keeps them with `--publish_bams true` (under `<WHALE_METHYLOMES>/<STUDY_ID>/bam/`,
  ~3-4 GB a library); for runs made without it they are in the Nextflow work
  directory until it is cleaned.
- **PILEUP**: `biscuit pileup` (BISCUIT 1.10.2) over all libraries at once
  (MAPQ >= 40, base quality >= 20). BISCUIT reads Bismark's strand tags and
  separates reads that inform a SNP from reads where bisulfite conversion makes
  C/T (or G/A) ambiguous. `pileup_sites.awk` keeps only BISCUIT's definite
  base calls, sums them over the libraries, and passes on sites with >= 2
  alternative reads.
- **BUILD_MASK** (`snp_mask.py`): a site is a SNP when, pooled over the
  libraries, it has **>= 3 definite alternative reads** and they are **>= 2%** of
  its definite reads (`--snp_min_alt`, `--snp_min_af`; changing them re-runs
  only this step with `-resume`). A reference CpG is masked when a SNP sits on
  its C or its G.

Pooling favours common variants: with ~5x per library, a SNP carried by one
heterozygous animal in 24 gives ~1-3 definite reads and is often below the
threshold. The mask removes the systematic bias of population SNPs; it does
not remove every private variant.

## Applying it

`nextflow run workflows/methylation ... --snp_mask <mask folder>` (or
`WHALE_SNP_MASK=<mask folder>`) leaves the masked CpGs out of every gene table;
`run_manifest.json` and `study_summary.yaml` record the mask id and checksum.
Re-running a finished study with `-resume` redoes only the gene step.

## Outputs

Published to `<WHALE_REFERENCE>/<ASSEMBLY>/snp_mask/<MASK_ID>/`:

| File | |
|---|---|
| `cpg_mask.bed.gz` | masked CpGs, BED6: chrom, C position (0-based), +2, first SNP (`1234C>T`), its alt reads, `.` — the same coordinates as the CpG files' start |
| `snps.tsv.gz` | every SNP: chrom, pos0, ref, alt, depth, alt reads, alt fraction, libraries with alt reads |
| `mask_manifest.json` | mask id, libraries and studies, thresholds, BISCUIT version, pipeline commit, genome and output SHA-256s, counts per sequence, alt-read distribution |

`whale record-mask --dir <that folder>` uploads the three files to the data
host (`reference/<ASSEMBLY>/snp_mask/<MASK_ID>/`) and adds the mask, with its
checksums, URLs, studies, thresholds and counts, to
`data/reference/snp_masks.yaml`, which is committed in a normal PR. Analyses
name the mask they used by its id. Masks are versioned, not replaced: a mask
built from more studies gets a new id.

## Running it on klone

List the libraries, one line each (`sample_id`, `study_id`, `bam`). With
`--publish_bams true` the BAMs are in `<WHALE_METHYLOMES>/<STUDY_ID>/bam/`;
for an earlier run, take them from its work directory and check the count
(a retried task can leave a partial BAM):

```bash
printf 'sample_id\tstudy_id\tbam\n' > bams.tsv
find /gscratch/srlab/sr320/whale/work -name 'SRR*.dedup.bam' -size +100M | sort | while read -r f; do
  printf '%s\tROBERTS2022_PH_PLOIDY_WGBS\t%s\n' "$(basename "$f" .dedup.bam)" "$f"
done >> bams.tsv
```

Then, from a head job set up like the methylation one (non-preemptible
partition, `whale` conda environment, the `WHALE_*` variables):

```bash
nextflow run workflows/snp_mask -profile klone -resume \
    --mask_id GCF_963853765.1.snpmask.v1 --bam_sheet bams.tsv
whale record-mask --dir $WHALE_REFERENCE/GCF_963853765.1/snp_mask/GCF_963853765.1.snpmask.v1
```

## Checks

`make real-run` (and CI) plant 20 C/T SNPs at CpGs in a copy of the fixture
reference (`tests/plant_snps.py`: 20 `TG` sites rewritten as `CG`, so the
reads' T looks exactly like an unmethylated C) and require the mask to hold at
least 80% of them and at most 5 other CpGs. On the fixtures it masks 18 of 20
and nothing else.
