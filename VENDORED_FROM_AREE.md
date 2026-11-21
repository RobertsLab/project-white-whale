# Code adapted from AREE

White Whale vendors code from [AREE](https://github.com/sr320/AREE) (MIT)
rather than depending on it (PLAN.md §10). Each entry names the AREE source
and the AREE commit it was taken from, so a fix in either project can be
carried across by hand.

| White Whale | AREE source | AREE commit | Changes |
|---|---|---|---|
| `src/whale/claims.py` | `src/intake/study_claims.py` | `44022f1` | BioProject-only accessions; per-branch `status.yaml`; compute profiles; part claims; `set_status`, `list_claims`, `release_study` |
| `src/whale/ena.py` (`whale fetch-samplesheet`) | `src/intake/ena_samplesheet.py` | `0a87b69` | Methylation strategies kept by default, others reported; taxid 29159 enforced; sheet in the methylation workflow's columns with `sample_id` = run accession and `library_id` = experiment; no methylKit sheet; `--show-attributes`; BioSample title as `sample_title` |
| `src/whale/reference.py` (`whale fetch-reference`) | `src/intake/reference_bundle.py` | `eb838d7` | Genome and GTF always, transcripts only with `--rnaseq`; adds `chrom.sizes` and gene/exon/intron/promoter/intergenic BEDs; no proteome or crosswalk check; `$WHALE_REFERENCE` |
| `src/whale/scope.py` | `src/validation/study_scope.py` | `44022f1` | Part ownership; test-file naming |
| `.github/PULL_REQUEST_TEMPLATE/study.md` | same path | `44022f1` | Methylome checklist |
| `config/base.config` (`klone` profile) | `config/base.config` | `44022f1` | `WHALE_*` variables; methylation processes only |
| `modules/methylation/*.nf` (containers, Bismark `--score_min L,0,-0.6`) | `modules/methylation/` | `44022f1` | Rewritten; container pins and Bismark settings kept |
| `modules/rnaseq/fastp.nf` | `modules/rnaseq/trim_fastp.nf` | `8b088a0` | Container pin, explicit single-end adapter, `--detect_adapter_for_pe`, `--tag_seq` poly-X kept; trimmed reads not published; outputs under `<STUDY_ID>/rnaseq/qc/` |
| `modules/rnaseq/salmon.nf` (`SALMON_QUANT`, `SALMON_INDEX`) | `modules/rnaseq/align_quant_salmon.nf` | `8b088a0` | Biocontainers pin of Salmon 1.10.3; genome decoys added to the index, which is cached with `storeDir` and checksummed; no `--geneMap` (genes from `tx2gene.tsv` in GENE_MATRICES); PE-only `--gcBias` and `--tag_seq` handling kept |
| `workflows/rnaseq/nextflow.config` (`se_adapter`, `tag_seq`) | `workflows/rnaseq/nextflow.config` | `8b088a0` | Same defaults |
| `data/reference/genome_assemblies.yaml` | same path | `44022f1` | Two assemblies; `default` flag |
| `docs/compute_profiles.md` (klone section) | `docs/running_on_klone.md` | `44022f1` | Summarized |

When vendoring more, add a row here with the AREE commit
(`git -C <AREE checkout> log -1 --format=%h -- <file>`).
