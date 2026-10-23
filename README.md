# White Whale

**Public DNA methylation of the Pacific oyster (*Crassostrea / Magallana gigas*),
reprocessed from raw reads on one reference, so studies can be compared.**

Two goals:

1. **Characterize DNA methylation across the genome, across studies**: what is
   shared and what varies by tissue, life stage, population, ploidy and
   treatment.
2. **Determine the role of DNA methylation in regulating gene expression**,
   using every study with paired RNA-seq.

The design, analyses and phases are in [PLAN.md](PLAN.md). The structure
follows [AREE](https://github.com/sr320/AREE), narrowed to one species and one
assay family.

> **Status: Phase 1 under way.** The registry, validation, the claim system
> for running studies on several machines, the compute profiles, a
> stub-tested methylation workflow, the ENA sample-sheet generator and the
> reference-bundle fetcher are in place. No study has been registered or run yet. See [PLAN.md §9](PLAN.md#9-phases).

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
whale --help
```

Requires Python ≥ 3.10. Running the workflow also needs Nextflow ≥ 24.10 and
Docker or Apptainer.

## How studies get processed

Every public study is one entry in [`registry/backlog.tsv`](registry/backlog.tsv).
A machine **checks a study out** by claiming it: `whale start-study` creates a
branch and worktree for the study and pushes the branch. The push is the claim,
so two machines cannot take the same study. Progress is recorded on the branch
and every machine can see it.

```bash
whale backlog                       # what is free, claimed (where, how far) or registered
whale start-study --accession PRJNA682817 --study-id ROBERTS2022_PH_PLOIDY_WGBS --profile raven
cd ../white-whale-ROBERTS2022_PH_PLOIDY_WGBS
whale fetch-samplesheet --bioproject PRJNA682817 --condition-attribute ploidy --condition-attribute water_pH
# ...fill the study YAML, stage FASTQ...
whale set-status --stage fastq_staged --note "105 GB, MD5 ok"
nextflow run workflows/methylation -c config/ROBERTS2022_PH_PLOIDY_WGBS.config -profile raven
whale set-status --stage called
whale claims                        # every claimed study, its machine and stage
```

The study's PR changes only that study's files, so study PRs merge in any order.
Details: [running studies on several machines](docs/running_studies_on_several_machines.md),
[compute profiles (desktop, raven, klone)](docs/compute_profiles.md),
[adding a study](docs/adding_a_study.md),
[the reference bundle](docs/reference_bundle.md).

## Commands

| Command | Purpose |
|---|---|
| `whale backlog` | The backlog with each study's state: free, registered, or claimed (machine and stage) |
| `whale start-study --accession <PRJ…> --study-id <ID> --profile <desktop\|raven\|klone> [--part <P>]` | Check a study (or one part of a registered study) out to this machine |
| `whale set-status --stage <stage> [--note …] [--run-name …]` | Record progress on the current study branch; commits and pushes |
| `whale claims` | Every claimed study, from the pushed branches |
| `whale release-study --branch <b> --reason …` | Give a claim up so another machine can take it |
| `whale fetch-samplesheet --bioproject <PRJ…> --condition-attribute <a> […] [--show-attributes]` | Sample sheet, FASTQ manifest (ENA MD5s) and provenance from the deposited metadata; downloads nothing |
| `whale fetch-reference [--rnaseq] [--plan] [--check]` | Fetch, MD5-verify and derive the reference bundle (genome, GTF, chrom.sizes, feature BEDs) |
| `whale validate-study [files] [--strict]` | Validate study YAMLs against the schema and vocabularies |
| `whale validate-status [files]` | Validate `status.yaml` files |
| `whale build-registry` | Rebuild `registry/study_registry.csv` (not committed) |
| `whale check-study-scope --branch <b>` | Fail if a study branch touches files it does not own (run by CI) |
| `whale profiles` | The compute profiles |

## Layout

```
registry/        backlog.tsv, studies/<STUDY_ID>.yaml, controlled_vocabularies/
schemas/         study and status JSON schemas
workflows/       methylation/ (Nextflow DSL2)
modules/         Nextflow processes
config/          base.config (profiles) and one config per study
data/            reference/ (assembly registry), studies/<STUDY_ID>/ (small committed outputs);
                 raw/, methylomes/ and reference bundles are per-machine symlinks
src/whale/       the `whale` CLI
tests/           pytest suite and workflow fixtures
docs/            how-tos
```

## The previous White Whale

Until October 2026 this repository was a catalogue of public *C. gigas*
RNA-seq and methylation datasets with a downloader. That catalogue is the seed
of `registry/backlog.tsv`; the earlier files are in the git history (commit
`c7c64dd`).

## License

MIT. See [LICENSE](LICENSE).
