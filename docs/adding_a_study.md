# Adding a study

1. **Pick it from the backlog.** `whale backlog` lists every known public
   *C. gigas* methylation BioProject and whether it is free. A study not in
   the backlog: add it to `registry/backlog.tsv` in a small PR to `main` first.
2. **Name it.** `FIRSTAUTHOR<YEAR>_<TOPIC>[_<CHEMISTRY>]`, e.g.
   `WANG2020_TIDAL_HEAT_WGBS`. If AREE has registered the same deposit, reuse
   its `study_id` and set `aree_study_id`.
3. **Claim it** for the machine that will run it (see
   [running studies on several machines](running_studies_on_several_machines.md)):

   ```bash
   whale start-study --accession <PRJ…> --study-id <STUDY_ID> --profile <desktop|raven|klone>
   ```

4. **Fill the YAML** in the new worktree, from the deposit first (ENA/SRA run
   and BioSample metadata) and the paper second. Leave anything neither states
   as `null`. Record pooling, shared BioSamples, confounded designs and missing
   spike-ins as quality flags. Iterate with

   ```bash
   whale validate-study registry/studies/<STUDY_ID>.yaml
   ```

   until it is clean, then `--strict`.
5. **Generate the sample sheet and FASTQ manifest** from ENA, never by hand.
   First see what the submitters deposited:

   ```bash
   whale fetch-samplesheet --bioproject <PRJ…> --show-attributes
   ```

   then name the attributes that define the groups, in order:

   ```bash
   whale fetch-samplesheet --bioproject PRJNA682817 \
     --condition-attribute ploidy --condition-attribute water_pH \
     --attribute replicate --attribute tissue
   ```

   This writes `samplesheet.tsv` (the columns the workflow reads:
   `sample_id`, `run_accession`, `fastq_1`, `fastq_2`, plus the design),
   `fastq_manifest.tsv` (URLs with ENA's MD5s) and `ena_provenance.json` into
   `data/studies/<STUDY_ID>/`; inside a study worktree `--study` is taken from
   the branch. It keeps only methylation runs (Bisulfite-Seq, MeDIP-Seq,
   MBD-Seq, MRE-Seq) unless `--library-strategy` says otherwise, refuses runs
   that are not taxid 29159, and reports, without keeping, any other runs in
   the project: often the paired RNA-seq, which goes under `paired_rnaseq`.
   Its summary prints the numbers for `samples.n_runs`, `samples.n_libraries`
   and `data_availability.fastq_size_gb`, and warns when a quality flag is
   probably needed (`shared_biosample`, `low_replication`).

   When the design is not in the BioSample attributes:
   - in a run-level string (`library_name`, `experiment_title`): add
     `--label-field library_name --label-pattern '<regex with named groups>'`
     and use the group names as condition attributes;
   - only in the authors' own table: commit it under `data/studies/<STUDY_ID>/`
     and add `--design-table <tsv> --design-join <attribute> --design-source <url>`;
     add the `design_from_paper` quality flag.

   With the YAML valid and these three files committed,
   `whale set-status --stage registered`.
6. **Stage FASTQ** with

   ```bash
   whale fetch-fastq --dry-run        # what is verified already, what is left, free space
   whale fetch-fastq --mark-staged    # download, verify, record, and move to fastq_staged
   ```

   It reads `data/studies/<STUDY_ID>/fastq_manifest.tsv` and writes into
   `$WHALE_RAW/<STUDY_ID>/`. Files already mirrored on gannet
   (`whale mirror-fastq`, [FASTQ mirror](fastq_mirror.md)) are copied from
   there and checked against the mirror ledger's MD5; the rest come from ENA
   and are checked against ENA's. (Copies converted from NCBI's `.sra` hold
   the same reads as ENA's files but different bytes, so the two MD5s differ.)
   A mirror copy that fails falls back to ENA. Each file downloads to
   `<name>.part` and resumes from where it stopped (re-run the same command
   after an interruption); a mismatch deletes it. Verified files get a
   `<name>.md5` sidecar so re-runs skip them, and are listed in
   `data/studies/<STUDY_ID>/fastq_staged.tsv` with their source, ENA's MD5,
   the machine and the date (committed). `--jobs` sets parallel downloads
   (default 2), `--run` limits it to some runs, `--source ena|mirror` forces
   one source. It refuses to start if free space is short, and refuses a
   second downloader into the same directory. On klone, run it as an
   `sbatch` job on the first non-preemptible partition with room, not on
   `ckpt-all` and not on a login node
   ([choosing a partition](compute_profiles.md#choosing-a-partition)).

   If a study is not mirrored yet, mirror it from raven first
   (`whale mirror-fastq --accession <PRJ…>`): ENA is slow from most machines
   ([speeds](compute_profiles.md#where-fastq-comes-from-and-how-fast)).
7. **Set `chemistry`** in `config/<STUDY_ID>.config` to match the YAML
   (`start-study` writes it as `null`; `whale validate-study` fails while a
   filled-in YAML and its config disagree) and
   run the workflow with the profile you claimed, from the study's worktree
   (`workflows/methylation`; MeDIP and MBD studies use `workflows/enrichment`,
   see [enrichment tier](enrichment.md)),
   with `whale` on the head job's `PATH` (or `WHALE_CMD` set to it).
   On a desktop or raven, `whale fetch-fastq --launch` does steps 6 and 7 in
   one go: it stages and marks the FASTQ as in step 6, then starts
   `nextflow run` for the claimed profile in the foreground (run it under
   `tmux` or `screen`), and exits with the run's status. It checks that it
   is in the study's worktree, that `nextflow` is on `PATH`, that the config
   sets the YAML's `chemistry` and that every file will be staged (no
   `--run`) before it downloads anything. On klone, submit the run as an `sbatch` head job
   instead ([compute profiles](compute_profiles.md)).

   **Check the library type after the first run** (or on one library first,
   for a deep study): the library kit is rarely in the deposit, and two kinds
   need a setting in the study config.
   - *Non-directional* (Zymo Pico Methyl-Seq and similar). FastQC's per-base
     content shows read 1 with as much C as G (~9% each); a directional
     library has under 1% C in read 1. Directional mode maps only half of
     such a library's pairs (ROBERTS2022: 30%). Add `--non_directional` to
     `bismark_align_args`.
   - *A synthetic tail on read 2* (adaptase-based kits). The M-bias table
     (`<sample>.dedup.M-bias.txt`) shows the first ~15 positions of read 2
     highly methylated in CpG **and** CHH context, and non-CpG methylation
     of 2-4% that looks like poor conversion (WANG2023). Clip it with
     `trim_args` (`--clip_r2 18 --three_prime_clip_r1 18`); do not lower
     `qc_min_conversion` for it.

   Every paired library also reads low at the first two bases of read 2
   (end repair); the workflow leaves them out of the calls
   (`params.ignore_r2 = 2`).
8. **Nothing to record by hand.** When the run succeeds it calls `whale
   record-run`, which commits `qc_summary.tsv`, `run_manifest.json`,
   `methylome_manifest.tsv`, `study_summary.yaml` and the gene x sample
   matrices (`genes/`, [gene-level methylation](gene_methylation.md)) under
   `data/studies/<STUDY_ID>/`, uploads the per-sample methylomes to the data
   host, moves `status.yaml` through `aligned`, `called`, `qc_passed` and
   `summarized`, pushes, and marks the PR ready (where the `gh` CLI is
   installed; raven has none, so mark a raven study's PR ready from another machine). The dashboard shows the
   study's finding from `study_summary.yaml` within three hours. If the
   upload failed (the machine has no key for gannet), the study stops at
   `qc_passed`; run `whale record-run --results <WHALE_METHYLOMES>/<STUDY_ID>`
   from a machine that can reach the host, or again once the key is in place.
   Then work through the PR checklist.
