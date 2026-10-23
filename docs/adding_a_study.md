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
6. **Stage FASTQ** under `$WHALE_RAW/<STUDY_ID>/`, check MD5s against the
   manifest, `whale set-status --stage fastq_staged`.
7. **Set `chemistry`** in `config/<STUDY_ID>.config` to match the YAML and
   run the workflow with the profile you claimed. Record each stage as it
   finishes.
8. **Commit the small outputs** (`run_manifest.json`, QC tables) under
   `data/studies/<STUDY_ID>/`; copy per-sample methylomes to the data host.
   Reach `summarized`, work through the PR checklist, and mark the PR ready.
