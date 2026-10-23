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

   until it is clean, then `--strict`. Then `whale set-status --stage registered`.
5. **Sample sheet and FASTQ manifest** go in `data/studies/<STUDY_ID>/`
   (`samplesheet.tsv`: `sample_id`, `run_accession`, `fastq_1`, `fastq_2`).
   They are generated from ENA metadata, never typed. *(The generator,
   `whale fetch-samplesheet`, adapted from AREE, is Phase 1 work.)*
6. **Stage FASTQ** under `$WHALE_RAW/<STUDY_ID>/`, check MD5s against the
   manifest, `whale set-status --stage fastq_staged`.
7. **Set `chemistry`** in `config/<STUDY_ID>.config` to match the YAML and
   run the workflow with the profile you claimed. Record each stage as it
   finishes.
8. **Commit the small outputs** (`run_manifest.json`, QC tables) under
   `data/studies/<STUDY_ID>/`; copy per-sample methylomes to the data host.
   Reach `summarized`, work through the PR checklist, and mark the PR ready.
