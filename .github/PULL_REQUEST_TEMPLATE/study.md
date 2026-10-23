<!-- Study PR: branch study/<ACCESSION>/<STUDY_ID>[.part-<PART>], opened by `whale start-study`.
     It may change only files this study owns; CI enforces that (check-study-scope).
     See docs/running_studies_on_several_machines.md. -->

## Study

- **study_id:**
- **Accession(s):**
- **Chemistry:**
- **Machine, profile and storage used:**

## Checklist

- [ ] `whale validate-study registry/studies/<STUDY_ID>.yaml --strict` passes (whole-study claims)
- [ ] Sample sheet and FASTQ manifest generated from ENA/SRA metadata, not hand-written
- [ ] FASTQ MD5s verified against the manifest before the run
- [ ] `run_manifest.json` committed under `data/studies/<STUDY_ID>/`, with parameters, tool versions and input checksums
- [ ] Per-sample QC table committed; no FASTQ, BAMs, `work/` or `results/`
- [ ] Per-sample methylomes copied to the data host and listed with checksums in `methylome_manifest.tsv`
- [ ] `status.yaml` is at `summarized` (or `excluded`, with the reason in `notes`)
- [ ] Quality flags and limitations record pooling, confounding, missing spike-ins and excluded libraries
- [ ] No edits to shared code, vocabularies or docs; those went in a separate PR

## Curation decisions

<!-- Excluded libraries, sample-to-group mapping, M-bias trimming, anything a
     reviewer should check against the deposit and the paper. -->
