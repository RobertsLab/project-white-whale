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
- [ ] The run recorded itself (`whale record-run`): `run_manifest.json`, `qc_summary.tsv`, `methylome_manifest.tsv` and `study_summary.yaml` committed under `data/studies/<STUDY_ID>/`; no FASTQ, BAMs, `work/` or `results/`
- [ ] `status.yaml` is at `summarized` (methylomes on the data host), or `excluded` with the reason in `notes`
- [ ] Libraries failing the QC gate (`qc_summary.tsv`, `qc_reasons`) checked against the deposit
- [ ] Quality flags and limitations record pooling, confounding, missing spike-ins and excluded libraries
- [ ] No edits to shared code, vocabularies or docs; those went in a separate PR

## Curation decisions

<!-- Excluded libraries, sample-to-group mapping, M-bias trimming, anything a
     reviewer should check against the deposit and the paper. -->
