# Code adapted from AREE

White Whale vendors code from [AREE](https://github.com/sr320/AREE) (MIT)
rather than depending on it (PLAN.md §10). Each entry names the AREE source
and the AREE commit it was taken from, so a fix in either project can be
carried across by hand.

| White Whale | AREE source | AREE commit | Changes |
|---|---|---|---|
| `src/whale/claims.py` | `src/intake/study_claims.py` | `44022f1` | BioProject-only accessions; per-branch `status.yaml`; compute profiles; part claims; `set_status`, `list_claims`, `release_study` |
| `src/whale/scope.py` | `src/validation/study_scope.py` | `44022f1` | Part ownership; test-file naming |
| `.github/PULL_REQUEST_TEMPLATE/study.md` | same path | `44022f1` | Methylome checklist |
| `config/base.config` (`klone` profile) | `config/base.config` | `44022f1` | `WHALE_*` variables; methylation processes only |
| `modules/methylation/*.nf` (containers, Bismark `--score_min L,0,-0.6`) | `modules/methylation/` | `44022f1` | Rewritten; container pins and Bismark settings kept |
| `data/reference/genome_assemblies.yaml` | same path | `44022f1` | Two assemblies; `default` flag |
| `docs/compute_profiles.md` (klone section) | `docs/running_on_klone.md` | `44022f1` | Summarized |

Still to vendor in Phase 1: `src/intake/ena_samplesheet.py` (as
`whale fetch-samplesheet`) and `src/intake/reference_bundle.py` (as
`whale fetch-reference`).

When vendoring more, add a row here with the AREE commit
(`git -C <AREE checkout> log -1 --format=%h -- <file>`).
