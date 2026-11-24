# Autopilot: running White Whale with limited human input

*Draft, 2026-11-18. Status: proposal, nothing built yet.*

## Goal

Move the backlog (Phase 4) and the analyses that follow it (Phases 5–6) from
"a person at a terminal drives each study" to "a scheduled controller drives
each study, calling the Claude API only where judgment is needed, and asks a
person only at a few named gates". A person should spend minutes a week, not
hours: read a digest, approve what is flagged, merge.

## Where we are

Much of the pipeline is already hands-off once a run starts:

| Step | Today | Already automatic? |
|---|---|---|
| Pick study, machine | Person reads `whale backlog`, runs `start-study` | No |
| Fill study YAML | Person reads ENA + paper | No (judgment) |
| Choose condition attributes | Person reads `--show-attributes` | No (judgment) |
| Set `chemistry` in config | Person edits (start-study writes `null`) | No |
| Mirror / stage FASTQ | `mirror-fastq` cron on raven; `fetch-fastq` by hand | Partly |
| Launch Nextflow | `fetch-fastq --launch` in tmux, or `sbatch` on klone | No |
| Recover failed runs | Person diagnoses (swap_mates, resume cache misses, preempted partitions) | No |
| Record results, QC, upload, PR ready | `whale record-run` | **Yes** |
| Review PR | Person reads `study_summary.yaml` / `qc_summary.tsv` | No |
| Merge | Person | No |
| Dashboard, atlas | `pages.yml` on push + every 3 h | **Yes** |

State today: 6 of 21 backlog studies `summarized`; 2 sit at `fastq_staged`
(DECICOMP2025 on klone, ROBERTSLAB2013 on Hete); 13 are unclaimed. The work
that remains is mostly the "No" rows above, repeated 13+ times.

## Design

### Principle: deterministic first, model second

Everything that can be a rule stays in the `whale` CLI and is tested. The
Claude API is called only for steps that currently need a person to *read*
something (deposit metadata, a paper, a Nextflow error log, a QC summary).
Each model call returns structured output that a deterministic validator
checks (`whale validate-study --strict`, the PR readiness check, the scope
check) before anything is committed. A model answer that fails validation is
retried once, then becomes a human gate.

### Components

1. **`whale autopilot tick`** — a new CLI command, idempotent, run on a
   schedule. For every study it reads the current stage (from `whale claims`
   and `status.yaml`) and attempts **one** transition. No long-lived
   process; a tick that crashes is simply retried by the next one.
2. **State machine** (per study), using the existing stages:

   ```
   free → claimed → registered → fastq_staged → running → summarized → reviewed → merged
                         ↑ model         ↑ deterministic   ↑ deterministic   ↑ model + rules
   ```

   plus `blocked` (needs a person) and `excluded` (with reason).
3. **Agent steps** (Claude API, tool use, structured JSON out):

   | Agent | Input | Output | Validated by |
   |---|---|---|---|
   | `register` | `fetch-samplesheet --show-attributes`, ENA/BioSample records, paper abstract/methods (PubMed tools) | Study YAML, condition attributes, label regex if needed, chemistry, quality flags | `validate-study --strict`; samplesheet groups have ≥2 replicates or `low_replication` flag |
   | `diagnose` | Failed run's `.nextflow.log`, `.command.err` of failed tasks, `status.yaml`, known-failure list | Classification (known remedy / transient / unknown) + proposed action | Only actions from an allow-list are executed |
   | `review` | `study_summary.yaml`, `qc_summary.tsv`, `snp_mask`, assembly, PR diff | Pass / concerns, with numbers quoted | Hard thresholds checked in code first (pass rate, mapping %, conversion %) |
   | `curate` | New SRA/ENA sweep results | New backlog rows | Backlog schema; PR to `main`, human-merged |

   Model choice: `claude-sonnet-5-5` for `register`, `review` and `curate`;
   `claude-opus-5-5` for `diagnose` of unknown failures. Use prompt caching
   for the stable system prompt (vocabularies, schema, readiness rules) so
   repeated ticks stay cheap.
4. **Known-failure playbook** (deterministic, no model): codify what we have
   already learned so it never needs a person again —
   - swap_mates + `-resume` reusing stale trim/align: move trim work dirs aside, relaunch;
   - Hete resume cache misses: export `WHALE_RAW/REFERENCE/METHYLOMES/WORKDIR` before relaunch;
   - klone preemption / wrong account: pick partition by the per-partition account table; never `scontrol` move Nextflow tasks;
   - gannet upload failed (no key): rerun `record-run --results` from raven;
   - transient network / ENA: retry with backoff, fall back to mirror.

   `diagnose` is called only when no playbook entry matches; a remedy it
   finds that works twice is added to the playbook by PR.
5. **GitHub through the REST API**, not `gh` (raven has no `gh`). A
   fine-grained token scoped to this repo: contents, pull requests, issues.
6. **Notification**: one daily digest (GitHub issue comment on a pinned
   "Autopilot log" issue, optionally mirrored to Slack): what advanced, what
   is running, what is blocked and the one question each blocked study needs
   answered.

### Where it runs

| Host | Role | Schedule |
|---|---|---|
| **raven** | Controller: `autopilot tick`, model calls, GitHub API, launches raven runs, `ssh klone` to submit `sbatch` head jobs, already runs the mirror cron | cron, hourly |
| **klone** | Executes deep studies; head job runs `record-run` on success (already) | via `sbatch` from raven |
| **Hete / desktop** | Optional worker for small studies; its own hourly tick limited to `--profile desktop` | cron / launchd |
| **GitHub Actions** | CI, scope check, dashboard/atlas (unchanged); add a `review` job on study PRs that posts the review agent's verdict | on PR |

The controller never runs two transitions on the same study at once (lock
file per study in `$WHALE_WORKDIR/autopilot/`), and respects per-machine
limits: max concurrent runs, free-disk floor (`fetch-fastq` already refuses
when short), and a klone core-hour budget.

### Human gates (the only things that wait for a person)

1. **Merging a study PR** — at first. After ~5 studies where the review
   agent and a person agreed, allow auto-merge for PRs that pass every hard
   threshold *and* the review agent finds no concerns; anything else stays
   manual.
2. **Any change to code, workflows, modules, schemas or the backlog** — always
   a PR for a person. The autopilot may open these, never merge them.
3. **Open decisions** in PLAN.md §10 (e.g. PRJEB60400 granularity) and any
   study the `register` agent marks ambiguous (design only in a paper table,
   pooled samples, non-*C. gigas* runs).
4. **Exclusions** — the agent proposes `excluded` with a reason; a person confirms.
5. **Budget overrun** — monthly API spend or klone hours above a set cap pause the controller.

### Guardrails

- The agent's write actions are an allow-list of `whale` subcommands plus
  "open PR / comment / label". No shell from model output, no force-push, no
  deleting branches or data, no writes outside the study's own files (the
  existing scope check enforces this on study branches).
- Every model call and its validated output is logged under
  `data/studies/<STUDY_ID>/autopilot_log.jsonl` (committed on the study
  branch) so the provenance of each YAML field is auditable.
- Model-filled YAML fields carry `source: autopilot` so a reader can tell
  them from hand-curated ones; the methods page reports how many.
- Kill switch: a `AUTOPILOT_PAUSED` file in the repo root (or label on the
  log issue) makes every tick a no-op.

## Rollout

| Step | Deliverable | Done when |
|---|---|---|
| **A. Deterministic tick** | `whale autopilot tick --dry-run` printing the next transition per study; playbook for the known failures; launch + relaunch from raven (incl. `ssh klone sbatch`) | Dry-run output matches what a person would do for every study in `whale claims` |
| **B. Unattended runs** | Tick live for `fastq_staged → running → summarized` only; digest issue | DECICOMP2025 and ROBERTSLAB2013 reach `summarized` with no terminal session |
| **C. Register agent** | `register` step with validators; human approves each generated YAML (PR comment "approve") | 3 studies registered where the person made no edits |
| **D. Review agent** | PR review job + hard thresholds; still human merge | Agent and person agree on 5 consecutive PRs |
| **E. Full loop** | Claim → merge for small/medium WGBS studies; auto-merge enabled under gate 1 | Backlog drains with only digests and gate approvals |
| **F. Curate + Phase 6** | Monthly SRA/ENA sweep agent; paired RNA-seq workflow (built by people, then driven by the same tick) | New deposits appear as backlog PRs; regulation pages rebuild as studies merge |

Steps A–B need no model calls and remove most of the current hands-on time;
C–E add the API where reading is required.

## Cost sketch (to be measured in step C)

Per study: one `register` call (tens of thousands of input tokens with the
paper's methods; cached system prompt), zero to a few `diagnose` calls, one
`review` call. For ~15 remaining studies plus reruns this is small next to
compute; set an initial monthly cap and report actual spend in the digest.

## Open questions for the PI

1. Auto-merge after step D, or keep every study merge manual?
2. Digest channel: GitHub issue only, or also Slack/email?
3. Should Hete keep running studies, or should raven + klone carry the backlog alone?
4. Monthly API budget cap.
5. Does the untracked `code/download_methylation_data*` (NCBI E-utilities
   downloader) get retired in favor of `fetch-fastq` / `mirror-fastq`?
