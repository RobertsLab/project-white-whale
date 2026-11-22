# Autopilot

The autopilot is a controller, run on a schedule, that moves each study one
step per run (a *tick*). It calls the Claude API only where something has to
be read, and asks a person only at a few named gates. The design and the
open questions for the PI are in the plan (`docs/autopilot_plan.md`, PR #59).
This page describes what is built: **step A**, a deterministic dry run.

## Step A: `whale autopilot tick --dry-run`

```bash
whale autopilot tick --dry-run                 # one line per study, then the Markdown digest
whale autopilot tick --dry-run --markdown      # the digest only
whale autopilot tick --dry-run --digest digest.md --no-fetch
```

The command changes nothing. It launches nothing, pushes nothing, comments
on no PR and makes no model call. It refuses to run without `--dry-run`.
If a file named `AUTOPILOT_PAUSED` is in the repository root, every tick
prints its contents and does nothing else (the kill switch). The one
network call is `git fetch` of main and the `study/` branches, which only
updates this checkout's remote-tracking refs (`--no-fetch` skips it).

For every study the tick reads:

| Source | What |
|---|---|
| `registry/backlog.tsv` on main | every study, its suggested profile, size and notes |
| `registry/studies/` on main | what is registered |
| each pushed `study/` branch | `status.yaml`, study YAML, `config/<ID>.config`, `data/studies/<ID>/*.run.config`, `qc_summary.tsv`, `study_summary.yaml`, `run_manifest.json` |
| `git merge-base --is-ancestor` | whether the branch is merged into main |
| the study's worktree, if it is on this machine | `.nextflow.log`, `.nextflow/history`, `.nextflow.pid`, the failed task's `.exitcode` and `.command.err`, other `*.out` / `*.log` files |
| `config/autopilot.yaml` | per-machine limits |
| `squeue -u $USER` on each sbatch machine (klone) | its jobs: run locally on that machine, otherwise `ssh -o BatchMode=yes <ssh_host>` (`src/whale/autopilot/slurm.py`) |

A Nextflow log on another machine is not visible, but klone's queue is. A
job belongs to a study when its name or the last component of its working
directory contains the study id (the head job runs in the study checkout,
`white-whale-<STUDY_ID>`), or when a `whale-head-<PREFIX>` job's prefix names
exactly one study on klone. A staged study with such a job is `running`, takes
a run slot, and is never relaunched. If the queue cannot be read (no SSH
ControlMaster session, so Duo would have to prompt), a staged klone study
waits instead of advancing: the tick will not propose a launch that might
duplicate a head job it could not see. For any other machine the decision
says what to check before launching ("once no Nextflow is running it").

`fastq_staged` does not mean a study has never run: a klone run that failed
(DECICOMP2025 hit a scrubbed disk quota on 2026-11-21) leaves an empty queue
and a log the tick cannot read. So every launch passes `-resume` (with the
session id when this machine has the log): on a first launch it caches
nothing, after a failure it keeps the finished work. When the config sets
`swap_mates`, the decision also warns that `-resume` would reuse unswapped
trim/align tasks of a run that predates it.

### The state machine

```
free → claimed → registered → fastq_staged → running → summarized → reviewed → merged
                                                             ↘ blocked (one question for a person)
excluded (terminal)
```

These states are coarser than `pipeline_stage.yaml`. `record-run` writes
`aligned`, `called` and `qc_passed` in one go, so a study at any of them is
still `running` until it reaches `summarized`. Each tick picks **one**
transition per study (`src/whale/autopilot/decide.py`):

| State | Next transition | Who |
|---|---|---|
| free | `whale start-study` on the first machine for the backlog's profile that has room; at most one claim per machine per tick | autopilot |
| free, unnamed / "claim in parts" / needs a run list | blocked: name it, split it, list its runs | person (register agent in step C) |
| free, enrichment chemistry | waiting for the enrichment workflow on main | — |
| claimed | blocked: registration (study YAML, condition attributes) | person (step C) |
| registered | `whale fetch-fastq --mark-staged` (sbatch on klone), after the chemistry check and the disk floor | autopilot |
| fastq_staged | launch the run (tmux locally, `sbatch` head job on klone) with `-resume`, if the machine has a free run slot | autopilot |
| running, live | wait | — |
| running, finished | `whale record-run --results …` | autopilot |
| running, failed | the playbook's remedy, else blocked "needs diagnosis" with the evidence | autopilot / person (step C: model) |
| summarized | readiness check → `reviewed`; any failed check → blocked with one question | autopilot |
| reviewed | merge | person (gate 1) |
| merged, excluded | nothing; library checks may propose a re-run | — |

Each decision is one of: *advance* (the autopilot would do it), *running*,
*blocked* (needs a person; it always carries the one question to answer),
*waiting* (for a slot, disk space or a workflow), *done*.

### Machine limits

`config/autopilot.yaml` gives each machine its profile,
`max_concurrent_runs`, `disk_floor_gb` / `disk_path`, how it launches
(`local` or `sbatch`), whether it has a gannet SSH key and `gh`, the
`WHALE_*` variables to export before any launch or `-resume`, extra `-c`
configs (Hete's `hete-limits.config`), klone's per-partition accounts, and
klone's `ssh_host` for reading its queue. Accounts are written in full
(`compute-coenv`, `cpu-g2-mem2x-srlab`, `cpu-g2-coenv`): `sbatch` also takes
the short group name, but `scontrol` and `squeue` use the full one.
Step A uses them only to decide what *would* be launched where. Free disk
is measured only on the machine the tick runs on.

### Playbook of known failures

`src/whale/autopilot/playbook.py`. The first entry that matches gives the remedy:

| Entry | Recognized by | Remedy |
|---|---|---|
| `swap_mates_resume` | a `-resume` run cached `TRIM_GALORE` / `BISMARK_ALIGN` for a library in `swap_mates` (task order is not in the hash) | cancel before record-run, move those work dirs aside (paths from the cached task hashes), relaunch. Raised even when the run succeeded |
| `resume_env_cache_miss` | a `-resume` run that cached nothing, with a work dir inside the checkout or other than the machine's `WHALE_WORKDIR`, or earlier runs of the session | export `WHALE_RAW`, `WHALE_REFERENCE`, `WHALE_METHYLOMES`, `WHALE_WORKDIR`; relaunch with the session id |
| `gannet_upload` | stage `qc_passed` with `methylomes.on_data_host: false`, or an SSH key error in the log | rerun `whale record-run --results …` from a machine with a gannet key (raven) |
| `raven_no_gh` | `gh: command not found` on raven | do the GitHub step through the REST/GraphQL API |
| `klone_partition` | Slurm account/partition errors, repeated preemption, or `scontrol update … Partition=` | per-process override setting account and partition together in `klone.run.config`; never move tasks with `scontrol` |
| `chemistry_null` | the config's `chemistry` is null or disagrees with the study YAML (checked before staging and launch too) | set it to the YAML's term; if the YAML has none either, blocked |
| `oom_killed` | exit 137 / `Killed` | cap the process with `maxForks = 1` in the machine's limits config (or, if already capped, more memory or a bigger machine); relaunch |
| `transient_network` | connection resets, timeouts, 5xx, curl network codes | retry `fetch-fastq` with backoff; use the gannet mirror; mirror on raven first if needed |
| `head_stopped` | the log has no end and the `.nextflow.pid` process is gone (power cut, reboot) | relaunch with `-resume` and the machine's exports |

A failed run that matches no entry is reported as blocked, "needs
diagnosis". The report carries the failed process, its exit code and the
tail of its stderr.

### Library-type and readiness checks

`src/whale/autopilot/library.py` reads a finished run's `qc_summary.tsv` and
`run_manifest.json` and **proposes** config changes (never applies them):

- **Non-directional libraries** (Zymo Pico): every library maps below 40%
  in directional mode, tightly grouped. Proposal: confirm read 1 is ~9% C in
  FastQC, add `--non_directional` to `bismark_align_args`, re-run.
- **Synthetic read-2 tail**: at least 80% of libraries sit at 1.5-5%
  non-CpG methylation within a 2.5-point spread. Proposal: confirm in the
  read-2 M-bias plot, set `trim_args = '--clip_r2 18 --three_prime_clip_r1 18'`,
  keep the 98% conversion gate, re-run.
- When the config already has the fix but the recorded run did not use it,
  the proposal is to re-run with the current config.

The readiness check at `summarized` adds three more checks: the QC pass
rate (at least 75%), the SNP mask (without it the atlas leaves the study
out) and the genome assembly (the same as the merged studies).

## What later steps add

| Step | Adds | Model calls |
|---|---|---|
| **B. Unattended runs** | live transitions for `fastq_staged → running → summarized` only: launch, relaunch with allow-listed playbook remedies, `record-run`; a lock file per study; the digest posted to the "Autopilot log" issue through the REST API | none |
| **C. Register agent** | `claimed → registered`: study YAML, condition attributes and chemistry from the deposit and paper, checked by `validate-study --strict`; a person approves each. "Needs diagnosis" failures go to a `diagnose` call | yes |
| **D. Review agent** | `summarized → reviewed` with a written review on the PR; the hard thresholds above stay in code | yes |
| **E. Full loop** | claim to merge for small and medium WGBS studies; auto-merge only if the PI agrees (gate 1) | yes |
| **F. Curate + Phase 6** | a monthly SRA/ENA sweep that proposes backlog rows; the RNA-seq workflow driven by the same tick | yes |

Each step adds live transitions behind the same `decide` function. The dry
run stays available, so a person can always see what the next tick would do.

## Adding to it

- A new known failure: a function in `playbook.py` returning a `Match`,
  added to `ENTRIES`, with a fake-log test in
  `tests/test_autopilot_playbook.py`.
- A new transition or gate: `decide.py`, with a test in
  `tests/test_autopilot_decide.py`.
- A new machine: an entry in `config/autopilot.yaml`.
