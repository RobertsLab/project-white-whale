# Running studies on several machines

Lab desktops, raven and klone work through the backlog at the same time. Each
job is one study on its own branch, checked out to the machine that will run
it and merged back as a pull request.

```
study/<ACCESSION>/<STUDY_ID>                 a whole study
study/<ACCESSION>/<STUDY_ID>.part-<PART>     one part of a registered study
```

The model is AREE's ([parallel_study_runs.md](https://github.com/sr320/AREE/blob/main/docs/parallel_study_runs.md)),
with a status file on every branch so the ledger shows where each study is.

## Checking a study out

From an up-to-date checkout on the machine that will run it:

```bash
whale backlog
whale start-study --accession PRJNA682817 --study-id ROBERTS2022_PH_PLOIDY_WGBS --profile raven
```

`start-study`:

1. fetches `origin/main` and refuses if the study or accession is already
   registered there, or already used by a pushed (or local) `study/` branch;
2. creates the branch in its own worktree, `../white-whale-<STUDY_ID>` by
   default (`--dir` to change it), so one machine can hold several studies;
3. writes `registry/studies/<STUDY_ID>.yaml` from the template,
   `config/<STUDY_ID>.config` for the chosen profile, and
   `data/studies/<STUDY_ID>/status.yaml` at stage `claimed`;
4. links this checkout's `data/raw`, `data/methylomes` and reference-bundle
   symlinks into the worktree;
5. commits, **pushes the branch (this is the claim)**, and opens a draft PR
   naming the machine and profile.

If two machines claim the same study at the same moment, only one push can
create the branch. The other is rejected and `start-study` stops before
anything is downloaded, printing the commands that remove its local worktree.

`--no-push` builds the worktree without claiming anything. `--no-pr` claims
without opening a PR (also the case when the `gh` CLI is absent).

## Recording progress

Inside the study's worktree:

```bash
whale set-status --stage fastq_staged --note "105 GB, MD5 ok"
whale set-status --stage aligned --run-name nostalgic_curie
```

Each call appends to the status file's history, commits and pushes. Stages
(`registry/controlled_vocabularies/pipeline_stage.yaml`) only move forward:

`claimed → registered → fastq_staged → aligned → called → qc_passed → summarized → merged`

plus `excluded` (will not be processed; say why in `--note`). `--force` moves
back, e.g. after a failed run is restarted from alignment. `--profile klone`
records a study moving to another machine mid-way.

## Seeing what is checked out

```bash
whale claims     # study, part, stage, profile, machine, who, when
whale backlog    # the backlog with each entry's state
```

Both read the status files from the remote branches; no server is involved.

## Giving a claim up

```bash
whale release-study --branch study/PRJNA682817/ROBERTS2022_PH_PLOIDY_WGBS --reason "moving to klone"
```

closes the PR with the reason and deletes the remote branch, so another
machine can claim the study. The local worktree is left alone in case it holds
downloaded data.

## Big studies in parts

A study too large for one machine (PRJEB60400, 246 runs) is claimed whole
first, registered (YAML, sample sheet, design), and that PR merged. Parts are
then claimed against the registered study:

```bash
whale start-study --accession PRJEB60400 --study-id POMS2023_ADAPTATION_WGBS --part gill --profile klone
whale start-study --accession PRJEB60400 --study-id POMS2023_ADAPTATION_WGBS --part mantle --profile raven
```

A part owns only `data/studies/<STUDY_ID>/parts/<PART>/` and
`config/<STUDY_ID>.part-<PART>.config`, so parts merge independently. A part
cannot be claimed while the whole study is still claimed.

## Why study PRs don't conflict

| Claim | Files it may change |
|---|---|
| whole study | `registry/studies/<STUDY_ID>.yaml`, `config/<STUDY_ID>.config`, `data/studies/<STUDY_ID>/**` except `parts/`, `tests/test_<study_id>.py` |
| part | `data/studies/<STUDY_ID>/parts/<PART>/**`, `config/<STUDY_ID>.part-<PART>.config` |

CI's `study-scope` job runs `whale check-study-scope` on every PR from a
`study/` branch and fails if the PR touches anything else, if the YAML's
`study_id` differs from the branch, or if the branch's accession is not in the
YAML. A code, vocabulary or docs change a study needs goes in its own PR,
merged first; then rebase the study branch onto `main`.

## After merging

On `main`, `whale build-registry`; then remove the worktree on the machine
that ran it:

```bash
git worktree remove ../white-whale-<STUDY_ID>
```
