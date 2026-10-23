# Compute profiles: desktop, raven, klone

A study is claimed for one profile (`whale start-study --profile …`), and the
workflow runs with the matching Nextflow profile from `config/base.config`.
No config file names a path on one machine; each machine sets environment
variables instead.

| Variable | Meaning | Default |
|---|---|---|
| `WHALE_RAW` | directory holding `<STUDY_ID>/` FASTQ folders | `data/raw` |
| `WHALE_REFERENCE` | directory holding `<ASSEMBLY>/genomic.fna` etc. | `data/reference` |
| `WHALE_METHYLOMES` | where per-sample outputs are published | `data/methylomes` |
| `WHALE_WORKDIR` | Nextflow work directory, on fast local disk | `work` (klone: scrubbed) |
| `WHALE_MAX_CPUS`, `WHALE_MAX_MEMORY` | local executor limits (desktop, raven) | 10 / 28 GB; 32 / 120 GB |
| `WHALE_SLURM_ACCOUNT`, `WHALE_SLURM_QUEUE` | klone | `srlab`, `ckpt-all` |
| `WHALE_APPTAINER_CACHEDIR` | klone | `/gscratch/scrubbed/$USER/whale/apptainer-cache` |

Alternatively symlink `data/raw`, `data/methylomes` and
`data/reference/GCF_963853765.1` to the machine's drives; `start-study`
carries those links into each new worktree. Avoid exFAT and paths with spaces.

## desktop

Lab Macs and workstations: the local executor, one Bismark alignment at a
time. Good for studies up to ~150 GB of FASTQ, and for QC, summaries and the
analysis notebooks.

```bash
export WHALE_WORKDIR=/Volumes/<fast-ssd>/whale/work
nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile desktop,docker
```

## raven

The Roberts Lab server: the local executor with more cores and memory, two
alignments at once. Good for ~150-500 GB studies. Its core count, memory and
scratch path are not yet recorded (PLAN.md §10); until they are, set
`WHALE_MAX_CPUS`, `WHALE_MAX_MEMORY` and `WHALE_WORKDIR` in the shell that
starts the run, and use `docker` or `apptainer` as the engine.

```bash
nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile raven,docker
```

Run long jobs inside `tmux` or `screen` so a dropped SSH session does not stop
them.

## klone

UW Hyak's checkpoint partition, for the deep studies. Tasks are preemptible
and capped below Hyak's ckpt limit (4 h 45 min each) with six retries;
Apptainer is enabled by the profile. Lessons carried from AREE
([running_on_klone.md](https://github.com/sr320/AREE/blob/main/docs/running_on_klone.md)):

- Run the Nextflow **head job** on a non-preemptible partition (e.g. `cpu-g2`
  under the group account), 2 CPUs, 8 GB, a few days. Never on a login node.
- **Do not export `TMPDIR`** in the head job; Slurm forwards it to every task
  and into containers. Use `NXF_TEMP` for Nextflow's own temporary files.
- `/gscratch/scrubbed` deletes files untouched for 21 days and is often >97%
  full. Check `hyakstorage /gscratch/scrubbed` first; copy what you keep.
- Read-sized intermediates go to the node's `/scr`.
- ENA downloads through the cluster proxy ran at ~12 MB/s aggregate in
  October 2026: plan ~15 h per 600 GB, as a job.
- Select a subset of libraries with a directory of links, not a brace pattern
  in a glob.
- Put Klone-only settings (e.g. chunk size, once chunked alignment lands) in a
  run-only config passed last with `-c`, so a study config's workstation time
  limits do not override the profile.

```bash
#SBATCH --account=srlab --partition=cpu-g2 --cpus-per-task=2 --mem=8G --time=4-00:00:00
export NXF_HOME=/gscratch/scrubbed/$USER/whale/nxf-home
export NXF_TEMP=/gscratch/scrubbed/$USER/whale/tmp
export WHALE_RAW=/gscratch/scrubbed/$USER/whale/raw
export WHALE_REFERENCE=/gscratch/scrubbed/$USER/whale/reference
export WHALE_METHYLOMES=/gscratch/scrubbed/$USER/whale/methylomes
nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile klone
```

## Smoke test, any machine

```bash
nextflow run workflows/methylation -profile test -stub
```

runs the whole DAG on two tiny fixture samples without executing any tool.
