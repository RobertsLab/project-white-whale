# Compute profiles: desktop, raven, klone

A study is claimed for one profile (`whale start-study --profile …`), and the
workflow runs with the matching Nextflow profile from `config/base.config`.
No config file names a path on one machine; each machine sets environment
variables instead.

| Variable | Meaning | Default |
|---|---|---|
| `WHALE_RAW` | directory holding `<STUDY_ID>/` FASTQ folders | `data/raw` |
| `WHALE_REFERENCE` | directory holding `<ASSEMBLY>/genomic.fna` etc. | `data/reference` |
| `WHALE_BISMARK_INDEX` | directory holding the cached `bismark_genome/` index | the reference bundle |
| `WHALE_METHYLOMES` | where per-sample outputs are published | `data/methylomes` |
| `WHALE_WORKDIR` | Nextflow work directory, on fast local disk | `work` (klone: scrubbed) |
| `WHALE_MAX_CPUS`, `WHALE_MAX_MEMORY` | local executor limits (desktop, raven) | 10 / 28 GB; 32 / 120 GB |
| `WHALE_SLURM_ACCOUNT`, `WHALE_SLURM_QUEUE` | klone | `srlab`, `ckpt-all` |
| `WHALE_APPTAINER_CACHEDIR` | klone | `/gscratch/scrubbed/$USER/whale/apptainer-cache` |

Alternatively symlink `data/raw`, `data/methylomes` and
`data/reference/GCF_963853765.1` to the machine's drives; `start-study`
carries those links into each new worktree. Avoid exFAT and paths with spaces.

## Every machine, once

```bash
whale fetch-reference     # into $WHALE_REFERENCE/GCF_963853765.1 (docs/reference_bundle.md)
```

## Where FASTQ comes from, and how fast

Speed depends on the path between the machine and the source, not on the
tools. Measured 2026-10-23 (MB/s):

| Source → machine | Speed | Measured by |
|---|---|---|
| ENA HTTPS → lab Mac | ~0.09 per connection | `curl`, 30 s |
| ENA HTTPS → raven | ~0.04–1 per connection | `curl`, 30 s; [FASTQ mirror](fastq_mirror.md) |
| NCBI's AWS bucket (`sra-pub-run-odp`, `.sra`) → raven | ~25 | [FASTQ mirror](fastq_mirror.md) |
| gannet mirror → lab Mac | ~2–3 total, 4 files at once | `whale fetch-fastq`, 530 MB in 4 min |
| ENA → klone | ~12 aggregate through the proxy | AREE, October 2026 |

NCBI's bucket is the fast source, but its runs arrive as `.sra` and must be
converted (`fasterq-dump`). So the lab's FASTQ reaches studies in two steps:

1. `whale mirror-fastq` on raven pulls each backlog study from NCBI's bucket
   (ENA as fallback), converts it, verifies it and puts it on gannet, ahead
   of any claim ([FASTQ mirror](fastq_mirror.md)).
2. `whale fetch-fastq` on the machine that runs the study copies the files
   from gannet, and takes from ENA only what the mirror does not have yet.

A study that is not mirrored yet is still faster to mirror on raven first
than to fetch from ENA on a desktop: Wang 2014's 36 GB is days from ENA on
the Mac.

## desktop

Lab Macs and workstations: the local executor, one Bismark alignment at a
time. Good for studies up to ~150 GB of FASTQ, and for QC, summaries and the
analysis notebooks.

```bash
export WHALE_WORKDIR=/Volumes/<fast-ssd>/whale/work
nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile desktop,docker
```

## raven

The Roberts Lab server (`raven.fish.washington.edu`), inspected 2026-10-23:

| | |
|---|---|
| CPU | Intel Xeon Gold 5220R, 24 cores / 48 threads |
| Memory | 247 GB (+ 99 GB swap) |
| OS | Ubuntu 24.04 |
| Fast disk | `/home/shared/16TB_HDD_01`: 14 TB NVMe RAID0 (8 × 1.7 TB) despite its name; ~1.5 TB free. Fast, **no redundancy** |
| Bulk disk | `/home/shared/8TB_HDD_03` (3.1 TB free), `/home/shared/8TB_HDD_01` (2.1 TB free), HDD |
| Containers | **Singularity CE 3.10** (`/usr/bin/singularity`, setuid): runs local images and pulls `docker://` images. Docker needs membership of the `docker` group; the Apptainer 1.2.4 in `/home/shared/apptainer` fails on unprivileged user namespaces |
| Java / Nextflow | OpenJDK 21; **Nextflow not installed** |

The machine is shared (load ~8 of 48 when inspected), so the profile leaves
headroom: 36 CPUs and 180 GB for the run, 12 CPUs per alignment, three
alignments at once. It enables Singularity itself and defaults the work
directory and image cache to the NVMe RAID:

| Setting | Default | Override |
|---|---|---|
| work dir | `/home/shared/16TB_HDD_01/$USER/whale/work` | `WHALE_WORKDIR` |
| image cache | `/home/shared/16TB_HDD_01/$USER/whale/singularity-cache` | `WHALE_SINGULARITY_CACHEDIR` |
| executor limits | 36 CPUs, 180 GB | `WHALE_MAX_CPUS`, `WHALE_MAX_MEMORY` |

Keep FASTQ and published methylomes on a bulk disk, not the RAID0:

```bash
# one-time, if Nextflow is not on PATH (Java 21 is already installed)
curl -s https://get.nextflow.io | bash && mkdir -p ~/bin && mv nextflow ~/bin/

export WHALE_RAW=/home/shared/8TB_HDD_03/$USER/whale/raw
export WHALE_METHYLOMES=/home/shared/8TB_HDD_03/$USER/whale/methylomes
export WHALE_REFERENCE=/home/shared/8TB_HDD_03/$USER/whale/reference
nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile raven
```

Good for studies of ~150-500 GB; the work directory needs about twice the
FASTQ size, so check `df -h /home/shared/16TB_HDD_01` before a large run. Run
long jobs inside `tmux` or `screen` so a dropped SSH session does not stop
them. Bismark 0.24.0 and BS-SNPer are also installed under `/home/shared`,
but the workflow uses its pinned containers, not those.

## klone

UW Hyak's checkpoint partition, for the deep studies. Tasks are preemptible
and capped at the ckpt window (4.5 h each) with six retries;
Apptainer is enabled by the profile. Lessons carried from AREE
([running_on_klone.md](https://github.com/sr320/AREE/blob/main/docs/running_on_klone.md)):

- Run the Nextflow **head job** on a non-preemptible partition (`coenv`
  `compute` first; see [Choosing a partition](#choosing-a-partition)),
  2 CPUs, 8 GB, a few days. Never on a login node.
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
#SBATCH --account=coenv --partition=compute --cpus-per-task=2 --mem=8G --time=4-00:00:00
export NXF_HOME=/gscratch/scrubbed/$USER/whale/nxf-home
export NXF_TEMP=/gscratch/scrubbed/$USER/whale/tmp
export WHALE_RAW=/gscratch/scrubbed/$USER/whale/raw
export WHALE_REFERENCE=/gscratch/scrubbed/$USER/whale/reference
export WHALE_METHYLOMES=/gscratch/scrubbed/$USER/whale/methylomes
nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile klone
```

### Choosing a partition

Use the checkpoint partition, `ckpt-all`, for as much as possible. It runs
on idle nodes across Hyak, so it is the cheapest place to queue many tasks,
but it can preempt a job at any time and has a **4.5 h window**: no task
there may ask for more. The profile defaults every task to `ckpt-all` at
4 h 30 min, retried six times.

Anything that needs longer than 4.5 h goes to a non-preemptible partition,
in this order:

| Order | Account | Partition | Lab's share (max) |
|---|---|---|---|
| 1 | `coenv` | `compute` | 40 CPUs, 188 GB |
| 2 | `srlab` | `cpu-g2-mem2x` | 32 CPUs, 503 GB |
| 3 | `coenv` | `cpu-g2` | 32 CPUs, 200 GB |
| 4 | `coenv` | `cpu-g2`, a second job of the same size | 32 CPUs, 200 GB |

Take the first one with room (`hyakalloc` shows what each account has in use
right now). Skip straight to `srlab cpu-g2-mem2x` when one task needs more
memory than 200 GB. The maxima are for everything running under that
account and partition at once, the head job included, so a study's overrides
must fit inside them alongside other lab work.

| Work | Where |
|---|---|
| FastQC, trimming, chunked alignment, extraction: every task that finishes in 4.5 h | `ckpt-all` (profile default) |
| Nextflow head job (2 CPUs, 8 GB, days; a preempted head job stops the run) | first partition in the order above with room |
| A task that cannot finish in 4.5 h | first partition with room, as a per-process override |
| `fetch-fastq` downloads, interactive tests | `ckpt-all` via `sbatch` / `salloc`, never a login node |

Before moving a task off `ckpt-all`, try splitting it (as alignment is
chunked) so it fits the window. Signs that a task needs to move:
`sacct -u $USER -S today --format=JobName,Partition,State,Elapsed` shows it
failing at the time limit, or being preempted again and again.

The default account and partition come from `WHALE_SLURM_ACCOUNT` and
`WHALE_SLURM_QUEUE`. Per-process overrides go in a run-only config passed
after the study config, never in the study config itself. Set the account
with the partition, since the profile's `--account=srlab` applies to every
task otherwise, and cap `maxForks` so the override fits the lab's share:

```groovy
// config/<STUDY_ID>.klone.config
// nextflow run workflows/methylation -c config/<STUDY_ID>.config -c config/<STUDY_ID>.klone.config -profile klone
process {
    withName: 'BISMARK_ALIGN' {
        queue          = 'compute'
        clusterOptions = '--account=coenv'
        time           = '24h'
        maxForks       = 3    // 3 x (12 CPUs, 40 GB) + head job fits 40 CPUs, 188 GB
    }
}
```

If `coenv compute` is full, change `queue` and `clusterOptions` to the next
row of the order (`cpu-g2-mem2x` with `--account=srlab`, then `cpu-g2` with
`--account=coenv`) and lower `maxForks` to fit 32 CPUs.

## Smoke test, any machine

```bash
nextflow run workflows/methylation -profile test -stub
```

runs the whole DAG on two tiny fixture samples without executing any tool.

```bash
make real-run      # nextflow run workflows/methylation -profile test,docker, then the check
```

runs every tool for real, in the pinned containers, on the same fixtures: a
30 kb synthetic genome and 1,500 bisulfite read pairs (S1) and single reads
(S2) simulated with known methylation by `tests/fixtures/make_fixtures.py`.
`tests/check_real_run.py` then checks that nearly every read maps, that
CpG methylation comes out near the simulated 83%, and that CHG/CHH stay near
0%. It takes a few minutes and needs Docker; CI runs both on every PR.
