# Mirroring the backlog's FASTQ to gannet

Downloading is the slowest step of a study and needs no compute, so it runs
ahead of the claims: `whale mirror-fastq`, from raven, a little every night,
until every backlog study's FASTQ is on gannet. A claimed study then starts
from a local copy instead of a multi-day download.

## Where the reads come from

Measured from raven on 2026-10-23, ENA served FASTQ at 0.04–1 MB/s per
connection. NCBI's AWS open-data bucket (`sra-pub-run-odp`) served the same
run at ~25 MB/s. Every backlog run sampled, the ENA-submitted `ERR` runs
included, is in that bucket as a full-quality `.sra` (not SRA Lite) with an
MD5 published by NCBI's locator service. So by default (`--source ncbi`) each
run is:

1. located with NCBI's locator (`locate.ncbi.nlm.nih.gov/sdl/2/retrieve`),
   which gives the S3 link, size and MD5 of the `.sra`;
2. downloaded into staging (resuming a partial file) and checked against that
   MD5;
3. converted on raven with `fasterq-dump --split-3` and gzipped (pigz if
   installed, else gzip with one process per mate) under **ENA's file names**
   (`<RUN>_1.fastq.gz`, `<RUN>_2.fastq.gz`), so sample sheets need no change;
4. checked: both mates have the same number of reads, and reads and bases
   equal ENA's `read_count` and `base_count` for the run. An unpaired leftover
   file counts towards those totals and is then discarded, as ENA's third
   file is.

A run falls back to ENA's FASTQ (downloaded over HTTPS and checked against
ENA's MD5) when NCBI has no full-quality copy, a check fails, or only some
of its files are missing. `--source ena` uses ENA for everything.

The file list itself always comes from ENA's read-run report. For a
registered study, the runs are those in its committed `fastq_manifest.tsv`.
Otherwise every run with a methylation library strategy is taken
(`--library-strategy` to change), taxid 29159 only. Other runs in the
project (often the RNA-seq) are reported, not mirrored.

## On gannet

```
raw/<ACCESSION>/<file>                  every FASTQ, named as ENA names it
raw/<ACCESSION>/mirror_manifest.tsv     the ledger, one row per verified file
```

served at `https://gannet.fish.washington.edu/v1_web/project-white-whale/raw/`.
These are public data; nothing private is mirrored. Ledger columns:

| Column | |
|---|---|
| `run_accession`, `mate`, `file` | as in the study's `fastq_manifest.tsv` |
| `url` | where it came from: the S3 `.sra` or ENA's FASTQ |
| `md5`, `bytes` | **of the file on gannet**, recomputed there after upload |
| `source` | `ncbi_sra`, `fastq_ftp` or `submitted_ftp` |
| `ena_md5` | ENA's MD5 for the file at the time it was mirrored |
| `sra_md5`, `read_count` | for converted runs: the `.sra` MD5 and reads per mate |
| `mirrored_on` | date |

A converted FASTQ holds the same reads and qualities as ENA's, but its read
names and gzip bytes differ, so **its MD5 is not ENA's**. Check copies taken
from gannet against the ledger's `md5`, not the study's `fastq_manifest.tsv`.

A file is skipped when the ledger has it against ENA's current MD5 and it is
on gannet at its recorded size. If ENA reprocesses a run (new MD5), the next
run fetches it again. A second copy started while one is running exits at
once, so overlapping cron runs are harmless. The exit code is 1 when any run
failed.

## How far along it is

```bash
whale mirror-fastq --plan
```

lists, per accession, the files wanted, how many are mirrored and verified,
and the GB left; it downloads nothing. On 2026-10-23 the whole backlog came
to **~4.3 TB** of gzipped FASTQ by ENA's byte counts (`fastq_gb` in
`backlog.tsv` underestimates it).

## Running it on raven

Once, in a checkout on raven (`pip install -e .` into a venv):

```bash
export WHALE_SRATOOLKIT=/home/shared/sratoolkit.3.0.2-ubuntu64
export WHALE_PIGZ=/home/shared/pigz-2.8/bin/pigz
export WHALE_MIRROR_STAGING=/home/shared/8TB_HDD_03/$USER/whale/mirror-staging
whale mirror-fastq --plan
```

Disk: staging holds the `.sra` and the gzipped FASTQ of the runs in flight.
`--tmp` holds fasterq-dump's scratch and the uncompressed FASTQ, roughly 3–4×
the `.sra`, so 100–200 GB for one of the deepest runs. Put `--tmp` on the
NVMe RAID0 (`/home/shared/16TB_HDD_01/$USER/whale/mirror-tmp`) and check
`df -h` there before raising `--jobs`, since the workflows use it too.

CPU: compression is the slow step after conversion. pigz uses `--threads`
cores per run (split between the mates); without it, gzip uses one core per
mate. On raven pigz 2.8 is unpacked from Ubuntu's package into
`/home/shared/pigz-2.8` (no root needed: `apt-get download pigz`,
`dpkg-deb -x`). It compressed a 1.65 GB FASTQ in 22 s with 4 threads against
123 s for gzip.

Without `WHALE_PIGZ`, pigz is used only if it is on `PATH`.

Nightly from cron (`crontab -e`), capped so one night cannot fill staging
or take over the machine:

```cron
# 20:00 every night: up to 500 GB, 3 runs at a time
0 20 * * * cd $HOME/project-white-whale && WHALE_SRATOOLKIT=/home/shared/sratoolkit.3.0.2-ubuntu64 WHALE_PIGZ=/home/shared/pigz-2.8/bin/pigz WHALE_MIRROR_STAGING=/home/shared/8TB_HDD_03/sr320/whale/mirror-staging .venv/bin/whale mirror-fastq --jobs 3 --threads 4 --tmp /home/shared/16TB_HDD_01/sr320/whale/mirror-tmp --max-gb 500 >> $HOME/whale-mirror.log 2>&1
```

Cron needs raven's SSH key for gannet to work without a passphrase, which it
does (BatchMode is used). Or run it by hand in `tmux`: the same command
without `--max-gb`.

| Option | Meaning |
|---|---|
| `--accession PRJ…` | Only these (repeatable), e.g. the next studies to be claimed |
| `--source ncbi\|ena` | `ncbi` (default) converts NCBI's `.sra`, with ENA as fallback per run; `ena` uses ENA's FASTQ only |
| `--jobs N` | Runs at once (default 4) |
| `--threads N` | fasterq-dump / pigz threads per run (default 4) |
| `--tmp DIR` | fasterq-dump scratch (default `<staging>/tmp`) |
| `--limit-rate MB/s` | Total download cap across all jobs |
| `--max-gb G` | Stop starting new runs once G GB (ENA's sizes) are queued in this run |
| `--no-verify-remote` | Skip the MD5 on gannet (rsync already checks what it transfers) |
| `--dest DIR` | Mirror into a local directory instead (a mounted share, or testing) |

## Using the mirror

`whale fetch-fastq` (step 6 of [adding a study](adding_a_study.md)) uses the
mirror by default. For the study's accession it reads
`raw/<ACCESSION>/mirror_manifest.tsv` from `url_base` in
`data/data_host.yaml`, and for each file in the study's `fastq_manifest.tsv`:

- **mirrored against the manifest's current ENA MD5:** copied from gannet over
  HTTPS (resumable) and checked against the ledger's `md5`;
- **not mirrored, or mirrored against an older ENA MD5** (ENA reprocessed the
  run since): downloaded from ENA and checked against ENA's MD5;
- **mirror copy fails** (download or MD5): deleted and fetched from ENA instead.

A partial file is only resumed from the source it was started from. Files
already verified from either source are kept. The study's
`fastq_staged.tsv` records, per file, the MD5 actually checked, the
`source` (`mirror` or `ena`) and ENA's MD5. `--source mirror` refuses to
start unless every file is usable from the mirror; `--source ena` ignores it;
`--mirror-url` points at another copy.
