# Mirroring the backlog's FASTQ to gannet

Downloading from ENA is the slowest step of a study and needs no compute, so
it runs ahead of the claims: `whale mirror-fastq`, from raven, a little every
night, until every backlog study's FASTQ sits on gannet. Claiming a study then
starts from a local copy instead of a multi-day download.

## What it does

For each backlog accession (or each `--accession`), in backlog order:

1. reads the FASTQ list and MD5s from ENA's read-run report. For a registered
   study, the runs are those in its committed `fastq_manifest.tsv`; otherwise
   every run with a methylation library strategy (`--library-strategy` to
   change), taxid 29159 only. Other runs in the project (often the RNA-seq)
   are reported, not mirrored;
2. skips each file already in the ledger with ENA's current MD5 and present on
   gannet at its byte count;
3. downloads the rest over HTTPS into a staging directory, `--jobs` at a time,
   resuming any partial file left by an earlier run or a dropped connection;
4. checks each against ENA's MD5 (a mismatch is fetched once more, then
   reported and discarded), rsyncs it to gannet, recomputes the MD5 there, and
   adds it to the ledger.

On gannet (`data/data_host.yaml`):

```
raw/<ACCESSION>/<file>                  every FASTQ, named as ENA names it
raw/<ACCESSION>/mirror_manifest.tsv     run, mate, file, ENA URL, MD5, bytes, date mirrored
```

served at `https://gannet.fish.washington.edu/v1_web/project-white-whale/raw/`.
These are public ENA data; nothing private is mirrored.

If ENA reprocesses a run (new MD5), the next run notices and fetches it again.
A second copy started while one is running exits at once, so overlapping cron
runs are harmless. The exit code is 1 when any file failed.

## How far along it is

```bash
whale mirror-fastq --plan
```

lists, per accession, the files wanted, how many are mirrored and verified,
and the GB left. It downloads nothing. Measured 2026-10-23: **~4.3 TB** for the
whole backlog (ENA's gzipped byte counts; `fastq_gb` in `backlog.tsv` is an
underestimate).

## Running it on raven

Once, in a checkout on raven (`pip install -e .` into a venv):

```bash
export WHALE_MIRROR_STAGING=/home/shared/8TB_HDD_03/$USER/whale/mirror-staging
whale mirror-fastq --plan
```

Staging holds only the files in flight (`--jobs` of them, up to a few tens of
GB for the deep studies), so a bulk disk is enough. Keep it off the NVMe
RAID0, which the workflows use.

Nightly from cron (`crontab -e`), with a cap so one night cannot fill
staging or saturate the link while people are working:

```cron
# 20:00 every night: up to 300 GB, 4 files at a time, at most 40 MB/s
0 20 * * * cd $HOME/project-white-whale && WHALE_MIRROR_STAGING=/home/shared/8TB_HDD_03/sr320/whale/mirror-staging .venv/bin/whale mirror-fastq --jobs 4 --limit-rate 40 --max-gb 300 >> $HOME/whale-mirror.log 2>&1
```

Cron needs raven's SSH key for gannet to work without a passphrase, which it
does (BatchMode is used). Or run it by hand in `tmux`, which is the same
command without `--max-gb`.

Useful options:

| Option | Meaning |
|---|---|
| `--accession PRJ…` | Only these (repeatable), e.g. the next studies to be claimed |
| `--jobs N` | Files at once (default 4). ENA throughput per connection is low and variable, so parallel files help |
| `--limit-rate MB/s` | Total cap across all jobs |
| `--max-gb G` | Stop starting new files once G GB are queued in this run |
| `--no-verify-remote` | Skip the MD5 on gannet (rsync already checks what it transfers) |
| `--dest DIR` | Mirror into a local directory instead (a mounted share, or testing) |

## Using the mirror

A study claimed for raven or klone can stage its FASTQ from
`raw/<ACCESSION>/` on gannet (rsync or HTTPS) instead of ENA, and check it
against the study's `fastq_manifest.tsv` as usual (step 6 of
[adding a study](adding_a_study.md)).
