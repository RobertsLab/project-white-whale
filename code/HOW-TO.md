# How-To Guide: Downloading Oyster DNA Methylation Data

A step-by-step guide for **first-time users**. No prior bioinformatics
experience needed. Follow the steps in order; each one tells you what to type
and what you should see.

---

## What this does

This tool downloads publicly available DNA methylation sequencing data for the
Pacific oyster (*Crassostrea gigas* / *Magallana gigas*) from NCBI, the U.S.
national sequence database. You pick a "dataset" by name, and the tool figures
out which sequencing files belong to it and downloads them for you.

You do **not** need to know any accession numbers. The tool looks them up live
from NCBI every time, so you always get the real, current files.

---

## Before you start: a 2-minute checklist

| You need | How to check | If missing |
|----------|--------------|------------|
| **Python 3.6+** | `python3 --version` | Install from [python.org](https://www.python.org/downloads/) |
| **Internet** | Open any website | Connect to a network |
| **SRA Toolkit** | `fasterq-dump --version` | See Step 1 below |
| **Disk space** | A real download can be 10s–100s of GB | Use a drive with room |

> 💡 You can do **everything except the final download** without the SRA
> Toolkit. So feel free to explore first (Steps 2–3) and install it later.

---

## Step 1 — Install the download tool (one time only)

From inside the `code/` folder, run:

```bash
./install_dependencies.sh
```

This installs the **SRA Toolkit** (the program that actually fetches files). It
works on macOS (Homebrew), Linux (apt/yum), or anywhere with conda. If the
script can't auto-install, it prints a link to do it manually.

Check it worked:

```bash
fasterq-dump --version
```

You should see a version number. If you see "command not found", re-open your
terminal and try again, or follow the manual link the installer printed.

---

## Step 2 — See what's available

```bash
python3 download_methylation_data.py --list
```

This prints the available datasets, for example `wgbs_ph_ploidy`,
`medip_development`, `wgbs_pesto`. Each one bundles one or more NCBI
BioProjects. Note the **Dataset ID** of one you're interested in — you'll use it
in the next step.

---

## Step 3 — Preview before you download (always do this first)

A "dry run" shows you exactly what *would* happen, **without downloading
anything or using any disk space**:

```bash
python3 download_methylation_data.py --dataset wgbs_ph_ploidy --dry-run
```

You'll see lines like:

```
Found 24 matching run(s) in PRJNA682817
DRY RUN: would download 24 run(s) from PRJNA682817: SRR13207073, SRR13207072, ...
```

This tells you the tool reached NCBI and found real oyster methylation files.
Before listing runs, the tool checks each one's NCBI metadata and keeps only
*Crassostrea/Magallana gigas* bisulfite/MeDIP runs; anything off-target is
reported as `Skipped N off-target run(s)` and dropped. If a project has no
matching data you'll see `Found 0 matching run(s)` and it is simply skipped.
**This is normal and safe** — the tool never makes up files and never downloads
non-oyster or non-methylation data.

---

## Step 4 — Do a small test download

Before committing to a huge download, grab just a couple of files to confirm
your setup works end-to-end:

```bash
python3 download_methylation_data.py \
    --dataset wgbs_ph_ploidy \
    --bioproject PRJNA682817 \
    --max-runs 2
```

- `--bioproject` limits to one project.
- `--max-runs 2` limits to the first 2 files.

When it finishes you'll find compressed `.fastq.gz` files under
`methylation_data/wgbs_ph_ploidy/PRJNA682817/`. The tool also writes a
`runs.txt` listing every accession it found, so your download is reproducible.

---

## Step 5 — Download a full dataset

Once the test works, download everything in a dataset:

```bash
python3 download_methylation_data.py --dataset wgbs_ph_ploidy
```

Before it starts, the tool **checks your free disk space** and warns you if
there might not be enough. Large datasets can take **hours** (the WGBS datasets
are hundreds of GB each). You can safely stop (Ctrl-C) and re-run later — files
already downloaded are skipped.

To save somewhere with more room:

```bash
python3 download_methylation_data.py --dataset wgbs_ph_ploidy \
    --output-dir /Volumes/BigDrive/oyster_data
```

---

## Common situations

**"SRA Toolkit not found"** → You skipped Step 1. Datasets still preview with
`--dry-run`, but you must install the toolkit to download.

**"Could not reach NCBI"** → Network problem. Check your connection and retry;
the tool already retries a few times on its own.

**"Only N GB free…"** → Not enough disk space. Free some up, point
`--output-dir` at a bigger drive, or add `--force` to override the check at your
own risk.

**Downloads feel slow / rate-limited** → Tell NCBI who you are to get a higher
limit:

```bash
export NCBI_EMAIL="you@example.com"
# optional, from a free NCBI account: export NCBI_API_KEY="..."
```

---

## Generating a standalone script (optional)

If you'd rather run the downloads as a plain shell script (e.g. on a cluster):

```bash
python3 download_methylation_data.py --create-script \
    --datasets wgbs_ph_ploidy medip_development
```

This writes `methylation_data/download_script.sh`. It discovers the real runs
from NCBI when you run it — no fake accessions baked in. Test it with a small
limit first:

```bash
cd methylation_data
MAX_RUNS=2 ./download_script.sh   # download only 2 runs per project
./download_script.sh              # download everything
```

---

## Where things end up

```
methylation_data/
└── <dataset_id>/
    ├── dataset_info.json          # what this dataset is
    └── <BioProject>/
        ├── runs.txt               # the accessions that were found
        └── SRRxxxxxxx_1.fastq.gz  # the sequencing data
```

A full log of everything that happened is saved to
`download_methylation_data.log`.

---

## Getting help

- Run any command with `--help` to see all options.
- Detailed reference: [`USAGE.md`](USAGE.md).
- SRA Toolkit problems: <https://github.com/ncbi/sra-tools/issues>.
