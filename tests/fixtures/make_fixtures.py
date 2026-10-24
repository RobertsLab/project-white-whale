"""Regenerate the workflow fixtures: a small genome and bisulfite reads that map to it.

    python3 tests/fixtures/make_fixtures.py

Deterministic (fixed seed), so the committed files only change when this
script does. The reads are a directional WGBS library simulated from the
genome with known methylation, so a real (not stub) run of the workflow can
be checked against the truth (tests/check_real_run.py):

* every CpG is methylated with probability CPG_METHYLATION, on both strands;
* every other C is unmethylated, so bisulfite converts it to T;
* S1 is paired-end, S2 single-end, matching tests/fixtures/samplesheet.tsv.
"""
from __future__ import annotations

import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED = 20261024
CHROMS = {"chr1": 20_000, "chr2": 10_000}
GC = 0.36                 # about the oyster genome's
CPG_METHYLATION = 0.8
READ_LEN = 100
FRAGMENT = (180, 400)
ERROR_RATE = 0.002
N_READS = {"S1": 1500, "S2": 1500}   # pairs for S1, reads for S2

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def genome(rng: random.Random) -> dict[str, str]:
    weights = [(1 - GC) / 2, GC / 2, GC / 2, (1 - GC) / 2]
    return {name: "".join(rng.choices("ACGT", weights, k=n)) for name, n in CHROMS.items()}


def methylated_cpgs(seqs: dict[str, str], rng: random.Random) -> dict[str, set[int]]:
    """Positions i (0-based, top strand) of the C of each methylated CpG."""
    return {name: {i for i in range(len(s) - 1) if s[i:i + 2] == "CG" and rng.random() < CPG_METHYLATION}
            for name, s in seqs.items()}


def convert(strand: str, methylated: set[int], to_top) -> str:
    """Bisulfite-convert one strand: C -> T unless it is the C of a methylated CpG.

    `to_top(i)` maps index i on this strand to the top-strand position of the
    CpG's C, or None when the C is not in a CpG.
    """
    out = []
    for i, base in enumerate(strand):
        if base == "C":
            pos = to_top(i)
            out.append("C" if pos is not None and pos in methylated else "T")
        else:
            out.append(base)
    return "".join(out)


def fragment(seqs, meth, rng):
    """One converted fragment from a directional library: original top or bottom strand."""
    name = rng.choice(list(seqs))
    s = seqs[name]
    length = rng.randint(*FRAGMENT)
    start = rng.randrange(0, len(s) - length)
    top = s[start:start + length]
    m = meth[name]
    # CpG context comes from the genome, not the fragment, so a C at a fragment's
    # end is still judged by the base beyond it.
    if rng.random() < 0.5:   # original top: a C followed by G
        return convert(top, m, lambda i: start + i if s[start + i:start + i + 2] == "CG" else None)
    bottom = revcomp(top)    # original bottom: its CpG C sits opposite the top strand's CpG G
    end = start + length - 1

    def to_top(i):
        g = end - i          # top-strand position of this base's partner
        return g - 1 if g >= 1 and s[g - 1:g + 1] == "CG" else None
    return convert(bottom, m, to_top)


def with_errors(read: str, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < ERROR_RATE else b for b in read)


def fastq(name: str, seq: str) -> str:
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def write_gz(path: Path, text: str) -> None:
    # mtime=0 keeps the gzip bytes identical between runs.
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write(text.encode())


def main() -> None:
    rng = random.Random(SEED)
    seqs = genome(rng)
    meth = methylated_cpgs(seqs, rng)

    ref = HERE / "reference" / "TEST" / "genomic.fna"
    ref.parent.mkdir(parents=True, exist_ok=True)
    ref.write_text("".join(f">{n}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n"
                           for n, s in seqs.items()))

    reads = HERE / "reads" / "TEST2026_FIXTURE_WGBS"
    reads.mkdir(parents=True, exist_ok=True)
    r1, r2 = [], []
    for k in range(N_READS["S1"]):
        f = fragment(seqs, meth, rng)
        r1.append(fastq(f"S1.{k}/1", with_errors(f[:READ_LEN], rng)))
        r2.append(fastq(f"S1.{k}/2", with_errors(revcomp(f)[:READ_LEN], rng)))
    write_gz(reads / "S1_1.fastq.gz", "".join(r1))
    write_gz(reads / "S1_2.fastq.gz", "".join(r2))
    se = [fastq(f"S2.{k}", with_errors(fragment(seqs, meth, rng)[:READ_LEN], rng)) for k in range(N_READS["S2"])]
    write_gz(reads / "S2_1.fastq.gz", "".join(se))

    n_cpg = sum(s.count("CG") for s in seqs.values())
    print(f"{sum(CHROMS.values()):,} bp, {n_cpg} CpGs ({sum(map(len, meth.values()))} methylated); "
          f"{N_READS['S1']} pairs (S1), {N_READS['S2']} single reads (S2)")


if __name__ == "__main__":
    main()
