"""Regenerate the enrichment-workflow fixtures: MeDIP-like reads from the TEST genome.

    python3 tests/fixtures/make_enrichment_fixtures.py

Deterministic (fixed seed). Uses the genome and features written for the
methylation fixtures (reference/TEST/), and simulates a mosaic methylome like
an oyster's: CpGs inside gene bodies (features/genes.bed) are methylated with
probability 0.9, CpGs elsewhere with 0.05. A fragment is captured with
probability proportional to its methylated CpGs (saturating at CAPTURE_AT),
plus a small uniform background, as an anti-5mC pull-down would. Reads are
not bisulfite-converted. About 5% of fragments are sequenced twice, so the
duplicate step has something to remove.

* E1 is paired-end, E2 single-end (tests/fixtures/samplesheet_enrichment.tsv).

tests/check_enrichment_run.py checks a real run against this: gene bodies
enriched over intergenic DNA, captured fragments CpG-rich (relH > 1).
"""
from __future__ import annotations

import random
from pathlib import Path

from make_fixtures import fastq, revcomp, with_errors, write_gz

HERE = Path(__file__).resolve().parent
SEED = 20261121
GENE_METH, OTHER_METH = 0.9, 0.05
CAPTURE_AT = 4            # methylated CpGs at which capture probability reaches 1
BACKGROUND = 0.03
DUPLICATES = 0.05
READ_LEN = 100
FRAGMENT = (160, 280)
N_FRAGMENTS = {"E1": 2500, "E2": 2500}


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def genes(path: Path) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return out


def methylation(seqs, gene_iv, rng) -> dict[str, set[int]]:
    out = {}
    for name, s in seqs.items():
        ivs = gene_iv.get(name, [])
        out[name] = {i for i in range(len(s) - 1) if s[i:i + 2] == "CG"
                     and rng.random() < (GENE_METH if any(a <= i < b for a, b in ivs) else OTHER_METH)}
    return out


def captured(seqs, meth, rng) -> str:
    """One pulled-down fragment, as top-strand sequence."""
    names = list(seqs)
    weights = [len(seqs[n]) for n in names]
    while True:
        name = rng.choices(names, weights)[0]
        s = seqs[name]
        length = rng.randint(*FRAGMENT)
        start = rng.randrange(0, len(s) - length)
        m = sum(1 for i in range(start, start + length - 1) if i in meth[name])
        if rng.random() < min(1.0, m / CAPTURE_AT) or rng.random() < BACKGROUND:
            frag = s[start:start + length]
            return frag if rng.random() < 0.5 else revcomp(frag)


def main() -> None:
    rng = random.Random(SEED)
    ref = HERE / "reference" / "TEST"
    seqs = read_fasta(ref / "genomic.fna")
    meth = methylation(seqs, genes(ref / "features" / "genes.bed"), rng)
    out = HERE / "reads" / "TEST2026_FIXTURE_MEDIP"
    out.mkdir(parents=True, exist_ok=True)

    frags = []
    for _ in range(N_FRAGMENTS["E1"]):
        f = captured(seqs, meth, rng)
        frags.append(f)
        if rng.random() < DUPLICATES:
            frags.append(f)
    r1 = [fastq(f"E1.{k}/1", with_errors(f[:READ_LEN], rng)) for k, f in enumerate(frags)]
    r2 = [fastq(f"E1.{k}/2", with_errors(revcomp(f)[:READ_LEN], rng)) for k, f in enumerate(frags)]
    write_gz(out / "E1_1.fastq.gz", "".join(r1))
    write_gz(out / "E1_2.fastq.gz", "".join(r2))

    se = []
    for _ in range(N_FRAGMENTS["E2"]):
        f = captured(seqs, meth, rng)
        se.append(f)
        if rng.random() < DUPLICATES:
            se.append(f)
    write_gz(out / "E2_1.fastq.gz", "".join(fastq(f"E2.{k}", with_errors(f[:READ_LEN], rng))
                                            for k, f in enumerate(se)))

    (HERE / "samplesheet_enrichment.tsv").write_text(
        "sample_id\trun_accession\tfastq_1\tfastq_2\n"
        "E1\tSRR0000011\tE1_1.fastq.gz\tE1_2.fastq.gz\n"
        "E2\tSRR0000012\tE2_1.fastq.gz\t\n")
    print(f"{len(frags)} pairs (E1), {len(se)} single reads (E2)")


if __name__ == "__main__":
    main()
