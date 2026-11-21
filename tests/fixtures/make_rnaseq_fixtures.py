"""Regenerate the RNA-seq workflow fixtures: transcripts, tx2gene and reads.

    python3 tests/fixtures/make_rnaseq_fixtures.py

Deterministic (fixed seed). Builds on the methylation fixtures' genome
(reference/TEST/genomic.fna, written by make_fixtures.py) and its four
hand-written genes (reference/TEST/features/exons.bed):

* ``reference/TEST/rna.fna``     the genes' transcripts, spliced from the genome
                                 exons, named and titled like NCBI's rna.fna; G1
                                 has two isoforms that differ in exon 1's end;
* ``reference/TEST/tx2gene.tsv`` transcript -> gene, as `whale fetch-reference
                                 --rnaseq` writes it (no header);
* ``reads/TEST2026_FIXTURE_WGBS/rnaseq/``  a stranded (dUTP, read 1 antisense)
                                 library per sample: R1 paired-end, R2 single-end;
* ``rnaseq/samplesheet.tsv``     the sample sheet the test profile reads;
* ``rnaseq/truth.tsv``           fragments simulated per sample and gene and
                                 transcript, plus two kinds that belong to no
                                 transcript: ``genomic`` (inside introns or
                                 intergenic DNA, matching nothing) and ``premrna``
                                 (unspliced, across an exon-intron boundary, which
                                 Salmon's genome decoys should take; ``premrna_<gene>``
                                 splits them by the gene whose boundary they cross).

tests/check_rnaseq_run.py compares a real (Docker) run's quantification with
truth.tsv.
"""
from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED = 20261121
READ_LEN = 75
FRAGMENT = (150, 300)
ERROR_RATE = 0.002
GENOMIC_FRACTION = 0.10   # fragments from introns / intergenic DNA: match no transcript
PREMRNA_FRACTION = 0.10   # unspliced fragments across an exon-intron boundary: the decoys' job
STUDY = "TEST2026_FIXTURE_WGBS"

# transcript -> (gene, strand, exons as 0-based half-open intervals), matching exons.bed.
TRANSCRIPTS = {
    "XM_900001.1": ("G1", "+", "chr1", [(2000, 3000), (6000, 8000)]),
    "XM_900002.1": ("G1", "+", "chr1", [(2000, 3500), (6000, 8000)]),
    "XM_900003.1": ("G2", "-", "chr1", [(12000, 13000), (14000, 15000)]),
    "XM_900004.1": ("G3", "+", "chr2", [(1000, 2000), (5000, 6000)]),
    "XR_900005.1": ("G4", "+", "chr2", [(9000, 9500)]),
}
TITLES = {
    "G1": "test gene one (G1), transcript variant X{v}, mRNA",
    "G2": "test gene two (G2), mRNA",
    "G3": "test gene three (G3), mRNA",
    "G4": "uncharacterized LOC900004 (G4), ncRNA",
}
# Not exonic in any transcript: introns of G1, G2, G3 and intergenic DNA.
GENOMIC_REGIONS = [("chr1", 3600, 5900), ("chr1", 8200, 11800), ("chr1", 15200, 19800),
                   ("chr2", 2100, 4900), ("chr2", 6200, 8800)]

# Exon-intron boundaries of the multi-exon genes, with the gene's strand.
BOUNDARIES = [("chr1", 3500, "+", "G1"), ("chr1", 6000, "+", "G1"), ("chr1", 13000, "-", "G2"),
              ("chr1", 14000, "-", "G2"), ("chr2", 2000, "+", "G3"), ("chr2", 5000, "+", "G3")]
MIN_OVERLAP = 35          # bases of a pre-mRNA fragment on each side of its boundary

# Share of transcript-derived fragments per transcript.
SAMPLES = {
    "R1": {"paired": True, "fragments": 3000,
           "abundance": {"XM_900001.1": 0.30, "XM_900002.1": 0.10, "XM_900003.1": 0.35,
                         "XM_900004.1": 0.25, "XR_900005.1": 0.0}},
    "R2": {"paired": False, "fragments": 3000,
           "abundance": {"XM_900001.1": 0.10, "XM_900002.1": 0.10, "XM_900003.1": 0.20,
                         "XM_900004.1": 0.50, "XR_900005.1": 0.10}},
}

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def transcript_seqs(genome: dict[str, str]) -> dict[str, str]:
    out = {}
    for tx, (_, strand, chrom, exons) in TRANSCRIPTS.items():
        s = "".join(genome[chrom][a:b] for a, b in exons)
        out[tx] = revcomp(s) if strand == "-" else s
    return out


def mutate(seq: str, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < ERROR_RATE else b for b in seq)


def fragment(seq: str, rng: random.Random) -> str:
    n = min(rng.randint(*FRAGMENT), len(seq))
    start = rng.randint(0, len(seq) - n)
    return seq[start:start + n]


def fastq(name: str, seq: str) -> str:
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def simulate(sample: str, spec: dict, txs: dict[str, str], genome: dict[str, str],
             rng: random.Random) -> tuple[list[str], list[str], dict[str, int]]:
    n_genomic = round(spec["fragments"] * GENOMIC_FRACTION)
    n_premrna = round(spec["fragments"] * PREMRNA_FRACTION)
    n_tx = spec["fragments"] - n_genomic - n_premrna
    names = list(spec["abundance"])
    weights = [spec["abundance"][t] for t in names]
    truth = {t: 0 for t in names}
    r1, r2 = [], []
    for i in range(spec["fragments"]):
        if i < n_tx:
            tx = rng.choices(names, weights)[0]
            truth[tx] += 1
            frag = fragment(txs[tx], rng)
        elif i < n_tx + n_premrna:
            chrom, b, strand, gene = rng.choice(BOUNDARIES)
            truth[f"premrna_{gene}"] = truth.get(f"premrna_{gene}", 0) + 1
            n = rng.randint(*FRAGMENT)
            start = b - rng.randint(MIN_OVERLAP, n - MIN_OVERLAP)
            frag = genome[chrom][start:start + n]
            if strand == "-":
                frag = revcomp(frag)
        else:
            chrom, a, b = rng.choice(GENOMIC_REGIONS)
            frag = fragment(genome[chrom][a:b], rng)
            if rng.random() < 0.5:
                frag = revcomp(frag)
        # dUTP: read 1 is the reverse complement of the RNA, read 2 the RNA itself.
        name = f"{sample}_{i + 1}"
        r1.append(fastq(f"{name}/1", mutate(revcomp(frag)[:READ_LEN], rng)))
        if spec["paired"]:
            r2.append(fastq(f"{name}/2", mutate(frag[:READ_LEN], rng)))
    order = list(range(len(r1)))
    rng.shuffle(order)
    r1 = [r1[i] for i in order]
    r2 = [r2[i] for i in order] if r2 else []
    truth["genomic"] = n_genomic
    truth["premrna"] = n_premrna
    return r1, r2, truth


def write_gz(path: Path, records: list[str]) -> None:
    with open(path, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write("".join(records).encode())


def main() -> None:
    rng = random.Random(SEED)
    ref = HERE / "reference" / "TEST"
    genome = read_fasta(ref / "genomic.fna")
    txs = transcript_seqs(genome)

    variant = 0
    with open(ref / "rna.fna", "w") as fh:
        for tx, seq in txs.items():
            gene = TRANSCRIPTS[tx][0]
            if gene == "G1":
                variant += 1
            title = TITLES[gene].format(v=variant)
            fh.write(f">{tx} PREDICTED: Magallana gigas {title}\n")
            fh.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
    (ref / "tx2gene.tsv").write_text("".join(f"{tx}\t{g[0]}\n" for tx, g in sorted(TRANSCRIPTS.items())))

    reads = HERE / "reads" / STUDY / "rnaseq"
    reads.mkdir(parents=True, exist_ok=True)
    sheet, truth_rows = [], []
    for sample, spec in SAMPLES.items():
        r1, r2, truth = simulate(sample, spec, txs, genome, rng)
        write_gz(reads / f"{sample}_1.fastq.gz", r1)
        if r2:
            write_gz(reads / f"{sample}_2.fastq.gz", r2)
        sheet.append({"sample_id": sample, "run_accession": f"SRR99{sample[1:]:0>5}",
                      "fastq_1": f"{sample}_1.fastq.gz", "fastq_2": f"{sample}_2.fastq.gz" if r2 else "",
                      "condition": "control" if sample == "R1" else "treated"})
        by_gene: dict[str, int] = {}
        for tx, n in truth.items():
            if tx.startswith("premrna_"):
                continue
            gene = TRANSCRIPTS[tx][0] if tx in TRANSCRIPTS else tx
            by_gene[gene] = by_gene.get(gene, 0) + n
        truth_rows += [{"sample_id": sample, "feature": k, "fragments": v} for k, v in sorted(by_gene.items())]
        truth_rows += [{"sample_id": sample, "feature": k, "fragments": v}
                       for k, v in sorted(truth.items()) if k in TRANSCRIPTS or k.startswith("premrna_")]

    out = HERE / "rnaseq"
    out.mkdir(exist_ok=True)
    for name, rows in (("samplesheet.tsv", sheet), ("truth.tsv", truth_rows)):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
