#!/usr/bin/env python3
"""Per-sample CpG methylation of every gene's body, promoter, exons and introns
(PLAN.md §5 SUMMARIES, §6.2 gene_features). Standard library only: it runs in
the pinned python container, and tests/test_methylome_features.py loads it
from this file.

    methylome_features.py SAMPLE COV_GZ FEATURES_DIR OUT_TSV_GZ

COV_GZ is one strand-merged CpG file from coverage2cytosine (chrom, start,
end, % methylated, count methylated, count unmethylated; covered CpGs only).
FEATURES_DIR is the reference bundle's features/ (genes.bed, promoters.bed,
exons.bed, introns.bed; BED6, name = NCBI GeneID).

The output has one row per gene in genes.bed and feature class, in
genes.bed order, so every sample of every study has the same rows; a gene
with no covered CpG in a class has zeros and empty percentages. Counts are
kept so that coverage cut-offs can be chosen downstream:

    feature            gene_body, promoter, exon or intron
    gene_id            NCBI GeneID
    cpgs               CpGs with at least one call
    cpgs_5x            CpGs with at least five calls
    meth, unmeth       methylated and unmethylated calls, summed over the CpGs
    meth_pct           100 * meth / (meth + unmeth)
    mean_meth_5x_pct   mean of the per-CpG methylation of the CpGs at >= 5x

A CpG in two overlapping exons (or introns) of one gene is counted once for
that gene; a CpG in two genes counts for both.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

CLASSES = [("gene_body", "genes.bed"), ("promoter", "promoters.bed"),
           ("exon", "exons.bed"), ("intron", "introns.bed")]
COLUMNS = ["feature", "gene_id", "cpgs", "cpgs_5x", "meth", "unmeth", "meth_pct", "mean_meth_5x_pct"]


def read_bed(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """chrom -> [(start, end, name)] sorted by start."""
    out: dict[str, list[tuple[int, int, str]]] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and not line.startswith(("#", "track")):
                out.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3]))
    for rows in out.values():
        rows.sort()
    return out


def gene_universe(features_dir: Path) -> list[str]:
    seen, order = set(), []
    with open(features_dir / "genes.bed") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 4 and f[3] not in seen:
                seen.add(f[3])
                order.append(f[3])
    return order


class Sweep:
    """The features of one class covering ascending positions along one chromosome."""

    def __init__(self, intervals: list[tuple[int, int, str]]):
        self.iv, self.next, self.active = intervals, 0, []

    def at(self, pos: int) -> set[str]:
        iv = self.iv
        while self.next < len(iv) and iv[self.next][0] <= pos:
            self.active.append(iv[self.next])
            self.next += 1
        if self.active:
            self.active = [a for a in self.active if a[1] > pos]
        return {a[2] for a in self.active}


def summarize(cov_gz: Path, features_dir: Path) -> dict[tuple[str, str], list]:
    """(class, gene) -> [cpgs, cpgs_5x, meth, unmeth, sum of per-CpG fractions at >= 5x]."""
    beds = {cls: read_bed(features_dir / name) for cls, name in CLASSES}
    acc: dict[tuple[str, str], list] = {}
    chrom, sweeps, last, done = None, {}, -1, set()
    with gzip.open(cov_gz, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 6:
                continue
            pos, m, u = int(f[1]), int(f[4]), int(f[5])
            if f[0] != chrom:
                if f[0] in done:
                    raise ValueError(f"{cov_gz}: {f[0]} is not contiguous; the CpG file must be sorted")
                if chrom is not None:
                    done.add(chrom)
                chrom, last = f[0], -1
                sweeps = {cls: Sweep(beds[cls].get(chrom, [])) for cls, _ in CLASSES}
            if pos < last:
                raise ValueError(f"{cov_gz}: {chrom}:{pos} comes after {last}; the CpG file must be sorted")
            last = pos
            d = m + u
            if d == 0:
                continue
            frac5 = m / d if d >= 5 else None
            for cls, _ in CLASSES:
                for gene in sweeps[cls].at(pos):
                    a = acc.get((cls, gene))
                    if a is None:
                        a = acc[(cls, gene)] = [0, 0, 0, 0, 0.0]
                    a[0] += 1
                    a[2] += m
                    a[3] += u
                    if frac5 is not None:
                        a[1] += 1
                        a[4] += frac5
    return acc


def rows(acc: dict[tuple[str, str], list], genes: list[str]):
    for cls, _ in CLASSES:
        for gene in genes:
            n, n5, m, u, s5 = acc.get((cls, gene), (0, 0, 0, 0, 0.0))
            yield [cls, gene, n, n5, m, u,
                   f"{100 * m / (m + u):.3f}" if m + u else "",
                   f"{100 * s5 / n5:.3f}" if n5 else ""]


def write(out: Path, acc: dict, genes: list[str]) -> None:
    # mtime=0: the same calls give byte-identical files (and the same SHA-256).
    with open(out, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(("\t".join(COLUMNS) + "\n").encode())
        for r in rows(acc, genes):
            gz.write(("\t".join(map(str, r)) + "\n").encode())


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    sample, cov, features, out = argv[0], Path(argv[1]), Path(argv[2]), Path(argv[3])
    genes = gene_universe(features)
    acc = summarize(cov, features)
    write(out, acc, genes)
    covered = sum(1 for g in genes if acc.get(("gene_body", g), [0])[0])
    print(f"{sample}: {covered} of {len(genes)} gene bodies have a covered CpG", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
