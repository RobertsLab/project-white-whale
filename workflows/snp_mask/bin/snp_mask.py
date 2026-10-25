#!/usr/bin/env python3
"""Build the C/T SNP mask (PLAN.md §6.1) from pooled bisulfite pileup sites.
Standard library only: it runs in the pinned python container, and
tests/test_snp_mask.py loads it from this file.

    snp_mask.py --genome genome.fa --params params.json [--min-alt 3] [--min-af 0.02] SITES...

SITES are pileup_sites.awk outputs (gzipped or plain): chrom, pos (1-based),
ref, alt bases, depth, alt reads, samples with alt. A site is a SNP when its
alt reads >= --min-alt and alt reads / depth >= --min-af. A reference CpG
is masked when a SNP sits on its C or its G: a C/T SNP at the C (or G/A at
the G) reads exactly like an unmethylated C, and any other change breaks the
CpG. Writes, in the current directory:

    snps.tsv.gz        every SNP: chrom, pos0, ref, alt, depth, alt_reads, alt_frac, samples_with_alt
    cpg_mask.bed.gz    masked CpGs, BED6: chrom, start (the C, 0-based), start + 2,
                       name (first SNP, e.g. 1234C>T), score (its alt reads, max 1000), "."
    mask_manifest.json params (mask id, inputs, thresholds), genome and output
                       checksums, and counts: SNPs, CpGs in the genome and masked,
                       per sequence, and the alt-read distribution
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

SNP_COLUMNS = ["chrom", "pos0", "ref", "alt", "depth", "alt_reads", "alt_frac", "samples_with_alt"]


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def read_sites(paths: list[Path], *, min_alt: int, min_af: float) -> dict[str, list[tuple]]:
    """chrom -> [(pos0, ref, alt, depth, alt_reads, af, samples)] passing the thresholds, sorted."""
    out: dict[str, list[tuple]] = {}
    for path in paths:
        with _open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 7:
                    continue
                depth, alt_reads = int(f[4]), int(f[5])
                af = alt_reads / depth if depth else 0.0
                if alt_reads >= min_alt and af >= min_af:
                    out.setdefault(f[0], []).append((int(f[1]) - 1, f[2], f[3], depth, alt_reads, af, int(f[6])))
    for rows in out.values():
        rows.sort()
    return out


def read_fasta(path: Path):
    """(name, sequence) per record, one at a time."""
    name, parts = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def count_cpgs(seq: str) -> int:
    return seq.upper().count("CG")


def mask_cpgs(seq: str, snps: list[tuple]) -> dict[int, tuple]:
    """C position of each masked CpG -> the first SNP that masks it."""
    up = seq.upper()
    masked: dict[int, tuple] = {}
    for snp in snps:
        p = snp[0]
        for c in (p, p - 1):      # the SNP is the CpG's C, or its G
            if 0 <= c and up[c:c + 2] == "CG" and c not in masked:
                masked[c] = snp
    return masked


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _gz_writer(path: Path):
    # No name or mtime in the header: the same mask gives the same bytes.
    raw = open(path, "wb")
    return raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0)


def build(genome: Path, sites: list[Path], *, min_alt: int, min_af: float, params: dict,
          out_dir: Path = Path(".")) -> dict:
    by_chrom = read_sites(sites, min_alt=min_alt, min_af=min_af)
    snps_path, bed_path = out_dir / "snps.tsv.gz", out_dir / "cpg_mask.bed.gz"
    per_seq, alt_hist = [], Counter()
    (snps_raw, snps_gz), (bed_raw, bed_gz) = _gz_writer(snps_path), _gz_writer(bed_path)
    snps_gz.write(("\t".join(SNP_COLUMNS) + "\n").encode())
    seen = set()
    for name, seq in read_fasta(genome):
        seen.add(name)
        snps = by_chrom.get(name, [])
        masked = mask_cpgs(seq, snps)
        for p0, ref, alt, depth, alt_reads, af, samples in snps:
            snps_gz.write(f"{name}\t{p0}\t{ref}\t{alt}\t{depth}\t{alt_reads}\t{af:.4f}\t{samples}\n".encode())
            alt_hist[min(alt_reads, 50)] += 1
        for c in sorted(masked):
            p0, ref, alt, _, alt_reads, _, _ = masked[c]
            bed_gz.write(f"{name}\t{c}\t{c + 2}\t{p0 + 1}{ref}>{alt}\t{min(alt_reads, 1000)}\t.\n".encode())
        n_cpg = count_cpgs(seq)
        per_seq.append({"sequence": name, "length": len(seq), "snps": len(snps), "cpgs": n_cpg,
                        "cpgs_masked": len(masked)})
    for gz, raw in ((snps_gz, snps_raw), (bed_gz, bed_raw)):
        gz.close()
        raw.close()
    unknown = sorted(set(by_chrom) - seen)
    if unknown:
        raise ValueError(f"sites on sequences absent from the genome: {unknown[:5]}")
    total = {k: sum(s[k] for s in per_seq) for k in ("length", "snps", "cpgs", "cpgs_masked")}
    manifest = {
        "params": {**params, "min_alt_reads": min_alt, "min_alt_frac": min_af},
        "genome": {"name": genome.name, "sha256": sha256(genome)},
        "outputs": {p.name: sha256(p) for p in (snps_path, bed_path)},
        "counts": {**total, "cpgs_masked_frac": round(total["cpgs_masked"] / total["cpgs"], 6) if total["cpgs"] else 0.0,
                   "snps_per_kb": round(1000 * total["snps"] / total["length"], 4) if total["length"] else 0.0},
        "alt_reads_hist": {str(k) if k < 50 else "50+": v for k, v in sorted(alt_hist.items())},
        "per_sequence": per_seq,
    }
    (out_dir / "mask_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--params", type=Path, required=True, help="JSON: mask id, samples, inputs (from the workflow).")
    ap.add_argument("--min-alt", type=int, default=3)
    ap.add_argument("--min-af", type=float, default=0.02)
    ap.add_argument("sites", nargs="*", type=Path)
    a = ap.parse_args(argv)
    m = build(a.genome, a.sites, min_alt=a.min_alt, min_af=a.min_af, params=json.loads(a.params.read_text()))
    c = m["counts"]
    print(f"{c['snps']:,} SNPs ({c['snps_per_kb']} per kb); {c['cpgs_masked']:,} of {c['cpgs']:,} CpGs masked "
          f"({100 * c['cpgs_masked_frac']:.2f}%)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
