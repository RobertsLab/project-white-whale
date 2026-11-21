#!/usr/bin/env python3
"""Per-sample enrichment of a MeDIP / MBD-seq library in 1 kb tiles, gene bodies
and promoters, with the library's enrichment QC (docs/enrichment.md).
Standard library only: it runs in the pinned python container, and
tests/test_enrichment_counts.py loads it from this file.

    enrichment_counts.py count SAMPLE FRAGMENTS_TSV_GZ GENOME_FASTA FEATURES_DIR OUTDIR \\
        [--window 1000] [--fragment-length 200] [--min-cpgs 5] [--seed 1]
    enrichment_counts.py qc-table --samplesheet SHEET --out qc_summary.tsv [thresholds] REPORT...

FRAGMENTS_TSV_GZ is `samtools view` of the deduplicated BAM, cut to FLAG,
RNAME, POS, CIGAR, TLEN (modules/enrichment/dedup.nf). A properly paired
fragment is taken from the mate with TLEN > 0 and spans the template; a
single-end read is extended in its own direction to --fragment-length (as
MEDIPS' `extend`). Each fragment is counted once, in the tile or feature
holding its midpoint.

Normalisation (recorded in run_manifest.json as `enrichment_method`):

    cpm      count * 1e6 / N, N = fragments counted in the whole library
    obs_exp  (count / N) / (cpgs / C), C = CpGs in the genome: reads per CpG
             in the region over reads per CpG in the genome. 1 = as enriched
             as the average CpG; empty where the region has < --min-cpgs CpGs.

obs_exp corrects for library size and for the region's CpG count, which is
what drives capture by an anti-5mC antibody or an MBD protein. It does not
convert enrichment into a methylation proportion: that needs a calibration
(MEDIPS' coupling factor, QSEA's sigmoid) fitted on a genome-wide CpG-density
curve that assumes mammalian, mostly methylated genomes, which an oyster's
mosaic methylome is not.

`count` writes, under OUTDIR:

    SAMPLE.tiles_1kb.tsv.gz             chrom, start, end, cpgs, count, cpm, obs_exp
    SAMPLE.enrichment_features.tsv.gz   feature, gene_id, cpgs, count, cpm, obs_exp
                                        (every gene in genes.bed order; gene_body, promoter)
    SAMPLE.enrichment_qc.tsv            one row, the QC columns in QC_COLUMNS
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import random
import re
import sys
from array import array
from pathlib import Path

CIGAR_REF = re.compile(r"(\d+)([MDN=X])")
CLASSES = [("gene_body", "genes.bed"), ("promoter", "promoters.bed")]
QC_COLUMNS = ["sample_id", "fragments", "fragment_len_median", "genome_cpgs", "cpgs_1x_pct", "cpgs_5x_pct",
              "relh", "goge", "saturation_r", "saturation_est", "tiles_with_reads_pct",
              "gene_body_obs_exp", "intergenic_obs_exp"]
COMPOSITION_MAX = 2_000_000   # fragments read for relH/GoGe (every k-th, so the sample is spread out)


# --------------------------------------------------------------------------- inputs

def read_fasta(path: Path) -> dict[str, str]:
    """name -> upper-case sequence (NCBI genomes are soft-masked)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, name, chunks = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def cpg_positions(seq: str) -> array:
    """0-based positions of the C of every CpG (top strand; strands merged)."""
    out = array("I")
    i = seq.find("CG")
    while i != -1:
        out.append(i)
        i = seq.find("CG", i + 1)
    return out


def read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else "."))
    return rows


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def ref_length(cigar: str) -> int:
    return sum(int(n) for n, op in CIGAR_REF.findall(cigar))


def fragment(fields: list[str], frag_len: int, chrom_len: int) -> tuple[int, int] | None:
    """(start, end) of the fragment a SAM record stands for, or None to skip the record."""
    flag, pos, cigar, tlen = int(fields[0]), int(fields[2]) - 1, fields[3], int(fields[4])
    if flag & 1:                                  # paired: one record per template
        if tlen <= 0:
            return None
        return pos, min(chrom_len, pos + tlen)
    span = ref_length(cigar)
    length = max(frag_len, span)
    if flag & 16:
        end = min(chrom_len, pos + span)
        return max(0, end - length), end
    return pos, min(chrom_len, pos + length)


def read_fragments(path: Path, genome: dict[str, str], frag_len: int):
    """chrom -> (starts, ends) arrays, in file order."""
    out: dict[str, tuple[array, array]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[1] not in genome:
                continue
            frag = fragment(f, frag_len, len(genome[f[1]]))
            if frag is None or frag[1] <= frag[0]:
                continue
            s, e = out.setdefault(f[1], (array("I"), array("I")))
            s.append(frag[0])
            e.append(frag[1])
    return out


# --------------------------------------------------------------------------- statistics

def pearson(x: list[float], y: list[float]) -> float | None:
    n = len(x)
    if n < 3:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else None


def obs_exp(count: int, cpgs: int, n: int, genome_cpgs: int, min_cpgs: int) -> float | None:
    if cpgs < min_cpgs or n == 0 or genome_cpgs == 0:
        return None
    return (count / n) / (cpgs / genome_cpgs)


def _fmt(v, digits=4) -> str:
    return "" if v is None else f"{v:.{digits}f}"


# --------------------------------------------------------------------------- count

def count(sample: str, fragments_gz: Path, genome_fa: Path, features_dir: Path, outdir: Path, *,
          window: int = 1000, frag_len: int = 200, min_cpgs: int = 5, seed: int = 1) -> dict:
    genome = read_fasta(genome_fa)
    frags = read_fragments(fragments_gz, genome, frag_len)
    n = sum(len(s) for s, _ in frags.values())
    rng = random.Random(seed)

    cpgs = {c: cpg_positions(s) for c, s in genome.items()}
    genome_cpgs = sum(len(p) for p in cpgs.values())
    genome_len = sum(len(s) for s in genome.values())
    genome_c = sum(s.count("C") for s in genome.values())
    genome_g = sum(s.count("G") for s in genome.values())

    mids: dict[str, list[int]] = {}
    tiles: dict[str, tuple[array, array, array, array]] = {}   # cpgs, count, half A, half B
    covered_1x = covered_5x = 0
    lengths: list[int] = []
    step = max(1, math.ceil(n / COMPOSITION_MAX))
    comp_len = comp_c = comp_g = comp_cg = 0
    k = 0
    for chrom, seq in genome.items():
        n_tiles = math.ceil(len(seq) / window)
        tc, tn, ta, tb = (array("I", [0]) * n_tiles for _ in range(4))
        for p in cpgs[chrom]:
            tc[p // window] += 1
        starts, ends = frags.get(chrom, (array("I"), array("I")))
        m = []
        for s, e in zip(starts, ends):
            mid = (s + e) // 2
            m.append(mid)
            t = mid // window
            tn[t] += 1
            (ta if rng.random() < 0.5 else tb)[t] += 1
            if k % step == 0:
                frag_seq = seq[s:e]
                comp_len += len(frag_seq)
                comp_c += frag_seq.count("C")
                comp_g += frag_seq.count("G")
                comp_cg += frag_seq.count("CG")
                if len(lengths) < 100_000:
                    lengths.append(e - s)
            k += 1
        m.sort()
        mids[chrom] = m
        tiles[chrom] = (tc, tn, ta, tb)
        # CpG coverage by whole fragments: starts <= p minus ends <= p.
        ss, es = sorted(starts), sorted(ends)
        for p in cpgs[chrom]:
            depth = bisect.bisect_right(ss, p) - bisect.bisect_right(es, p)
            covered_1x += depth >= 1
            covered_5x += depth >= 5

    outdir.mkdir(parents=True, exist_ok=True)
    # Tiles
    half_a, half_b = [], []
    with_reads = with_cpgs = 0
    with gzip.open(outdir / f"{sample}.tiles_1kb.tsv.gz", "wt", compresslevel=6) as out:
        out.write("chrom\tstart\tend\tcpgs\tcount\tcpm\tobs_exp\n")
        for chrom, seq in genome.items():
            tc, tn, ta, tb = tiles[chrom]
            for i in range(len(tc)):
                start, end = i * window, min(len(seq), (i + 1) * window)
                c = tn[i]
                cpm = c * 1e6 / n if n else None
                out.write(f"{chrom}\t{start}\t{end}\t{tc[i]}\t{c}\t{_fmt(cpm, 3)}\t"
                          f"{_fmt(obs_exp(c, tc[i], n, genome_cpgs, min_cpgs))}\n")
                if tc[i] >= 1:
                    with_cpgs += 1
                    with_reads += c > 0
                    half_a.append(math.log1p(ta[i]))
                    half_b.append(math.log1p(tb[i]))

    # Genes: gene bodies and promoters, every gene in genes.bed order.
    genes_bed = read_bed(features_dir / "genes.bed")
    order = list(dict.fromkeys(g for *_, g in genes_bed))
    per_class: dict[str, dict[str, tuple[int, int]]] = {}
    for cls, fname in CLASSES:
        by_gene: dict[str, dict[str, list[tuple[int, int]]]] = {}
        path = features_dir / fname
        for chrom, s, e, g in read_bed(path) if path.exists() else []:
            by_gene.setdefault(g, {}).setdefault(chrom, []).append((s, e))
        vals = {}
        for g, chroms in by_gene.items():
            nc = nf = 0
            for chrom, ivs in chroms.items():
                m, p = mids.get(chrom, []), cpgs.get(chrom, array("I"))
                for s, e in merge(ivs):
                    nf += bisect.bisect_left(m, e) - bisect.bisect_left(m, s)
                    nc += bisect.bisect_left(p, e) - bisect.bisect_left(p, s)
            vals[g] = (nc, nf)
        per_class[cls] = vals
    with gzip.open(outdir / f"{sample}.enrichment_features.tsv.gz", "wt") as out:
        out.write("feature\tgene_id\tcpgs\tcount\tcpm\tobs_exp\n")
        for cls, _ in CLASSES:
            for g in order:
                nc, nf = per_class[cls].get(g, (0, 0))
                cpm = nf * 1e6 / n if n else None
                out.write(f"{cls}\t{g}\t{nc}\t{nf}\t{_fmt(cpm, 3)}\t"
                          f"{_fmt(obs_exp(nf, nc, n, genome_cpgs, min_cpgs))}\n")

    # Region classes for QC: all gene bodies vs intergenic, as unions.
    def region(fname: str) -> float | None:
        path = features_dir / fname
        if not path.exists():
            return None
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, s, e, _ in read_bed(path):
            by_chrom.setdefault(chrom, []).append((s, e))
        nc = nf = 0
        for chrom, ivs in by_chrom.items():
            m, p = mids.get(chrom, []), cpgs.get(chrom, array("I"))
            for s, e in merge(ivs):
                nf += bisect.bisect_left(m, e) - bisect.bisect_left(m, s)
                nc += bisect.bisect_left(p, e) - bisect.bisect_left(p, s)
        return obs_exp(nf, nc, n, genome_cpgs, 1)

    relh = goge = None
    if comp_len and genome_len:
        relh = (comp_cg / comp_len) / (genome_cpgs / genome_len)
        if comp_c and comp_g and genome_c and genome_g:
            goge = (comp_cg * comp_len / (comp_c * comp_g)) / (genome_cpgs * genome_len / (genome_c * genome_g))
    r = pearson(half_a, half_b)
    lengths.sort()
    qc = {
        "sample_id": sample,
        "fragments": n,
        "fragment_len_median": lengths[len(lengths) // 2] if lengths else "",
        "genome_cpgs": genome_cpgs,
        "cpgs_1x_pct": _fmt(100 * covered_1x / genome_cpgs if genome_cpgs else None, 2),
        "cpgs_5x_pct": _fmt(100 * covered_5x / genome_cpgs if genome_cpgs else None, 2),
        "relh": _fmt(relh, 3),
        "goge": _fmt(goge, 3),
        # MEDIPS.saturation's idea: two random halves of the library; Spearman-Brown for the whole.
        "saturation_r": _fmt(r, 3),
        "saturation_est": _fmt(2 * r / (1 + r) if r is not None and r > -1 else None, 3),
        "tiles_with_reads_pct": _fmt(100 * with_reads / with_cpgs if with_cpgs else None, 2),
        "gene_body_obs_exp": _fmt(region("genes.bed"), 3),
        "intergenic_obs_exp": _fmt(region("intergenic.bed"), 3),
    }
    with open(outdir / f"{sample}.enrichment_qc.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=QC_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(qc)
    return qc


# --------------------------------------------------------------------------- QC table

BOWTIE2_READS = re.compile(r"^(\d+) reads; of these:", re.M)
BOWTIE2_RATE = re.compile(r"^([\d.]+)% overall alignment rate", re.M)


def parse_bowtie2(text: str) -> dict:
    reads, rate = BOWTIE2_READS.search(text), BOWTIE2_RATE.search(text)
    return {"reads": int(reads.group(1)) if reads else None, "mapping_pct": float(rate.group(1)) if rate else None}


def parse_markdup(text: str) -> dict:
    """samtools markdup -f stats: duplicates among the reads it examined."""
    vals = {}
    for line in text.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            try:
                vals[k.strip()] = int(v.strip())
            except ValueError:
                pass
    examined, dups = vals.get("EXAMINED"), vals.get("DUPLICATE TOTAL")
    pct = 100 * dups / examined if examined and dups is not None else None
    return {"duplication_pct": pct, "library_size_est": vals.get("ESTIMATED_LIBRARY_SIZE")}


def _float(v) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


TABLE_COLUMNS = ["sample_id", "condition", "reads", "mapping_pct", "duplication_pct", "fragments",
                 "fragment_len_median", "relh", "goge", "cpgs_1x_pct", "cpgs_5x_pct", "saturation_r",
                 "saturation_est", "tiles_with_reads_pct", "gene_body_obs_exp", "intergenic_obs_exp",
                 "qc_pass", "qc_reasons"]


def qc_table(reports: list[Path], samplesheet: Path, out: Path, *, min_mapping: float = 0.0,
             min_fragments: int = 0, min_relh: float = 0.0) -> list[dict]:
    """One row per library in the sample sheet; qc_pass false (with reasons) below any threshold."""
    with open(samplesheet, newline="") as fh:
        sheet = list(csv.DictReader(fh, delimiter="\t"))
    by: dict[str, dict] = {r["sample_id"]: {} for r in sheet}
    for path in reports:
        name = path.name
        for suffix, kind in ((".bowtie2.log", "bt2"), (".markdup.txt", "dup"), (".enrichment_qc.tsv", "enr")):
            if name.endswith(suffix):
                sid = name[: -len(suffix)]
                if sid not in by:
                    continue
                text = path.read_text()
                if kind == "bt2":
                    by[sid].update(parse_bowtie2(text))
                elif kind == "dup":
                    by[sid].update(parse_markdup(text))
                else:
                    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
                    if rows:
                        by[sid].update({k: v for k, v in rows[0].items() if k != "sample_id"})
    rows = []
    for r in sheet:
        sid, v = r["sample_id"], by[r["sample_id"]]
        reasons = []
        mapping, frags, relh = _float(v.get("mapping_pct")), _float(v.get("fragments")), _float(v.get("relh"))
        if mapping is None or mapping < min_mapping:
            reasons.append(f"mapping {mapping}% < {min_mapping}%")
        if frags is None or frags < min_fragments:
            reasons.append(f"fragments {frags and int(frags)} < {min_fragments}")
        if relh is None or relh < min_relh:
            reasons.append(f"CpG enrichment relH {relh} < {min_relh}")
        dup = v.get("duplication_pct")
        rows.append({
            "sample_id": sid, "condition": r.get("condition", ""),
            "reads": v.get("reads") if v.get("reads") is not None else "",
            "mapping_pct": _fmt(mapping, 2), "duplication_pct": _fmt(dup, 2),
            **{c: v.get(c, "") for c in TABLE_COLUMNS[5:16]},
            "qc_pass": "false" if reasons else "true", "qc_reasons": "; ".join(reasons),
        })
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TABLE_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------------------- CLI

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("count")
    c.add_argument("sample")
    c.add_argument("fragments", type=Path)
    c.add_argument("genome", type=Path)
    c.add_argument("features", type=Path)
    c.add_argument("outdir", type=Path)
    c.add_argument("--window", type=int, default=1000)
    c.add_argument("--fragment-length", type=int, default=200)
    c.add_argument("--min-cpgs", type=int, default=5)
    c.add_argument("--seed", type=int, default=1)
    q = sub.add_parser("qc-table")
    q.add_argument("--samplesheet", type=Path, required=True)
    q.add_argument("--out", type=Path, required=True)
    q.add_argument("--min-mapping", type=float, default=0.0)
    q.add_argument("--min-fragments", type=int, default=0)
    q.add_argument("--min-relh", type=float, default=0.0)
    q.add_argument("reports", type=Path, nargs="*")
    a = ap.parse_args(argv)
    if a.cmd == "count":
        qc = count(a.sample, a.fragments, a.genome, a.features, a.outdir, window=a.window,
                   frag_len=a.fragment_length, min_cpgs=a.min_cpgs, seed=a.seed)
        print("\t".join(f"{k}={qc[k]}" for k in QC_COLUMNS), file=sys.stderr)
    else:
        qc_table(a.reports, a.samplesheet, a.out, min_mapping=a.min_mapping,
                 min_fragments=a.min_fragments, min_relh=a.min_relh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
