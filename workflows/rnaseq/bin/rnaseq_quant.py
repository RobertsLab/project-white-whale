#!/usr/bin/env python3
"""Gene-level matrices and the QC table for the RNA-seq workflow. Standard
library only: it runs in the pinned python container, and
tests/test_rnaseq_quant.py loads it from this file.

    rnaseq_quant.py matrices --samplesheet S --tx2gene T --out-dir D QUANT_DIR...
    rnaseq_quant.py qc-table --samplesheet S --matrices D --out qc_summary.tsv \
        [--min-mapping 30] [--min-mapped 1000000] [--min-genes 1000] FILE_OR_DIR...

`matrices` reads each sample's Salmon output folder (named after the sample,
holding quant.sf) and writes, in sample-sheet order:

* transcript_tpm.tsv.gz, transcript_counts.tsv.gz   transcript x sample, as Salmon
  estimated them (TPM, NumReads);
* gene_tpm.tsv.gz, gene_counts.tsv.gz, gene_length.tsv.gz   gene x sample,
  aggregated with tx2gene the way tximport does with countsFromAbundance =
  "no": a gene's TPM and estimated count are the sums over its transcripts,
  and its length is the TPM-weighted mean of their effective lengths (the
  plain mean where the gene has no expression in that sample). The length
  matrix is what DESeq2's tximport import uses as an offset;
* matrices_report.json   how many transcripts matched tx2gene, and how many
  genes of tx2gene had no transcript in the index.

Gene IDs are tx2gene's (NCBI GeneIDs in a real bundle), the same IDs as the
reference bundle's features/genes.bed and the methylation gene matrices.

`qc-table` joins, per sample, fastp's JSON (<S>.fastp.json), Salmon's
aux_info/meta_info.json and lib_format_counts.json (in the folder <S>/) and
the gene TPM matrix, and applies the gate:

* Salmon mapping rate >= --min-mapping (%), decoy hits not counted as mapped;
* fragments mapped >= --min-mapped;
* genes with TPM >= 1 >= --min-genes.

A library that fails is kept with qc_pass = false and the reasons. The
table also reports each library's observed strandedness (from Salmon's
per-format fragment counts), so a study that mixes library preparations shows it.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import sys
from pathlib import Path

QC_COLUMNS = ["sample_id", "condition", "layout", "reads_in", "reads_after_trim_pct", "q30_pct",
              "fragments_processed", "fragments_mapped", "mapping_pct", "decoy_pct",
              "library_type", "strandedness", "sense_pct", "compatible_fragment_pct", "genes_tpm_ge1",
              "qc_pass", "qc_reasons"]
DETECTED_TPM = 1.0
STRANDED_PCT = 80.0   # sense (or antisense) share of strand-informative fragments to call a library stranded


class QuantError(RuntimeError):
    pass


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt(x: float) -> str:
    s = f"{x:.4f}".rstrip("0").rstrip(".")
    return "0" if s in ("", "-0") else s


def gene_key(g: str):
    return (0, int(g), "") if g.isdigit() else (1, 0, g)


def read_tx2gene(path: Path) -> dict[str, str]:
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or not f[0] or f[0] in ("transcript_id", "tx", "TXNAME"):
                continue
            if f[0] in out and out[f[0]] != f[1]:
                raise QuantError(f"{path}: {f[0]} maps to both {out[f[0]]} and {f[1]}")
            out[f[0]] = f[1]
    if not out:
        raise QuantError(f"{path} has no transcript -> gene pairs")
    return out


def read_quant(path: Path) -> list[dict]:
    """quant.sf rows: Name, Length, EffectiveLength, TPM, NumReads (numbers as floats)."""
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({"Name": r["Name"], "Length": float(r["Length"]),
                         "EffectiveLength": float(r["EffectiveLength"]),
                         "TPM": float(r["TPM"]), "NumReads": float(r["NumReads"])})
    return rows


def aggregate(quant: list[dict], tx2gene: dict[str, str]) -> dict[str, dict]:
    """gene -> {tpm, counts, length} for one sample (tximport, countsFromAbundance = "no")."""
    acc: dict[str, dict] = {}
    for r in quant:
        g = tx2gene.get(r["Name"])
        if g is None:
            continue
        a = acc.setdefault(g, {"tpm": 0.0, "counts": 0.0, "wlen": 0.0, "lens": []})
        a["tpm"] += r["TPM"]
        a["counts"] += r["NumReads"]
        a["wlen"] += r["TPM"] * r["EffectiveLength"]
        a["lens"].append(r["EffectiveLength"])
    out = {}
    for g, a in acc.items():
        length = a["wlen"] / a["tpm"] if a["tpm"] > 0 else statistics.fmean(a["lens"])
        out[g] = {"tpm": a["tpm"], "counts": a["counts"], "length": length}
    return out


def _write_gz(path: Path, rows) -> None:
    # mtime=0: the same inputs give byte-identical files, so a re-run commits nothing.
    with open(path, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
        for row in rows:
            gz.write(("\t".join(row) + "\n").encode())


def build_matrices(samples: list[str], quant_dirs: dict[str, Path], tx2gene: dict[str, str],
                   out_dir: Path) -> dict:
    missing = [s for s in samples if s not in quant_dirs]
    if missing:
        raise QuantError(f"no Salmon output for sample(s) {missing}")
    quants = {s: read_quant(quant_dirs[s] / "quant.sf") for s in samples}
    names = [r["Name"] for r in quants[samples[0]]]
    for s in samples[1:]:
        if [r["Name"] for r in quants[s]] != names:
            raise QuantError(f"{s}: its quant.sf lists other transcripts than {samples[0]}'s "
                             "(quantified against another index?)")
    matched = [t for t in names if t in tx2gene]
    if not matched:
        raise QuantError("no transcript in quant.sf is in tx2gene; are rna.fna and tx2gene.tsv from the same bundle?")
    unmatched = [t for t in names if t not in tx2gene]

    out_dir.mkdir(parents=True, exist_ok=True)
    header = ["transcript_id", *samples]
    for metric, col in (("tpm", "TPM"), ("counts", "NumReads")):
        _write_gz(out_dir / f"transcript_{metric}.tsv.gz",
                  [header, *([t, *(fmt(quants[s][i][col]) for s in samples)] for i, t in enumerate(names))])

    per_sample = {s: aggregate(quants[s], tx2gene) for s in samples}
    genes = sorted({tx2gene[t] for t in matched}, key=gene_key)
    for metric in ("tpm", "counts", "length"):
        _write_gz(out_dir / f"gene_{metric}.tsv.gz",
                  [["gene_id", *samples], *([g, *(fmt(per_sample[s][g][metric]) for s in samples)] for g in genes)])

    report = {
        "samples": samples,
        "transcripts_quantified": len(names),
        "transcripts_in_tx2gene": len(matched),
        "transcripts_not_in_tx2gene": len(unmatched),
        "transcripts_not_in_tx2gene_examples": unmatched[:10],
        "genes": len(genes),
        "tx2gene_genes_without_quantified_transcript": len(set(tx2gene.values()) - set(genes)),
        "aggregation": "tximport countsFromAbundance='no': TPM and counts summed per gene; "
                       "length = TPM-weighted mean effective length (plain mean when TPM is 0)",
    }
    (out_dir / "matrices_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def _pct(num, den) -> float | None:
    return round(100.0 * num / den, 3) if num is not None and den else None


def parse_fastp(data: dict) -> dict:
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    paired = "read2_before_filtering" in data
    total = before.get("total_reads")
    return {"layout": "PE" if paired else "SE",
            "reads_in": (total // 2 if paired else total) if total is not None else None,
            "reads_after_trim_pct": _pct(after.get("total_reads"), total),
            "q30_pct": round(100.0 * after["q30_rate"], 3) if after.get("q30_rate") is not None else None}


def strandedness(lib: dict) -> tuple[str | None, float | None]:
    """Observed strandedness from lib_format_counts.json's per-format counts.

    Salmon's `expected_format` is what -l A settled on, which falls back to an
    unstranded type on small libraries; the counts say what the reads did.
    Sense = read 1 on the transcript's strand (ISF, SF). Returns ("reverse" /
    "forward" / "unstranded", sense %) or (None, None) without counts.
    """
    sense = sum(lib.get(k, 0) for k in ("ISF", "SF", "MSF", "OSF"))
    anti = sum(lib.get(k, 0) for k in ("ISR", "SR", "MSR", "OSR"))
    if not sense + anti:
        return None, None
    pct = round(100.0 * sense / (sense + anti), 3)
    call = "forward" if pct >= STRANDED_PCT else "reverse" if pct <= 100 - STRANDED_PCT else "unstranded"
    return call, pct


def parse_salmon(meta: dict, lib: dict | None) -> dict:
    processed, mapped = meta.get("num_processed"), meta.get("num_mapped")
    types = meta.get("library_types") or []
    lib = lib or {}
    ratio = lib.get("compatible_fragment_ratio")
    stranded, sense_pct = strandedness(lib)
    return {"strandedness": stranded, "sense_pct": sense_pct,
            "fragments_processed": processed, "fragments_mapped": mapped,
            "mapping_pct": round(float(meta["percent_mapped"]), 3) if meta.get("percent_mapped") is not None
            else _pct(mapped, processed),
            "decoy_pct": _pct(meta.get("num_decoy_fragments"), processed),
            "library_type": lib.get("expected_format") or (types[0] if types else None),
            "compatible_fragment_pct": round(100.0 * ratio, 3) if ratio is not None else None}


def genes_detected(matrix_dir: Path) -> dict[str, int]:
    path = matrix_dir / "gene_tpm.tsv.gz"
    with gzip.open(path, "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    samples = rows[0][1:]
    return {s: sum(float(r[i + 1]) >= DETECTED_TPM for r in rows[1:]) for i, s in enumerate(samples)}


def qc_table(sheet: list[dict], files: list[Path], matrix_dir: Path, *, min_mapping: float,
             min_mapped: int, min_genes: int) -> list[dict]:
    fastp = {p.name[: -len(".fastp.json")]: p for p in files if p.name.endswith(".fastp.json")}
    salmon = {p.name: p for p in files if p.is_dir()}
    detected = genes_detected(matrix_dir)
    rows = []
    for s in sheet:
        sid = s["sample_id"]
        row = {"sample_id": sid, "condition": s.get("condition", "")}
        if sid in fastp:
            row.update(parse_fastp(json.loads(fastp[sid].read_text())))
        if sid in salmon:
            meta = json.loads((salmon[sid] / "aux_info" / "meta_info.json").read_text())
            libp = salmon[sid] / "lib_format_counts.json"
            row.update(parse_salmon(meta, json.loads(libp.read_text()) if libp.exists() else None))
        row["genes_tpm_ge1"] = detected.get(sid)
        reasons = []
        if row.get("mapping_pct") is None or row["mapping_pct"] < min_mapping:
            reasons.append(f"mapping {row.get('mapping_pct')}% < {min_mapping}%")
        if row.get("fragments_mapped") is None or row["fragments_mapped"] < min_mapped:
            reasons.append(f"mapped fragments {row.get('fragments_mapped')} < {min_mapped}")
        if row["genes_tpm_ge1"] is None or row["genes_tpm_ge1"] < min_genes:
            reasons.append(f"genes with TPM >= 1 {row['genes_tpm_ge1']} < {min_genes}")
        row["qc_pass"] = "false" if reasons else "true"
        row["qc_reasons"] = "; ".join(reasons)
        rows.append(row)
    return rows


def write_table(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=QC_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: "" if r.get(k) is None else r[k] for k in QC_COLUMNS})


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("matrices")
    m.add_argument("--samplesheet", required=True, type=Path)
    m.add_argument("--tx2gene", required=True, type=Path)
    m.add_argument("--out-dir", required=True, type=Path)
    m.add_argument("quant_dirs", nargs="+", type=Path)
    q = sub.add_parser("qc-table")
    q.add_argument("--samplesheet", required=True, type=Path)
    q.add_argument("--matrices", required=True, type=Path)
    q.add_argument("--out", required=True, type=Path)
    q.add_argument("--min-mapping", type=float, default=30.0)
    q.add_argument("--min-mapped", type=int, default=1_000_000)
    q.add_argument("--min-genes", type=int, default=1000)
    q.add_argument("files", nargs="+", type=Path)
    a = ap.parse_args(argv)
    sheet = read_tsv(a.samplesheet)
    try:
        if a.cmd == "matrices":
            build_matrices([r["sample_id"] for r in sheet], {p.name: p for p in a.quant_dirs},
                           read_tx2gene(a.tx2gene), a.out_dir)
        else:
            write_table(qc_table(sheet, a.files, a.matrices, min_mapping=a.min_mapping,
                                 min_mapped=a.min_mapped, min_genes=a.min_genes), a.out)
    except QuantError as exc:
        print(f"rnaseq_quant.py: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
