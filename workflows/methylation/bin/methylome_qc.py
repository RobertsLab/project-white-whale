#!/usr/bin/env python3
"""Per-sample CpG summaries and the QC gate (PLAN.md §5). Standard library only:
it runs in the pinned python container, and tests/test_methylome_qc.py loads
it from this file.

    methylome_qc.py cpg-summary SAMPLE COV_GZ OUT_TSV
    methylome_qc.py qc-table --samplesheet S --out qc_summary.tsv \
        [--min-mapping 20] [--min-conversion 98] [--min-cpgs-5x 1000000] FILE...

`cpg-summary` reads one strand-merged CpG file from coverage2cytosine
(chrom, start, end, % methylated, count methylated, count unmethylated) and
writes one row: CpGs covered at 1x, 5x and 10x, their mean depth, the
methylation of all calls pooled, and the mean per-CpG methylation at >= 5x.

`qc-table` joins, per sample in the sample sheet, the Bismark alignment
report (<S>_PE_report.txt or <S>_SE_report.txt), the deduplication report
(<S>.dedup.deduplication_report.txt; none for RRBS), the extractor's
splitting report and the CpG summary, and applies the gate:

* mapping efficiency >= --min-mapping (%)
* bisulfite conversion >= --min-conversion (%), estimated as 100 minus the
  methylation of non-CpG cytosines, which are almost never methylated in
  oyster
* CpGs covered at >= 5x >= --min-cpgs-5x

A library that fails is kept in the table with qc_pass = false and the
reasons, and left out of every summary built on it.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

CPG_COLUMNS = ["sample_id", "cpgs_1x", "cpgs_5x", "cpgs_10x", "mean_depth",
               "cpg_meth_pooled_pct", "cpg_meth_mean_5x_pct"]
QC_COLUMNS = ["sample_id", "condition", "pairs_or_reads", "mapping_pct", "duplication_pct",
              "cpg_meth_pct", "non_cpg_meth_pct", "conversion_pct",
              "cpgs_1x", "cpgs_5x", "cpgs_10x", "mean_depth", "cpg_meth_mean_5x_pct",
              "qc_pass", "qc_reasons"]


def _number(text: str, label: str) -> float | None:
    m = re.search(rf"^{re.escape(label)}:?\s+([\d.]+)", text, re.MULTILINE)
    return float(m.group(1)) if m else None


def _pct(num: float | None, den: float | None) -> float | None:
    return round(100.0 * num / den, 3) if num is not None and den else None


def parse_alignment_report(text: str) -> dict:
    """Bismark's <S>_PE_report.txt / <S>_SE_report.txt."""
    total = _number(text, "Sequence pairs analysed in total") or _number(text, "Sequences analysed in total")
    return {"pairs_or_reads": int(total) if total is not None else None,
            "mapping_pct": _number(text, "Mapping efficiency")}


def parse_dedup_report(text: str) -> dict:
    m = re.search(r"Total number duplicated alignments removed:\s+\d+\s+\(([\d.]+)%\)", text)
    return {"duplication_pct": float(m.group(1)) if m else None}


def parse_splitting_report(text: str) -> dict:
    """Methylation from the extractor's counts (after deduplication)."""
    def ctx(name):
        meth = _number(text, f"Total methylated C's in {name} context")
        conv = _number(text, f"Total C to T conversions in {name} context")
        return meth, conv
    cpg_m, cpg_u = ctx("CpG")
    chg_m, chg_u = ctx("CHG")
    chh_m, chh_u = ctx("CHH")
    non_m = None if chg_m is None or chh_m is None else chg_m + chh_m
    non_u = None if chg_u is None or chh_u is None else chg_u + chh_u
    cpg = _pct(cpg_m, None if cpg_m is None or cpg_u is None else cpg_m + cpg_u)
    non = _pct(non_m, None if non_m is None or non_u is None else non_m + non_u)
    return {"cpg_meth_pct": cpg, "non_cpg_meth_pct": non,
            "conversion_pct": round(100.0 - non, 3) if non is not None else None}


def cpg_summary(sample_id: str, cov_gz: Path) -> dict:
    n1 = n5 = n10 = 0
    depth = meth = 0
    mean5 = 0.0
    with gzip.open(cov_gz, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 6:
                continue
            m, u = int(f[4]), int(f[5])
            d = m + u
            if d == 0:
                continue
            n1 += 1
            depth += d
            meth += m
            if d >= 5:
                n5 += 1
                mean5 += m / d
                if d >= 10:
                    n10 += 1
    return {"sample_id": sample_id, "cpgs_1x": n1, "cpgs_5x": n5, "cpgs_10x": n10,
            "mean_depth": round(depth / n1, 3) if n1 else 0.0,
            "cpg_meth_pooled_pct": _pct(meth, depth),
            "cpg_meth_mean_5x_pct": round(100.0 * mean5 / n5, 3) if n5 else None}


def _find(files: dict[str, Path], *names: str) -> Path | None:
    for n in names:
        if n in files:
            return files[n]
    return None


def gate(row: dict, *, min_mapping: float, min_conversion: float, min_cpgs_5x: int) -> list[str]:
    reasons = []
    if row.get("mapping_pct") is None or row["mapping_pct"] < min_mapping:
        reasons.append(f"mapping {row.get('mapping_pct')}% < {min_mapping}%")
    if row.get("conversion_pct") is None or row["conversion_pct"] < min_conversion:
        reasons.append(f"conversion {row.get('conversion_pct')}% < {min_conversion}%")
    if (row.get("cpgs_5x") or 0) < min_cpgs_5x:
        reasons.append(f"{row.get('cpgs_5x') or 0} CpGs at >=5x < {min_cpgs_5x}")
    return reasons


def qc_table(samplesheet: Path, paths: list[Path], *, min_mapping: float, min_conversion: float,
             min_cpgs_5x: int) -> list[dict]:
    files = {p.name: p for p in paths}
    with open(samplesheet, newline="") as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    rows = []
    for s in samples:
        sid = s["sample_id"]
        row = {"sample_id": sid, "condition": s.get("condition", "")}
        missing = []
        aln = _find(files, f"{sid}_PE_report.txt", f"{sid}_SE_report.txt")
        if aln:
            row.update(parse_alignment_report(aln.read_text()))
        else:
            missing.append("alignment report")
        dedup = _find(files, f"{sid}.dedup.deduplication_report.txt")
        row.update(parse_dedup_report(dedup.read_text()) if dedup else {"duplication_pct": None})
        split = _find(files, f"{sid}.dedup_splitting_report.txt", f"{sid}_splitting_report.txt")
        if split:
            row.update(parse_splitting_report(split.read_text()))
        else:
            missing.append("splitting report")
        cpg = _find(files, f"{sid}.cpg_summary.tsv")
        if cpg:
            with open(cpg, newline="") as fh:
                c = next(csv.DictReader(fh, delimiter="\t"))
            row.update({k: c[k] for k in CPG_COLUMNS[1:] if k != "cpg_meth_pooled_pct"})
            for k in ("cpgs_1x", "cpgs_5x", "cpgs_10x"):
                row[k] = int(row[k])
        else:
            missing.append("CpG summary")
        reasons = [f"no {m}" for m in missing] + gate(row, min_mapping=min_mapping,
                                                         min_conversion=min_conversion, min_cpgs_5x=min_cpgs_5x)
        row["qc_pass"] = "false" if reasons else "true"
        row["qc_reasons"] = "; ".join(reasons)
        rows.append(row)
    return rows


def write_tsv(rows: list[dict], columns: list[str], out: Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r[k]) for k in columns})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cpg-summary")
    c.add_argument("sample_id")
    c.add_argument("cov_gz", type=Path)
    c.add_argument("out", type=Path)
    q = sub.add_parser("qc-table")
    q.add_argument("--samplesheet", type=Path, required=True)
    q.add_argument("--out", type=Path, required=True)
    q.add_argument("--min-mapping", type=float, default=20.0)
    q.add_argument("--min-conversion", type=float, default=98.0)
    q.add_argument("--min-cpgs-5x", type=int, default=1_000_000)
    q.add_argument("files", nargs="*", type=Path)
    a = ap.parse_args(argv)
    if a.cmd == "cpg-summary":
        write_tsv([cpg_summary(a.sample_id, a.cov_gz)], CPG_COLUMNS, a.out)
    else:
        rows = qc_table(a.samplesheet, a.files, min_mapping=a.min_mapping,
                        min_conversion=a.min_conversion, min_cpgs_5x=a.min_cpgs_5x)
        write_tsv(rows, QC_COLUMNS, a.out)
        failed = [r["sample_id"] for r in rows if r["qc_pass"] != "true"]
        print(f"{len(rows) - len(failed)} of {len(rows)} libraries pass the QC gate"
              + (f"; failed: {', '.join(failed)}" if failed else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
