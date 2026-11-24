"""Check a real (not stub) workflow run on the fixtures against the methylation they were made with.

    python3 tests/check_real_run.py results/test/TEST2026_FIXTURE_WGBS

The fixture reads (tests/fixtures/make_fixtures.py) are simulated with every
CpG methylated with probability 0.8 (83% of CpGs in the generated genome) and
every other C unmethylated. A run that aligns, deduplicates and calls them
correctly maps nearly every read, reports CpG methylation close to that, and
almost no methylation outside CpGs. Exits non-zero, listing every failure.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from pathlib import Path

SAMPLES = {"S1": "PE", "S2": "SE", "S3": "PE"}   # S3: S1 deposited with mates swapped (swap_mates)
MIN_MAPPING = 95.0          # % of reads with a unique best hit
CPG_RANGE = (75.0, 90.0)    # % methylated at CpGs; the truth is 83.4%
MAX_NON_CPG = 1.0           # % methylated in CHG and CHH contexts
MIN_CPG_SITES = 900         # of 980 CpGs, covered by at least one read
GENE_RANGE = (70.0, 95.0)   # % methylated in each fixture gene body (features/genes.bed)


def report_value(text: str, label: str) -> float:
    m = re.search(rf"^{re.escape(label)}:\s*([\d.]+)%", text, re.MULTILINE)
    if not m:
        raise ValueError(f"no '{label}' line")
    return float(m.group(1))


def check(run: Path) -> list[str]:
    problems = []
    for sample, layout in SAMPLES.items():
        report = run / "qc" / f"{sample}_{layout}_report.txt"
        if not report.exists():
            problems.append(f"{sample}: {report} is missing")
            continue
        text = report.read_text()
        try:
            mapping = report_value(text, "Mapping efficiency")
            cpg = report_value(text, "C methylated in CpG context")
            non_cpg = {c: report_value(text, f"C methylated in {c} context") for c in ("CHG", "CHH")}
        except ValueError as exc:
            problems.append(f"{sample}: {report.name} has {exc}; not a Bismark report from a real run")
            continue
        if mapping < MIN_MAPPING:
            problems.append(f"{sample}: mapping efficiency {mapping}% < {MIN_MAPPING}%")
        if not CPG_RANGE[0] <= cpg <= CPG_RANGE[1]:
            problems.append(f"{sample}: CpG methylation {cpg}% outside {CPG_RANGE}")
        for context, level in non_cpg.items():
            if level > MAX_NON_CPG:
                problems.append(f"{sample}: {context} methylation {level}% > {MAX_NON_CPG}% "
                                "(bisulfite conversion not modelled correctly?)")

        cov = run / "cpg" / f"{sample}.CpG_report.merged_CpG_evidence.cov.gz"
        if not cov.exists():
            problems.append(f"{sample}: {cov} is missing")
            continue
        with gzip.open(cov, "rt") as fh:
            rows = [line.split("\t") for line in fh if line.strip()]
        if len(rows) < MIN_CPG_SITES:
            problems.append(f"{sample}: {len(rows)} CpG sites in {cov.name}, expected at least {MIN_CPG_SITES}")
        meth = sum(int(r[4]) for r in rows)
        total = meth + sum(int(r[5]) for r in rows)
        if total and not CPG_RANGE[0] <= 100 * meth / total <= CPG_RANGE[1]:
            problems.append(f"{sample}: {cov.name} gives {100 * meth / total:.1f}% CpG methylation")
        print(f"{sample}: mapping {mapping}%, CpG {cpg}%, {len(rows)} CpG sites")

    qc = run / "qc_summary.tsv"
    if not qc.exists():
        problems.append(f"{qc} is missing")
    else:
        with open(qc, newline="") as fh:
            gate = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for sample in SAMPLES:
            row = gate.get(sample)
            if row is None:
                problems.append(f"{sample}: not in qc_summary.tsv")
            elif row["qc_pass"] != "true":
                problems.append(f"{sample}: fails the QC gate ({row['qc_reasons']})")

    for sample in SAMPLES:
        table = run / "features" / f"{sample}.features.tsv.gz"
        if not table.exists():
            problems.append(f"{sample}: {table} is missing")
            continue
        with gzip.open(table, "rt") as fh:
            bodies = [r for r in csv.DictReader(fh, delimiter="\t") if r["feature"] == "gene_body"]
        if [r["gene_id"] for r in bodies] != ["G1", "G2", "G3", "G4"]:
            problems.append(f"{sample}: {table.name} lists gene bodies {[r['gene_id'] for r in bodies]}")
        for r in bodies:
            # Each gene holds 20-180 CpGs, so allow more spread than the whole-genome range.
            if not r["meth_pct"] or not GENE_RANGE[0] <= float(r["meth_pct"]) <= GENE_RANGE[1]:
                problems.append(f"{sample}: gene {r['gene_id']} body methylation {r['meth_pct'] or 'missing'}%")

    manifest = run / "run_manifest.json"
    if not manifest.exists():
        problems.append(f"{manifest} is missing")
    else:
        versions = json.loads(manifest.read_text())["tool_versions"]
        stubbed = sorted(p for p, tools in versions.items() if "stub" in tools.values())
        if stubbed:
            problems.append(f"run_manifest.json records stub versions for {stubbed}; was this run with -stub?")
    return problems


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    problems = check(Path(sys.argv[1]))
    for p in problems:
        print(f"FAIL {p}")
    print("ok: the real run recovers the fixtures' methylation" if not problems else
          f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
