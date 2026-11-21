"""Check a real (not stub) RNA-seq workflow run on the fixtures against the reads it was given.

    python3 tests/check_rnaseq_run.py results/test/TEST2026_FIXTURE_WGBS/rnaseq

The fixture reads (tests/fixtures/make_rnaseq_fixtures.py) are a stranded
(dUTP) library per sample with known fragments per transcript, plus
fragments that belong to no transcript: some inside introns or intergenic DNA,
and some unspliced across an exon-intron boundary, which the genome decoys
should take. A run that trims, quantifies and aggregates them correctly gives
every gene about its simulated count, both G1 isoforms the right total, no
expression to the gene with no reads, and calls the libraries reverse-stranded.
Exits non-zero, listing every failure.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

TRUTH = Path(__file__).resolve().parent / "fixtures" / "rnaseq" / "truth.tsv"
GENES = ["G1", "G2", "G3", "G4"]
# A gene's estimate may exceed its simulated fragments by the unspliced fragments
# across its own exon-intron boundaries that the decoys did not take: with
# 35-75 exonic bases, many still pass Salmon's alignment score (the decoys take
# about a third of them paired-end and a sixth single-end on these reads).
SLACK = 0.05                  # below simulated, or above simulated + pre-mRNA, as a share of simulated
MAX_ABSENT = 5.0              # estimated fragments for a gene simulated with none
MIN_DECOY_SHARE = 0.10        # of the pre-mRNA fragments, discarded as best-mapped to decoys


def read_matrix(path: Path) -> dict[str, dict[str, float]]:
    with gzip.open(path, "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    samples = rows[0][1:]
    return {s: {r[0]: float(r[i + 1]) for r in rows[1:]} for i, s in enumerate(samples)}


def check(run: Path) -> list[str]:
    problems = []
    truth: dict[str, dict[str, int]] = {}
    with open(TRUTH, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            truth.setdefault(r["sample_id"], {})[r["feature"]] = int(r["fragments"])
    for name in ("run_manifest.json", "qc_summary.tsv", "multiqc_report.html",
                 "matrices/gene_tpm.tsv.gz", "matrices/gene_counts.tsv.gz", "matrices/gene_length.tsv.gz",
                 "matrices/transcript_tpm.tsv.gz", "matrices/transcript_counts.tsv.gz"):
        if not (run / name).exists():
            problems.append(f"{run / name} is missing")
    if problems:
        return problems

    counts = read_matrix(run / "matrices" / "gene_counts.tsv.gz")
    tpm = read_matrix(run / "matrices" / "gene_tpm.tsv.gz")
    tx_counts = read_matrix(run / "matrices" / "transcript_counts.tsv.gz")
    with open(run / "qc_summary.tsv", newline="") as fh:
        qc = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    for sample, want in truth.items():
        if sample not in counts:
            problems.append(f"{sample}: not in gene_counts.tsv.gz")
            continue
        for gene in GENES:
            got, exp = counts[sample].get(gene), want[gene]
            if got is None:
                problems.append(f"{sample} {gene}: not in gene_counts.tsv.gz")
            elif exp >= 100 and not exp * (1 - SLACK) <= got <= (exp + want.get(f"premrna_{gene}", 0)) * (1 + SLACK):
                problems.append(f"{sample} {gene}: {got:.0f} fragments estimated, {exp} simulated "
                                f"(+{want.get(f'premrna_{gene}', 0)} unspliced across its boundaries)")
            elif exp == 0 and got > MAX_ABSENT:
                problems.append(f"{sample} {gene}: {got:.0f} fragments estimated for a gene with none")
        g1 = sum(v for t, v in tx_counts[sample].items() if t in ("XM_900001.1", "XM_900002.1"))
        if abs(g1 - counts[sample]["G1"]) > 0.01:
            problems.append(f"{sample}: G1's count {counts[sample]['G1']} is not the sum of its isoforms ({g1})")
        if abs(sum(tpm[sample].values()) - 1e6) > 10:
            problems.append(f"{sample}: gene TPM sums to {sum(tpm[sample].values()):.1f}, not 1e6")

        meta = json.loads((run / "salmon" / sample / "aux_info" / "meta_info.json").read_text())
        decoys = meta.get("num_decoy_fragments", 0)
        if decoys < MIN_DECOY_SHARE * want["premrna"]:
            problems.append(f"{sample}: {decoys} fragments went to decoys, of {want['premrna']} simulated pre-mRNA")
        total = sum(want[k] for k in (*GENES, "genomic", "premrna"))
        tx_share = 100.0 * sum(want[g] for g in GENES) / total
        max_share = 100.0 * (sum(want[g] for g in GENES) + want["premrna"] - decoys) / total
        row = qc.get(sample, {})
        mapping = float(row.get("mapping_pct") or 0)
        if not tx_share - 2 <= mapping <= max_share + 0.5:
            problems.append(f"{sample}: mapping {mapping}%, expected {tx_share:.1f}-{max_share:.1f}% "
                            "(transcript fragments, plus pre-mRNA the decoys did not take)")
        if row.get("strandedness") != "reverse":
            problems.append(f"{sample}: strandedness {row.get('strandedness')!r}, simulated reverse (dUTP)")
        if row.get("qc_pass") != "true":
            problems.append(f"{sample}: failed the QC gate ({row.get('qc_reasons')})")
    return problems


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    problems = check(Path(sys.argv[1]))
    for p in problems:
        print(f"FAIL  {p}")
    if not problems:
        print("ok: gene and isoform counts, decoys, strandedness and the QC gate match the simulated reads")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
