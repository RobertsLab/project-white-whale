"""Plant C/T SNPs for a real run of workflows/snp_mask on the fixtures.

    python3 tests/plant_snps.py results/test/snp_reference      # before the mask run
    python3 tests/plant_snps.py --check results/test/snp_reference results/test/snp_mask

The fixture reads come from tests/fixtures/reference/TEST/genomic.fna and have
no SNPs. This writes a copy of that genome in which N_SITES `TG` dinucleotides
read `CG`, plus truth.tsv. Against that reference the reads' T at each planted
site is a C/T SNP at a CpG, which bisulfite alone could explain: exactly the
case the mask exists for. --check then requires the mask to hold nearly all
planted CpGs and almost nothing else.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENOME = HERE / "fixtures" / "reference" / "TEST" / "genomic.fna"
N_SITES = 20
MIN_FOUND = 0.8      # of planted CpGs masked (some sites have too few bottom-strand reads)
MAX_EXTRA = 5        # masked CpGs that were not planted (sequencing errors)


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def plant(out: Path) -> None:
    seqs = read_fasta(GENOME)
    sites = []
    for name, seq in seqs.items():
        # Evenly spaced TG sites, away from the ends where coverage is thin.
        cands = [i for i in range(500, len(seq) - 500) if seq[i:i + 2] == "TG"]
        step = max(1, len(cands) // (N_SITES // len(seqs)))
        for i in cands[::step][: N_SITES // len(seqs)]:
            sites.append((name, i))
    for name, i in sites:
        s = seqs[name]
        seqs[name] = s[:i] + "C" + s[i + 1:]
    d = out / "TEST"
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "genomic.fna", "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name}\n")
            fh.writelines(s[j:j + 80] + "\n" for j in range(0, len(s), 80))
    (out / "truth.tsv").write_text("".join(f"{n}\t{i}\n" for n, i in sites))
    print(f"planted {len(sites)} C/T SNPs at CpGs in {d / 'genomic.fna'}")


def check(ref: Path, mask_dir: Path) -> int:
    truth = {tuple([n, int(i)]) for n, i in (ln.split("\t") for ln in (ref / "truth.tsv").read_text().split("\n") if ln)}
    with gzip.open(mask_dir / "cpg_mask.bed.gz", "rt") as fh:
        masked = {(f[0], int(f[1])) for f in (ln.split("\t") for ln in fh)}
    found, extra = len(truth & masked), len(masked - truth)
    print(f"masked {found} of {len(truth)} planted CpGs; {extra} other CpG(s) masked")
    ok = found >= MIN_FOUND * len(truth) and extra <= MAX_EXTRA
    print("ok: the mask finds the planted SNPs" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    if sys.argv[1:2] == ["--check"]:
        sys.exit(check(Path(sys.argv[2]), Path(sys.argv[3])))
    plant(Path(sys.argv[1]))
