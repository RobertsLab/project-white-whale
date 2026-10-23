"""Fetch, verify and derive the reference bundle every sample is processed against.

A bundle lives in ``$WHALE_REFERENCE/<ACCESSION>/`` (default
``data/reference/<ACCESSION>/``, not committed) and holds:

* ``genomic.fna``   genome FASTA (Bismark), from NCBI
* ``genomic.gtf``   RefSeq annotation, from NCBI
* ``rna.fna``       transcript FASTA, only with ``--rnaseq`` (Salmon, Goal 2)
* ``tx2gene.tsv``   transcript -> NCBI GeneID, only with ``--rnaseq``
* ``chrom.sizes``   sequence name and length, from the genome
* ``features/``     BED6 tracks the per-sample summaries and analyses use
                    (PLAN.md §6.1): genes, exons, introns, promoters, intergenic

plus the compressed NCBI originals, NCBI's ``md5checksums.txt`` and
``reference_provenance.json``: where every file came from, its checksums, the
rules the features were derived with, and whether the GTF's annotation
release is the one registered in data/reference/genome_assemblies.yaml (NCBI
can re-annotate an assembly in place).

Every download is checked against NCBI's published MD5 before use. Bismark
indexes are not built here; the workflow builds them. The CpG index and the
C/T SNP mask are later Phase 1 work.

Adapted from AREE's src/intake/reference_bundle.py (github.com/sr320/AREE).
"""
from __future__ import annotations

import gzip
import json
import os
import re
import shutil
from datetime import date
from pathlib import Path

import yaml

from whale.checksums import md5_file, sha256_file
from whale.net import get_text

NCBI_GENOMES_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
TIMEOUT = 300
PROVENANCE_NAME = "reference_provenance.json"
BUNDLER_VERSION = "1.0.0"

CORE_FILES = {"_genomic.fna.gz": "genomic.fna", "_genomic.gtf.gz": "genomic.gtf"}
RNASEQ_FILES = {"_rna.fna.gz": "rna.fna"}

# Promoter window around each gene's TSS, strand-aware (PLAN.md §6.1).
PROMOTER_UPSTREAM = 1000
PROMOTER_DOWNSTREAM = 200
FEATURE_RULES = {
    "genes": "GTF `gene` features; name = NCBI GeneID; score column = gene_biotype",
    "exons": "GTF `exon` features of every transcript, de-duplicated by (chrom, start, end, strand, gene)",
    "introns": "gaps between consecutive exons of each transcript, de-duplicated like exons",
    "promoters": f"per gene, {PROMOTER_UPSTREAM} bp upstream to {PROMOTER_DOWNSTREAM} bp downstream of the gene's TSS, strand-aware, clipped to the sequence",
    "intergenic": "sequence not covered by any gene, either strand",
}


class ReferenceError(RuntimeError):
    """A bundle could not be fetched, or does not match what was recorded."""


# --- assemblies ----------------------------------------------------------------

def load_assemblies(repo: Path) -> list[dict]:
    doc = yaml.safe_load((Path(repo) / "data" / "reference" / "genome_assemblies.yaml").read_text())
    return doc["assemblies"]


def resolve_assembly(repo: Path, value: str | None = None) -> dict:
    """The registered assembly named by id or accession; the default one when value is None."""
    assemblies = load_assemblies(repo)
    if value is None:
        defaults = [a for a in assemblies if a.get("default")]
        if len(defaults) != 1:
            raise ReferenceError("genome_assemblies.yaml must mark exactly one assembly `default: true`.")
        record = defaults[0]
    else:
        matches = [a for a in assemblies
                   if value in (a["assembly_id"], a["ncbi_assembly_accession"], a.get("genbank_accession"))]
        if not matches:
            raise ReferenceError(f"{value!r} is not in data/reference/genome_assemblies.yaml; register it first.")
        record = matches[0]
    if not record["ncbi_assembly_accession"].startswith("GCF_") or not record.get("annotation_release"):
        raise ReferenceError(f"{record['assembly_id']} has no NCBI RefSeq annotation, so NCBI publishes no GTF for it.")
    return record


def default_bundle_dir(repo: Path, accession: str) -> Path:
    base = os.environ.get("WHALE_REFERENCE")
    return (Path(base) if base else Path(repo) / "data" / "reference") / accession


# --- NCBI genomes FTP ----------------------------------------------------------

def accession_dir_url(accession: str) -> str:
    """Parent directory holding one assembly's folder, e.g. .../GCF/963/853/765/."""
    m = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+", accession)
    if not m:
        raise ReferenceError(f"not an NCBI assembly accession: {accession!r}")
    prefix, a, b, c = m.groups()
    return f"{NCBI_GENOMES_BASE}/{prefix}/{a}/{b}/{c}/"


def find_assembly_folder(listing_html: str, accession: str) -> str:
    """The folder ``<accession>_<assembly name>`` in a directory listing (read, not reconstructed)."""
    names = set(re.findall(rf'href="({re.escape(accession)}_[^"/]+)/?"', listing_html))
    if len(names) != 1:
        raise ReferenceError(f"expected one folder for {accession} in the NCBI listing, found {sorted(names)}")
    return names.pop()


def parse_md5_listing(text: str) -> dict[str, str]:
    """``md5checksums.txt`` -> {file name: md5}."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[1].removeprefix("./")] = parts[0]
    return out


def _http_get_text(url: str) -> str:
    return get_text(url, timeout=TIMEOUT)


def _http_download(url: str, dest: Path) -> None:
    """Download via <dest>.part, retrying and resuming a partial file (the genome is ~0.2 GB)."""
    from whale.mirror import MirrorError, download

    try:
        download(url, dest)
    except MirrorError as exc:
        raise ReferenceError(f"download failed: {exc}") from None


def _decompress(src: Path, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)
    tmp.replace(dest)


# --- annotation ------------------------------------------------------------------

def gtf_header(gtf_path: Path) -> dict[str, str]:
    """The ``#!key value`` lines NCBI writes at the top of a GTF."""
    header = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("#!"):
                key, _, value = line[2:].strip().partition(" ")
                header[key] = value.strip()
    return header


def annotation_release_matches(registered: str | None, annotation_source: str | None) -> bool:
    """Whether the GTF's ``annotation-source`` names the registered release (ignoring the word RefSeq)."""
    if not registered or not annotation_source:
        return False

    def norm(text: str) -> str:
        return " ".join(w for w in text.split() if w.lower() != "refseq")

    return norm(registered) in norm(annotation_source)


_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
_GENEID_RE = re.compile(r"GeneID:(\d+)")


def _attrs(field: str) -> dict[str, str]:
    out = {}
    for k, v in _ATTR_RE.findall(field):
        if k == "db_xref":
            m = _GENEID_RE.search(v)
            if m:
                out["GeneID"] = m.group(1)
        else:
            out.setdefault(k, v)
    return out


def chrom_sizes(fasta: Path) -> list[tuple[str, int]]:
    """(sequence name, length) in FASTA order."""
    out, name, n = [], None, 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, n))
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        out.append((name, n))
    return out


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def build_features(gtf: Path, sizes: list[tuple[str, int]]) -> dict[str, list[tuple]]:
    """BED6 rows (chrom, start0, end, name, score, strand) per feature class; sorted.

    Gene name is the NCBI GeneID (falling back to the GTF gene_id).
    """
    length = dict(sizes)
    order = {c: i for i, (c, _) in enumerate(sizes)}
    genes, exons = [], set()
    tx_exons: dict[str, list[tuple]] = {}
    unknown_chroms = set()
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            chrom, start0, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            if chrom not in length:
                unknown_chroms.add(chrom)
                continue
            a = _attrs(f[8])
            gene = a.get("GeneID") or a.get("gene_id", "")
            if f[2] == "gene":
                genes.append((chrom, start0, end, gene, a.get("gene_biotype", "."), strand))
            else:
                exons.add((chrom, start0, end, gene, ".", strand))
                tx_exons.setdefault(a.get("transcript_id") or gene, []).append((chrom, start0, end, gene, strand))
    if unknown_chroms:
        raise ReferenceError(f"GTF names sequences absent from the genome: {sorted(unknown_chroms)[:5]}")

    introns = set()
    for parts in tx_exons.values():
        parts.sort(key=lambda p: p[1])
        for prev, nxt in zip(parts, parts[1:], strict=False):
            if nxt[1] > prev[2]:
                introns.add((prev[0], prev[2], nxt[1], prev[3], ".", prev[4]))

    promoters = []
    for chrom, s, e, gene, _, strand in genes:
        if strand == "-":
            ps, pe = e - PROMOTER_DOWNSTREAM, e + PROMOTER_UPSTREAM
        else:
            ps, pe = s - PROMOTER_UPSTREAM, s + PROMOTER_DOWNSTREAM
        ps, pe = max(0, ps), min(length[chrom], pe)
        if pe > ps:
            promoters.append((chrom, ps, pe, gene, ".", strand))

    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e, *_ in genes:
        by_chrom.setdefault(chrom, []).append((s, e))
    intergenic = []
    for chrom, size in sizes:
        pos = 0
        for s, e in _merge(by_chrom.get(chrom, [])):
            if s > pos:
                intergenic.append((chrom, pos, s, "intergenic", ".", "."))
            pos = max(pos, e)
        if pos < size:
            intergenic.append((chrom, pos, size, "intergenic", ".", "."))

    def key(row):
        return (order[row[0]], row[1], row[2], row[3], row[5])

    return {name: sorted(rows, key=key) for name, rows in {
        "genes": genes, "exons": list(exons), "introns": list(introns),
        "promoters": promoters, "intergenic": intergenic}.items()}


def write_bed(rows: list[tuple], path: Path) -> None:
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def build_tx2gene(gtf: Path) -> list[tuple[str, str]]:
    """(transcript_id, NCBI GeneID) for every GTF transcript feature, sorted, unique."""
    pairs = set()
    with open(gtf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 9 or f[2] != "transcript":
                continue
            a = _attrs(f[8])
            if a.get("transcript_id") and a.get("GeneID"):
                pairs.add((a["transcript_id"], a["GeneID"]))
    return sorted(pairs)


# --- plan, fetch, check ----------------------------------------------------------

def plan_bundle(repo: Path, assembly: str | None = None, *, rnaseq: bool = False,
                get_text=_http_get_text) -> dict:
    """What would be fetched, reading only NCBI's directory listing and md5checksums.txt."""
    record = resolve_assembly(repo, assembly)
    accession = record["ncbi_assembly_accession"]
    parent = accession_dir_url(accession)
    folder = find_assembly_folder(get_text(parent), accession)
    base = f"{parent}{folder}/"
    md5s = parse_md5_listing(get_text(f"{base}md5checksums.txt"))
    files = []
    for suffix, local in {**CORE_FILES, **(RNASEQ_FILES if rnaseq else {})}.items():
        remote = f"{folder}{suffix}"
        if remote not in md5s:
            raise ReferenceError(f"{remote} is not listed in {base}md5checksums.txt")
        files.append({"remote": remote, "url": f"{base}{remote}", "ncbi_md5": md5s[remote], "local": local})
    return {"assembly": record, "base_url": base, "folder": folder, "files": files}


def fetch_bundle(repo: Path, assembly: str | None = None, *, dest: Path | None = None,
                 rnaseq: bool = False, get_text=_http_get_text, download=_http_download,
                 today: str | None = None, log=lambda msg: None) -> dict:
    """Download, verify, derive and record a bundle. Returns the provenance dict.

    Compressed files already present with NCBI's MD5 are reused, so a re-run
    (or adding --rnaseq later) downloads only what is missing.
    """
    plan = plan_bundle(repo, assembly, rnaseq=rnaseq, get_text=get_text)
    record = plan["assembly"]
    dest = Path(dest) if dest else default_bundle_dir(repo, record["ncbi_assembly_accession"])
    dest.mkdir(parents=True, exist_ok=True)
    md5_path = dest / "md5checksums.txt"
    md5_path.write_text(get_text(f"{plan['base_url']}md5checksums.txt"))

    files = []
    for f in plan["files"]:
        gz = dest / f["remote"]
        if gz.exists() and md5_file(gz) == f["ncbi_md5"]:
            log(f"reuse    {f['remote']} (MD5 matches NCBI)")
        else:
            log(f"download {f['url']}")
            download(f["url"], gz)
        got = md5_file(gz)
        if got != f["ncbi_md5"]:
            gz.unlink(missing_ok=True)
            raise ReferenceError(f"{f['remote']}: MD5 {got} does not match NCBI's {f['ncbi_md5']}; deleted, re-run")
        local = dest / f["local"]
        log(f"unpack   {f['local']}")
        _decompress(gz, local)
        files.append({"file": f["local"], "sha256": sha256_file(local), "source_file": f["remote"],
                      "source_url": f["url"], "source_md5_ncbi": f["ncbi_md5"], "source_sha256": sha256_file(gz)})

    log("derive   chrom.sizes and features/")
    sizes = chrom_sizes(dest / "genomic.fna")
    (dest / "chrom.sizes").write_text("".join(f"{c}\t{n}\n" for c, n in sizes))
    files.append({"file": "chrom.sizes", "sha256": sha256_file(dest / "chrom.sizes"), "derived_from": "genomic.fna"})

    features = build_features(dest / "genomic.gtf", sizes)
    (dest / "features").mkdir(exist_ok=True)
    feature_counts = {}
    for name, rows in features.items():
        path = dest / "features" / f"{name}.bed"
        write_bed(rows, path)
        feature_counts[name] = len(rows)
        files.append({"file": f"features/{name}.bed", "sha256": sha256_file(path), "derived_from": "genomic.gtf"})

    tx2gene = None
    if rnaseq:
        pairs = build_tx2gene(dest / "genomic.gtf")
        if not pairs:
            raise ReferenceError("no transcript -> GeneID pairs in genomic.gtf")
        (dest / "tx2gene.tsv").write_text("".join(f"{t}\t{g}\n" for t, g in pairs))
        files.append({"file": "tx2gene.tsv", "sha256": sha256_file(dest / "tx2gene.tsv"), "derived_from": "genomic.gtf"})
        tx2gene = {"transcripts": len(pairs), "genes": len({g for _, g in pairs})}

    header = gtf_header(dest / "genomic.gtf")
    provenance = {
        "bundle": dest.name,
        "bundler": "whale fetch-reference (src/whale/reference.py)",
        "bundler_version": BUNDLER_VERSION,
        "date_generated": today or date.today().isoformat(),
        "assembly_id": record["assembly_id"],
        "species": record.get("species"),
        "ncbi_taxid": record.get("ncbi_taxid"),
        "ncbi_assembly_accession": record["ncbi_assembly_accession"],
        "registered_annotation_release": record.get("annotation_release"),
        "gtf_header": header,
        "annotation_release_matches_registry": annotation_release_matches(
            record.get("annotation_release"), header.get("annotation-source")),
        "source_base_url": plan["base_url"],
        "md5checksums_sha256": sha256_file(md5_path),
        "sequences": len(sizes),
        "genome_length": sum(n for _, n in sizes),
        "feature_counts": feature_counts,
        "feature_rules": FEATURE_RULES,
        "files": files,
        **({"tx2gene": tx2gene} if tx2gene else {}),
        "not_included": "Bismark index (built by the workflow); CpG index and C/T SNP mask (later Phase 1).",
    }
    (dest / PROVENANCE_NAME).write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    return provenance


def check_bundle(dest: Path) -> list[str]:
    """Problems with an existing bundle: missing files, checksum drift, wrong release. Empty = OK."""
    dest = Path(dest)
    prov_path = dest / PROVENANCE_NAME
    if not prov_path.exists():
        return [f"{prov_path} is missing; fetch the bundle with `whale fetch-reference`"]
    provenance = json.loads(prov_path.read_text())
    problems = []
    for entry in provenance["files"]:
        path = dest / entry["file"]
        if not path.exists():
            problems.append(f"{entry['file']}: missing")
        elif sha256_file(path) != entry["sha256"]:
            problems.append(f"{entry['file']}: SHA-256 differs from {PROVENANCE_NAME}")
    if not provenance.get("annotation_release_matches_registry"):
        problems.append(f"GTF annotation-source {provenance.get('gtf_header', {}).get('annotation-source')!r} "
                        f"does not name the registered release {provenance.get('registered_annotation_release')!r}")
    return problems
