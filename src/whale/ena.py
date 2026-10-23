"""Build a study's sample sheet and FASTQ manifest from its ENA/SRA BioProject.

The design of a study (which run belongs to which group, how many replicates
each group has) is read from the deposited record, the ENA read-run report
plus each BioSample's attributes, never transcribed from a paper. This writes,
into data/studies/<STUDY_ID>/:

* ``samplesheet.tsv``     one row per methylation run, in the column contract
                          workflows/methylation reads (sample_id, run_accession,
                          fastq_1, fastq_2), plus the design attributes;
* ``fastq_manifest.tsv``  every FASTQ URL with ENA's MD5 and byte count, so a
                          download is verified rather than trusted;
* ``ena_provenance.json`` the queries, the date, what was kept and dropped,
                          group sizes, and a checksum of each file written.

Nothing here downloads sequence data.

Differences from AREE (VENDORED_FROM_AREE.md):

* Runs are filtered to methylation library strategies by default, and every
  other strategy in the project (often the paired RNA-seq) is reported, not
  silently dropped.
* Every run must be Magallana / Crassostrea gigas (NCBI taxid 29159).
* sample_id is the run accession, which is unique; runs that share an ENA
  experiment (one library sequenced on several lanes) share a library_id.
* No methylKit sheet: differential analysis is downstream (PLAN.md §6.3).

Adapted from AREE's src/intake/ena_samplesheet.py (github.com/sr320/AREE).
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

from whale.net import get_text

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml"

RUN_FIELDS = [
    "run_accession", "sample_accession", "sample_alias", "experiment_accession",
    "experiment_title", "library_name", "library_strategy", "library_layout",
    "library_selection", "instrument_model", "scientific_name", "tax_id",
    "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
    "submitted_ftp", "submitted_md5", "submitted_bytes", "submitted_format",
]

# ENA library_strategy values that carry methylation data. EM-seq is usually
# deposited as Bisulfite-Seq; check the library protocol.
METHYLATION_STRATEGIES = ["Bisulfite-Seq", "MeDIP-Seq", "MBD-Seq", "MRE-Seq"]
EXPECTED_TAXID = "29159"

MATES_PER_LAYOUT = {"SINGLE": 1, "PAIRED": 2}
_FASTQ_NAME = re.compile(r"\.f(ast)?q(\.gz|\.bz2)?$", re.IGNORECASE)
_MATE_SUFFIX = re.compile(r"_R?([12])(_\d+)?\.f(ast)?q(\.gz|\.bz2)?$", re.IGNORECASE)

TIMEOUT = 120

SHEET_FIXED = ["sample_id", "run_accession", "fastq_1", "fastq_2", "library_id", "label",
               "sample_accession", "sample_alias", "condition"]
SHEET_TAIL = ["library_strategy", "library_layout", "instrument_model", "read_count"]


class ENAError(RuntimeError):
    """A query returned nothing usable, or the project is not what was expected."""


def _get(url: str) -> str:
    return get_text(url, timeout=TIMEOUT)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _run_fields(extra: list[str] | None = None) -> list[str]:
    return list(dict.fromkeys(RUN_FIELDS + (extra or [])))


def run_report_url(bioproject: str, extra_fields: list[str] | None = None) -> str:
    query = urllib.parse.urlencode(
        {"accession": bioproject, "result": "read_run",
         "fields": ",".join(_run_fields(extra_fields)), "format": "tsv"})
    return f"{ENA_PORTAL}?{query}"


def fetch_run_report(bioproject: str, extra_fields: list[str] | None = None) -> list[dict]:
    """One record per sequencing run, from ENA's read_run report."""
    text = _get(run_report_url(bioproject, extra_fields))
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ENAError(
            f"ENA returned no runs for {bioproject}. Either the accession is wrong, or the "
            "project has no public reads (some BioProjects register metadata only).")
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=False)) for ln in lines[1:]]


def fetch_sample_attributes(sample_accessions: list[str]) -> dict:
    """Attribute block per sample: the experimental factors, as deposited."""
    out: dict[str, dict] = {}
    for i in range(0, len(sample_accessions), 40):  # URL length limit
        batch = sample_accessions[i : i + 40]
        root = ET.fromstring(_get(f"{ENA_XML}/{','.join(batch)}"))
        for sample in root.iter("SAMPLE"):
            primary = sample.findtext(".//PRIMARY_ID")
            ext = [e.text for e in sample.iter("EXTERNAL_ID") if e.text]
            attrs = {a.findtext("TAG"): a.findtext("VALUE")
                     for a in sample.iter("SAMPLE_ATTRIBUTE") if a.findtext("TAG")}
            attrs["_alias"] = sample.get("alias")
            # The BioSample title often carries the design when attributes do not.
            attrs.setdefault("sample_title", sample.findtext("TITLE"))
            for key in filter(None, [primary, *ext]):
                out[key] = attrs
    return out


def _slug(value) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(value or "NA")).strip("_")


def _split(value) -> list[str]:
    return [v for v in (value or "").split(";") if v]


def run_fastq_files(run: dict) -> tuple[str | None, list[tuple[str, str, str]], list[str]]:
    """The FASTQ files for one run: (source, [(url, md5, bytes)], dropped urls).

    ENA's generated FASTQ (`fastq_ftp`) is preferred; when there is none, the
    submitter's files (`submitted_ftp`) are used if every one is a FASTQ.
    Anything else returns source None. A paired run with a third, unpaired
    file has it dropped and returned. Mates are ordered by their _1/_2 suffix.
    """
    def files(prefix):
        urls, md5s, sizes = (_split(run.get(f"{prefix}_{k}")) for k in ("ftp", "md5", "bytes"))
        return [(u, md5s[i] if i < len(md5s) else "", sizes[i] if i < len(sizes) else "")
                for i, u in enumerate(urls)]

    source, found = "fastq_ftp", files("fastq")
    if not found:
        formats = {f.upper() for f in _split(run.get("submitted_format"))}
        submitted = files("submitted")
        if submitted and formats <= {"FASTQ"} and all(_FASTQ_NAME.search(u) for u, _, _ in submitted):
            source, found = "submitted_ftp", submitted
        else:
            return None, [], []

    dropped: list[str] = []
    if run.get("library_layout") == "PAIRED" and len(found) == 3:
        unpaired = [f for f in found if not _MATE_SUFFIX.search(f[0])]
        if len(unpaired) == 1:
            dropped = [unpaired[0][0]]
            found = [f for f in found if f is not unpaired[0]]
    mates = [_MATE_SUFFIX.search(f[0]) for f in found]
    if len(found) > 1 and all(mates) and len({m.group(1) for m in mates}) == len(found):
        found = [f for _, f in sorted(zip((m.group(1) for m in mates), found, strict=True))]
    return source, found, dropped


def parse_run_labels(runs: list[dict], label_field: str, label_pattern: str) -> dict[str, dict]:
    """Per-run design attributes parsed from a run-level field.

    For projects that deposit every library under one BioSample, the design
    lives only in a run-level string (library_name, experiment_title). The
    regex must match every run; its named groups become attributes and the
    matched text the run's label. A run that does not match is an error.
    """
    rx = re.compile(label_pattern)
    if not rx.groupindex:
        raise ENAError("--label-pattern must define at least one named group, e.g. (?P<family>\\d+).")
    out, unmatched = {}, []
    for r in runs:
        m = rx.search(r.get(label_field) or "")
        if not m:
            unmatched.append(f"{r['run_accession']} ({label_field}={r.get(label_field)!r})")
            continue
        out[r["run_accession"]] = {**m.groupdict(), "_label": m.group(0)}
    if unmatched:
        raise ENAError(f"{len(unmatched)} run(s) do not match --label-pattern: {unmatched[:5]}")
    return out


def read_design_table(path: Path, join_attribute: str, attrs: dict) -> dict[str, dict]:
    """Per-sample design attributes from a curated table, joined on one attribute.

    For projects whose BioSamples do not carry the design, the groups come from
    the authors' own table, committed with the study. Every sample must match
    exactly one row and every row a sample.
    """
    with open(path, newline="") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    if not table or join_attribute not in table[0]:
        raise ENAError(f"--design-table {path} has no '{join_attribute}' column to join on.")
    keys = [row[join_attribute] for row in table]
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    if dupes:
        raise ENAError(f"--design-table {path} repeats {join_attribute} values: {dupes[:5]}")
    by_key = {row[join_attribute]: {k: v for k, v in row.items() if k != join_attribute} for row in table}

    out, unmatched = {}, []
    for sample, a in attrs.items():
        key = a.get(join_attribute)
        if key not in by_key:
            unmatched.append(f"{sample} ({join_attribute}={key!r})")
            continue
        out[sample] = by_key[key]
    if unmatched:
        raise ENAError(f"{len(unmatched)} sample(s) have no row in --design-table: {unmatched[:5]}")
    unused = sorted(set(by_key) - {attrs[s].get(join_attribute) for s in out})
    if unused:
        raise ENAError(f"--design-table rows match no sample in the project: {unused[:5]}")
    return out


def _rel(path: Path, repo: Path | None) -> str:
    path = Path(path).resolve()
    if repo is not None and path.is_relative_to(Path(repo).resolve()):
        return str(path.relative_to(Path(repo).resolve()))
    return str(path)


def build_samplesheet(
    bioproject: str,
    *,
    study_id: str,
    out_dir: Path,
    condition_attributes: list[str],
    extra_attributes: list[str] | None = None,
    label_field: str | None = None,
    label_pattern: str | None = None,
    library_strategies: list[str] | None = None,
    design_table: Path | None = None,
    design_join: str | None = None,
    design_source: str | None = None,
    repo: Path | None = None,
) -> dict:
    """Write samplesheet + FASTQ manifest + provenance for one BioProject.

    `condition_attributes` name the sample attributes that jointly define a
    group; their slugified concatenation is the `condition` column. Naming
    them keeps the design decision in the curation record. Named groups of
    `label_pattern` and columns of `design_table` are available as attributes
    too and override deposited ones of the same name.

    `library_strategies` (default: the methylation strategies) selects the runs
    kept. Every strategy in the project, kept or not, is reported with its run
    count, so paired RNA-seq deposited in the same project is visible.

    Returns the provenance record plus report-only keys for the CLI.
    """
    out_dir = Path(out_dir)
    extra_attributes = extra_attributes or []
    if bool(label_field) != bool(label_pattern):
        raise ENAError("--label-field and --label-pattern must be given together.")
    if bool(design_table) != bool(design_join):
        raise ENAError("--design-table and --design-join must be given together.")
    keep = list(library_strategies or METHYLATION_STRATEGIES)
    explicit = bool(library_strategies)

    all_runs = fetch_run_report(bioproject, [label_field] if label_field else None)
    strategies_present: dict[str, int] = {}
    for r in all_runs:
        strategies_present[r["library_strategy"]] = strategies_present.get(r["library_strategy"], 0) + 1
    if explicit:
        unknown = set(keep) - set(strategies_present)
        if unknown:
            raise ENAError(f"--library-strategy {sorted(unknown)} not in this project; it has "
                           f"{dict(sorted(strategies_present.items()))}")
    runs = [r for r in all_runs if r["library_strategy"] in keep]
    if not runs:
        raise ENAError(
            f"{bioproject} has no runs with library_strategy in {keep}; it has "
            f"{dict(sorted(strategies_present.items()))}. If its methylation runs are labelled "
            "otherwise (e.g. OTHER for EM-seq), pass --library-strategy explicitly.")

    wrong_taxon = sorted(f"{r['run_accession']} ({r.get('scientific_name')}, taxid {r.get('tax_id')})"
                         for r in runs if r.get("tax_id") != EXPECTED_TAXID)
    if wrong_taxon:
        raise ENAError(f"{len(wrong_taxon)} run(s) are not Magallana / Crassostrea gigas "
                       f"(taxid {EXPECTED_TAXID}): {wrong_taxon[:5]}")

    attrs = fetch_sample_attributes(sorted({r["sample_accession"] for r in runs}))
    run_labels = parse_run_labels(runs, label_field, label_pattern) if label_field else {}
    design = read_design_table(Path(design_table), design_join, attrs) if design_table else {}

    all_attr_cols = list(dict.fromkeys(condition_attributes + extra_attributes))
    attr_cols = [_slug(c).lower() for c in all_attr_cols]
    clash = sorted(set(attr_cols) & set(SHEET_FIXED + SHEET_TAIL))
    if clash:
        raise ENAError(f"attribute name(s) {clash} collide with sample-sheet columns; rename in a design table.")
    header = SHEET_FIXED + attr_cols + SHEET_TAIL

    rows, manifest_rows = [], []
    missing, no_fastq, wrong_count = [], [], []
    sources: dict[str, str] = {}
    unpaired_dropped: dict[str, list[str]] = {}
    for r in sorted(runs, key=lambda x: x["run_accession"]):
        a = {**attrs.get(r["sample_accession"], {}), **design.get(r["sample_accession"], {}),
             **run_labels.get(r["run_accession"], {})}
        values = [a.get(c) for c in condition_attributes]
        if any(v is None for v in values):
            missing.append(f"{r['run_accession']} ({r['sample_accession']})")
        condition = "_".join(_slug(v) for v in values)

        source, files, dropped = run_fastq_files(r)
        if source is None:
            no_fastq.append(f"{r['run_accession']} (submitted_format={r.get('submitted_format') or ''!r})")
            continue
        expected = MATES_PER_LAYOUT.get(r["library_layout"])
        if len(files) != expected:
            wrong_count.append(f"{r['run_accession']} ({r['library_layout']}: {len(files)} file(s) "
                               f"from {source}, expected {expected})")
            continue
        sources[r["run_accession"]] = source
        if dropped:
            unpaired_dropped[r["run_accession"]] = [f"ftp://{u}" for u in dropped]
        names = [u.rsplit("/", 1)[-1] for u, _, _ in files]
        for idx, (url, md5, size) in enumerate(files):
            manifest_rows.append({"run_accession": r["run_accession"], "mate": str(idx + 1),
                                  "file": names[idx], "url": f"ftp://{url}", "md5": md5,
                                  "bytes": size, "source": source})
        rows.append({
            "sample_id": r["run_accession"], "run_accession": r["run_accession"],
            "fastq_1": names[0], "fastq_2": names[1] if len(names) > 1 else "",
            "library_id": r.get("experiment_accession") or r["run_accession"],
            "label": a.get("_label") or r.get("library_name") or r.get("sample_alias") or "",
            "sample_accession": r["sample_accession"], "sample_alias": r.get("sample_alias") or "",
            "condition": condition,
            **{col: str(a.get(c) or "") for col, c in zip(attr_cols, all_attr_cols, strict=True)},
            "library_strategy": r["library_strategy"], "library_layout": r["library_layout"],
            "instrument_model": r.get("instrument_model") or "", "read_count": r.get("read_count") or "",
        })

    if missing:
        raise ENAError(
            f"{len(missing)} run(s) lack one of the condition attributes {condition_attributes}: "
            f"{missing[:5]}. Attribute names are submitter-defined; list them with --show-attributes.")
    if no_fastq:
        raise ENAError(f"{len(no_fastq)} run(s) have no usable FASTQ (fastq_ftp empty, submitted_ftp "
                       f"missing or not FASTQ): {no_fastq}")
    if wrong_count:
        raise ENAError(f"{len(wrong_count)} run(s) do not have one FASTQ per mate: {wrong_count}")

    out_dir.mkdir(parents=True, exist_ok=True)
    sheet_path = out_dir / "samplesheet.tsv"
    with sheet_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    manifest_path = out_dir / "fastq_manifest.tsv"
    with manifest_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["run_accession", "mate", "file", "url", "md5", "bytes", "source"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(manifest_rows)

    group_sizes: dict[str, int] = {}
    group_libraries: dict[str, set] = {}
    for r in rows:
        group_sizes[r["condition"]] = group_sizes.get(r["condition"], 0) + 1
        group_libraries.setdefault(r["condition"], set()).add(r["library_id"])
    runs_per_library: dict[str, int] = {}
    for r in rows:
        runs_per_library[r["library_id"]] = runs_per_library.get(r["library_id"], 0) + 1
    multi_run = {k: v for k, v in sorted(runs_per_library.items()) if v > 1}
    samples_per_biosample: dict[str, int] = {}
    for lib in {(r["sample_accession"], r["library_id"]) for r in rows}:
        samples_per_biosample[lib[0]] = samples_per_biosample.get(lib[0], 0) + 1
    shared_biosamples = {k: v for k, v in sorted(samples_per_biosample.items()) if v > 1}

    provenance = {
        "study_id": study_id,
        "bioproject": bioproject,
        "date_generated": date.today().isoformat(),
        "generated_by": "whale fetch-samplesheet (src/whale/ena.py)",
        "queries": {"run_report": run_report_url(bioproject, [label_field] if label_field else None),
                    "sample_xml": f"{ENA_XML}/<sample_accessions>"},
        "condition_attributes": condition_attributes,
        "library_strategies_kept": sorted(keep if explicit else set(keep) & set(strategies_present)),
        "runs_by_strategy_in_project": dict(sorted(strategies_present.items())),
        "n_runs": len(rows),
        "n_libraries": len(runs_per_library),
        "n_biosamples": len({r["sample_accession"] for r in rows}),
        "library_layouts": sorted({r["library_layout"] for r in rows}),
        "instrument_models": sorted({r["instrument_model"] for r in rows}),
        "group_sizes_runs": dict(sorted(group_sizes.items())),
        "group_sizes_libraries": {k: len(v) for k, v in sorted(group_libraries.items())},
        "min_group_libraries": min(len(v) for v in group_libraries.values()),
        "total_fastq_bytes": sum(int(m["bytes"]) for m in manifest_rows if m["bytes"]),
        "fastq_sources": {k: sum(v == k for v in sources.values()) for k in ("fastq_ftp", "submitted_ftp")},
        "outputs": {_rel(p, repo): sha256_file(p) for p in (sheet_path, manifest_path)},
    }
    if multi_run:
        provenance["libraries_with_several_runs"] = multi_run
    if shared_biosamples:
        provenance["biosamples_with_several_libraries"] = shared_biosamples
    if label_field:
        provenance["run_labels"] = {"field": label_field, "pattern": label_pattern,
                                    "labels": {k: v["_label"] for k, v in sorted(run_labels.items())}}
    if design_table:
        provenance["design_table"] = {"path": _rel(Path(design_table), repo),
                                      "sha256": sha256_file(Path(design_table)),
                                      "join_attribute": design_join, "source": design_source}
    if unpaired_dropped:
        provenance["unpaired_files_dropped"] = dict(sorted(unpaired_dropped.items()))
    other = {k: v for k, v in strategies_present.items() if k not in keep}
    if other:
        provenance["runs_not_kept"] = dict(sorted(other.items()))

    prov_path = out_dir / "ena_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return {**provenance, "provenance_file": str(prov_path)}


def sample_attribute_summary(bioproject: str, library_strategies: list[str] | None = None) -> dict:
    """attribute -> sorted distinct values, across the project's kept runs' BioSamples.

    For choosing --condition-attribute: shows what the submitters deposited.
    """
    keep = list(library_strategies or METHYLATION_STRATEGIES)
    runs = [r for r in fetch_run_report(bioproject) if r["library_strategy"] in keep]
    attrs = fetch_sample_attributes(sorted({r["sample_accession"] for r in runs}))
    seen = {}
    for s in {r["sample_accession"] for r in runs}:
        for k, v in attrs.get(s, {}).items():
            if not k.startswith("_"):
                seen.setdefault(k, set()).add(v or "")
    summary = {k: sorted(v) for k, v in sorted(seen.items())}
    for field in ("library_name", "experiment_title", "sample_alias"):
        summary[f"(run) {field}"] = sorted({r.get(field) or "" for r in runs})
    return summary
