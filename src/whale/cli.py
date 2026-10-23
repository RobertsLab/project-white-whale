"""The `whale` command line."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import click
import yaml

from whale import __version__
from whale.claims import (
    ClaimError,
    list_claims,
    parse_branch,
    registered_on,
    release_study,
    set_status,
    start_study,
)
from whale.paths import find_repo
from whale.registry import build_index, load_backlog, load_studies, study_files
from whale.scope import check_scope
from whale.validate import validate_status, validate_study_file
from whale.vocab import load_vocab


def _repo() -> Path:
    try:
        return find_repo()
    except FileNotFoundError as exc:
        raise click.ClickException(str(exc)) from None


def _table(rows: list[dict], columns: list[str]) -> str:
    if not rows:
        return "(none)"
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in columns}
    lines = ["  ".join(c.ljust(widths[c]) for c in columns),
             "  ".join("-" * widths[c] for c in columns)]
    lines += ["  ".join(str(r.get(c, "")).ljust(widths[c]) for c in columns) for r in rows]
    return "\n".join(lines)


@click.group()
@click.version_option(__version__, prog_name="whale")
def main():
    """White Whale: public DNA methylation of the Pacific oyster, reprocessed on one footing."""


# --- registry ----------------------------------------------------------------

@main.command("validate-study")
@click.argument("files", nargs=-1, type=click.Path(exists=True, dir_okay=False))
@click.option("--strict", is_flag=True, help="Treat warnings (fields not yet filled in) as errors.")
def validate_study_cmd(files, strict):
    """Validate study YAMLs (default: every file in registry/studies/)."""
    repo = _repo()
    paths = [Path(f) for f in files] or study_files(repo)
    failed = False
    for path in paths:
        res = validate_study_file(path, repo)
        bad = res.errors or (strict and res.warnings)
        click.echo(f"{'FAIL' if bad else 'ok  '}  {path}")
        for e in res.errors:
            click.echo(click.style(f"      error: {e}", fg="red"))
        for w in res.warnings:
            click.echo(click.style(f"      warning: {w}", fg="yellow"))
        failed |= bool(bad)
    if not paths:
        click.echo("No study YAMLs in registry/studies/ yet.")
    sys.exit(1 if failed else 0)


@main.command("validate-status")
@click.argument("files", nargs=-1, type=click.Path(exists=True, dir_okay=False))
def validate_status_cmd(files):
    """Validate status.yaml files (default: every one under data/studies/)."""
    repo = _repo()
    paths = [Path(f) for f in files] or sorted((repo / "data" / "studies").rglob("status.yaml"))
    failed = False
    for path in paths:
        errors = validate_status(yaml.safe_load(path.read_text()) or {}, repo)
        click.echo(f"{'FAIL' if errors else 'ok  '}  {path}")
        for e in errors:
            click.echo(click.style(f"      error: {e}", fg="red"))
        failed |= bool(errors)
    sys.exit(1 if failed else 0)


@main.command("build-registry")
@click.option("--out", type=click.Path(dir_okay=False), default=None)
def build_registry_cmd(out):
    """Rebuild registry/study_registry.csv (not committed) from every study YAML."""
    repo = _repo()
    path = build_index(repo, Path(out) if out else None)
    click.echo(f"Wrote {path} ({len(load_studies(repo))} studies)")


@main.command("backlog")
@click.option("--no-fetch", is_flag=True, help="Use the remote branches as last fetched.")
@click.option("--remote", default="origin", show_default=True)
def backlog_cmd(no_fetch, remote):
    """The study backlog, with what is registered, claimed (where, at what stage) or free."""
    repo = _repo()
    registered = load_studies(repo)
    registered_accs = {(s.get("accessions") or {}).get("bioproject"): sid for sid, s in registered.items()}
    try:
        claims = list_claims(repo, remote=remote, fetch=not no_fetch)
    except ClaimError as exc:
        click.echo(click.style(f"Could not read claims from {remote}: {exc}", fg="yellow"), err=True)
        claims = []
    by_acc: dict[str, list[dict]] = {}
    for c in claims:
        by_acc.setdefault(c["accession"], []).append(c)

    rows = []
    for b in load_backlog(repo):
        acc = b["accession"]
        held = by_acc.get(acc, [])
        if held:
            state = "; ".join(f"{c['stage']} on {c['machine']}" + (f" ({c['part']})" if c["part"] else "")
                              for c in held)
        elif acc in registered_accs:
            state = f"registered as {registered_accs[acc]}"
        else:
            state = "free"
        rows.append({"accession": acc, "study": b["suggested_study_id"], "chem": b["chemistry"],
                     "runs": b["runs"], "GB": b["fastq_gb"], "rnaseq": b["paired_rnaseq"],
                     "profile": b["suggested_profile"], "state": state})
    click.echo(_table(rows, ["accession", "study", "chem", "runs", "GB", "rnaseq", "profile", "state"]))


@main.command("fetch-samplesheet")
@click.option("--bioproject", required=True, help="BioProject accession, e.g. PRJNA682817.")
@click.option("--study", "study_id", default=None,
              help="study_id the sheet is for (default: the study of the current study branch).")
@click.option("--show-attributes", is_flag=True,
              help="Only list the deposited sample attributes and run labels, to choose --condition-attribute.")
@click.option("--condition-attribute", "condition_attributes", multiple=True,
              help="Sample attribute defining an experimental group (repeatable; combined in order into "
                   "`condition`). May name a --label-pattern group, a --design-table column, or sample_title.")
@click.option("--attribute", "extra_attributes", multiple=True,
              help="Additional sample attribute to carry into the sheet (repeatable).")
@click.option("--label-field", default=None,
              help="Run-report field holding per-library design labels (e.g. library_name, experiment_title).")
@click.option("--label-pattern", default=None,
              help="Regex with named groups applied to --label-field; each group becomes an attribute.")
@click.option("--library-strategy", "library_strategies", multiple=True,
              help="Keep runs with this ENA library_strategy (repeatable). Default: the methylation "
                   "strategies (Bisulfite-Seq, MeDIP-Seq, MBD-Seq, MRE-Seq).")
@click.option("--design-table", default=None, type=click.Path(exists=True, dir_okay=False),
              help="TSV of design attributes for projects whose BioSamples lack the design; committed with the study.")
@click.option("--design-join", default=None, help="Deposited sample attribute matching the design table's column of that name.")
@click.option("--design-source", default=None, help="Where the design table came from; recorded in provenance.")
@click.option("--out", "out_dir", default=None, type=click.Path(),
              help="Output directory (default: data/studies/<study_id>/).")
def fetch_samplesheet_cmd(bioproject, study_id, show_attributes, condition_attributes, extra_attributes,
                          label_field, label_pattern, library_strategies, design_table, design_join,
                          design_source, out_dir):
    """Build a study's sample sheet and FASTQ manifest from its deposited ENA metadata.

    Downloads no sequence data. Writes samplesheet.tsv, fastq_manifest.tsv (ENA's
    MD5s) and ena_provenance.json into data/studies/<study_id>/.
    """
    from whale.claims import ClaimError as _CE
    from whale.claims import current_branch
    from whale.ena import ENAError, build_samplesheet, sample_attribute_summary

    repo = _repo()
    try:
        if show_attributes:
            for k, vals in sample_attribute_summary(bioproject, list(library_strategies) or None).items():
                shown = ", ".join(v or "''" for v in vals[:8]) + (f", … ({len(vals)} values)" if len(vals) > 8 else "")
                click.echo(f"{k:28} {shown}")
            return
        if not condition_attributes:
            raise click.UsageError("--condition-attribute is required (use --show-attributes to see the choices).")
        if study_id is None:
            try:
                parsed = parse_branch(current_branch(repo))
            except _CE:
                parsed = None
            if parsed is None:
                raise click.UsageError("not on a study branch; pass --study.")
            if parsed[0] != bioproject:
                raise click.UsageError(f"this branch is for {parsed[0]}, not {bioproject}; pass --study to override.")
            study_id = parsed[1]
        report = build_samplesheet(
            bioproject, study_id=study_id,
            out_dir=Path(out_dir) if out_dir else repo / "data" / "studies" / study_id,
            condition_attributes=list(condition_attributes), extra_attributes=list(extra_attributes),
            label_field=label_field, label_pattern=label_pattern,
            library_strategies=list(library_strategies) or None,
            design_table=Path(design_table) if design_table else None,
            design_join=design_join, design_source=design_source, repo=repo)
    except ENAError as exc:
        raise click.ClickException(str(exc)) from None

    click.echo(f"{report['bioproject']} -> {study_id}: {report['n_runs']} runs, "
               f"{report['n_libraries']} libraries, {report['n_biosamples']} BioSamples")
    click.echo(f"  kept       : {', '.join(report['library_strategies_kept'])}")
    click.echo(f"  layouts    : {', '.join(report['library_layouts'])}")
    click.echo(f"  platforms  : {', '.join(report['instrument_models'])}")
    click.echo(f"  download   : {report['total_fastq_bytes'] / 1e9:.1f} GB of FASTQ")
    click.echo("  groups (libraries / runs):")
    for group, n in report["group_sizes_libraries"].items():
        click.echo(f"    {group:32} {n} / {report['group_sizes_runs'][group]}")
    if report.get("runs_not_kept"):
        click.echo(click.style(f"\nNot kept: {report['runs_not_kept']}. If these are RNA-seq of the same "
                               "animals, record them under paired_rnaseq in the study YAML.", fg="yellow"))
    if report.get("libraries_with_several_runs"):
        click.echo(click.style(f"\n{len(report['libraries_with_several_runs'])} libraries were sequenced in "
                               "several runs. The workflow processes each run separately for now.", fg="yellow"))
    if report.get("biosamples_with_several_libraries"):
        click.echo(click.style(f"\n{len(report['biosamples_with_several_libraries'])} BioSample(s) carry more "
                               "than one library: check whether the design is in a run-level field "
                               "(--label-field) and add the shared_biosample quality flag.", fg="yellow"))
    if report["min_group_libraries"] < 3:
        click.echo(click.style(f"\nSmallest group has {report['min_group_libraries']} libraries: "
                               "consider the low_replication quality flag.", fg="yellow"))
    click.echo(f"\nWrote {', '.join(report['outputs'])} and {report['provenance_file']}")
    click.echo(f"Fill in the study YAML: samples.n_runs={report['n_runs']}, "
               f"samples.n_libraries={report['n_libraries']}, "
               f"data_availability.fastq_size_gb={report['total_fastq_bytes'] / 1e9:.1f}")


@main.command("fetch-reference")
@click.option("--assembly", default=None,
              help="Assembly id or accession from data/reference/genome_assemblies.yaml (default: the default one).")
@click.option("--dest", default=None, type=click.Path(file_okay=False),
              help="Bundle directory (default: $WHALE_REFERENCE/<ACCESSION>, else data/reference/<ACCESSION>).")
@click.option("--rnaseq", is_flag=True, help="Also fetch the transcript FASTA and build tx2gene.tsv (Goal 2).")
@click.option("--plan", "plan_only", is_flag=True, help="Show what would be fetched, with NCBI MD5s; download nothing.")
@click.option("--check", "check_only", is_flag=True, help="Verify an existing bundle against its provenance.")
def fetch_reference_cmd(assembly, dest, rnaseq, plan_only, check_only):
    """Fetch, MD5-verify and derive the reference bundle (genome, GTF, chrom.sizes, feature BEDs)."""
    from whale.reference import (
        ReferenceError,
        check_bundle,
        default_bundle_dir,
        fetch_bundle,
        plan_bundle,
        resolve_assembly,
    )

    repo = _repo()
    try:
        record = resolve_assembly(repo, assembly)
        target = Path(dest) if dest else default_bundle_dir(repo, record["ncbi_assembly_accession"])
        if check_only:
            problems = check_bundle(target)
            for p in problems:
                click.echo(click.style(f"problem: {p}", fg="red"))
            if problems:
                sys.exit(1)
            click.echo(f"ok: {target} matches its provenance")
            return
        if plan_only:
            plan = plan_bundle(repo, assembly, rnaseq=rnaseq)
            click.echo(f"{record['assembly_id']} ({record['ncbi_assembly_accession']}) -> {target}")
            for f in plan["files"]:
                click.echo(f"  {f['local']:12} {f['ncbi_md5']}  {f['url']}")
            return
        prov = fetch_bundle(repo, assembly, dest=target, rnaseq=rnaseq, log=click.echo)
    except ReferenceError as exc:
        raise click.ClickException(str(exc)) from None

    click.echo(f"\n{prov['assembly_id']}: {prov['sequences']} sequences, {prov['genome_length'] / 1e6:.1f} Mb")
    for name, n in prov["feature_counts"].items():
        click.echo(f"  {name:11} {n:>9,}")
    if not prov["annotation_release_matches_registry"]:
        click.echo(click.style(
            f"\nThe GTF says {prov['gtf_header'].get('annotation-source')!r}, not the registered "
            f"{prov['registered_annotation_release']!r}. Update genome_assemblies.yaml before using this bundle.",
            fg="red"))
        sys.exit(1)
    click.echo(f"\nWrote {target}/reference_provenance.json")


@main.command("fetch-fastq")
@click.option("--study", "study_id", default=None, help="Default: the study of the current study branch.")
@click.option("--manifest", default=None, type=click.Path(exists=True, dir_okay=False),
              help="Default: data/studies/<study>/fastq_manifest.tsv.")
@click.option("--dest", default=None, type=click.Path(file_okay=False),
              help="Default: $WHALE_RAW/<study>, else data/raw/<study>.")
@click.option("--run", "runs", multiple=True, help="Only this run accession (repeatable).")
@click.option("--jobs", default=2, show_default=True, help="Files downloaded at once.")
@click.option("--dry-run", is_flag=True, help="Report what is staged and what is left; download nothing.")
@click.option("--ignore-space", is_flag=True, help="Start even if free space looks too small.")
@click.option("--mark-staged", is_flag=True,
              help="When every file is verified, set the study's stage to fastq_staged (commits and pushes).")
@click.option("--source", type=click.Choice(["auto", "mirror", "ena"]), default="auto", show_default=True,
              help="auto: the gannet mirror where it has a file, ENA otherwise; mirror: refuse unless every "
                   "file is mirrored; ena: ENA only.")
@click.option("--accession", default=None,
              help="BioProject whose mirror to use (default: the study branch's, else ena_provenance.json's).")
@click.option("--mirror-url", default=None,
              help="Mirror folder URL (default: data_host.yaml url_base + /raw/<ACCESSION>).")
def fetch_fastq_cmd(study_id, manifest, dest, runs, jobs, dry_run, ignore_space, mark_staged,
                    source, accession, mirror_url):
    """Download a study's FASTQ from its manifest, resuming and MD5-verifying every file.

    Files already on the gannet mirror are copied from there (checked against the
    mirror ledger's MD5); the rest, and any mirror copy that fails, come from ENA.
    """
    import json
    import os

    from whale.claims import ClaimError as _CE
    from whale.claims import current_branch
    from whale.fastq import (
        FastqError,
        apply_mirror,
        fetch_all,
        load_ledger,
        plan,
        read_manifest,
        record_staged,
    )
    from whale.mirror import LEDGER_NAME, load_data_host

    repo = _repo()
    try:
        parsed = parse_branch(current_branch(repo))
    except _CE:
        parsed = None
    if study_id is None:
        if parsed is None:
            raise click.UsageError("not on a study branch; pass --study.")
        study_id = parsed[1]
    study_dir = repo / "data" / "studies" / study_id
    manifest_path = Path(manifest) if manifest else study_dir / "fastq_manifest.tsv"
    if not manifest_path.exists():
        raise click.ClickException(f"{manifest_path} is missing; run whale fetch-samplesheet first.")
    raw = os.environ.get("WHALE_RAW")
    target = Path(dest) if dest else (Path(raw) if raw else repo / "data" / "raw") / study_id

    try:
        files = read_manifest(manifest_path, list(runs) or None)
    except FastqError as exc:
        raise click.ClickException(str(exc)) from None

    if source != "ena":
        if accession is None and parsed is not None and parsed[1] == study_id:
            accession = parsed[0]
        prov = manifest_path.parent / "ena_provenance.json"
        if accession is None and prov.exists():
            accession = json.loads(prov.read_text()).get("bioproject")
        if mirror_url is None and accession:
            try:
                base = load_data_host(repo).get("url_base")
            except FileNotFoundError:
                base = None
            mirror_url = f"{base.rstrip('/')}/raw/{accession}" if base else None
        if mirror_url is None:
            if source == "mirror":
                raise click.ClickException("no mirror to use: pass --accession or --mirror-url.")
            click.echo("No mirror location known (no accession); using ENA.")
        else:
            try:
                ledger = load_ledger(f"{mirror_url}/{LEDGER_NAME}")
            except (OSError, ValueError) as exc:
                if source == "mirror":
                    raise click.ClickException(f"could not read the mirror ledger: {exc}") from None
                click.echo(click.style(f"Could not read the mirror ledger ({exc}); using ENA.", fg="yellow"))
                ledger = None
            files, counts = apply_mirror(files, ledger, mirror_url)
            click.echo(f"Mirror {mirror_url}: {counts['mirror']} file(s) from the mirror, "
                       f"{counts['not_mirrored']} not mirrored, {counts['stale']} stale (ENA has changed them)")
            if source == "mirror" and counts["mirror"] < len(files):
                raise click.ClickException(f"{len(files) - counts['mirror']} file(s) are not usable from the "
                                           "mirror; use --source auto to take them from ENA.")

    p = plan(files, target)
    gb = 1e9
    click.echo(f"{study_id}: {p['files']} files, {p['total_bytes'] / gb:.1f} GB -> {target}")
    click.echo(f"  verified already : {p['staged']}")
    click.echo(f"  to download      : {len(p['todo'])} ({p['remaining_bytes'] / gb:.1f} GB left)")
    click.echo(f"  free space       : {p['free_bytes'] / gb:.1f} GB")
    short = p["remaining_bytes"] * 1.02 > p["free_bytes"]
    if short:
        click.echo(click.style("  Not enough free space for what is left.", fg="red"))
    if dry_run:
        return
    if short and not ignore_space:
        raise click.ClickException("not enough free space; free some or pass --ignore-space.")

    try:
        result = fetch_all(files, target, jobs=jobs, log=click.echo)
    except FastqError as exc:
        raise click.ClickException(str(exc)) from None
    # Beside the manifest, so a --manifest outside the repo never writes into it.
    record = manifest_path.parent / "fastq_staged.tsv"
    record_staged(files, target, record)
    click.echo(f"\n{result['downloaded']} downloaded, {result['skipped']} already verified, "
               f"{len(result['failed'])} failed. Recorded in {record}.")
    if result["failed"]:
        click.echo(click.style("Re-run the same command to resume the failed files.", fg="yellow"))
        sys.exit(1)
    if mark_staged:
        try:
            total = sum(f.bytes for f in files)
            # Commit the verification record with the stage change.
            if record.resolve().is_relative_to(repo.resolve()):
                subprocess.run(["git", "add", str(record.resolve().relative_to(repo.resolve()))],
                               cwd=repo, check=True)
            new = set_status("fastq_staged", repo=repo, note=f"{len(files)} files, {total / gb:.1f} GB, MD5 verified")
        except ClaimError as exc:
            raise click.ClickException(f"downloads verified, but the stage was not set: {exc}") from None
        click.echo(f"Stage: {new['stage']} on {new['machine']}")


@main.command("mirror-fastq")
@click.option("--accession", "accessions", multiple=True,
              help="BioProject to mirror (repeatable; default: every backlog entry, in backlog order).")
@click.option("--dest", default=None, type=click.Path(file_okay=False),
              help="Mirror into this local directory instead of the data host's upload_target.")
@click.option("--staging", default=None, type=click.Path(file_okay=False),
              help="Where files are downloaded and checked before upload "
                   "(default: $WHALE_MIRROR_STAGING, else $WHALE_RAW/.mirror-staging).")
@click.option("--source", type=click.Choice(["ncbi", "ena"]), default="ncbi", show_default=True,
              help="ncbi: .sra from NCBI's S3 bucket converted with fasterq-dump, ENA as fallback per run; "
                   "ena: ENA's FASTQ only.")
@click.option("--jobs", default=4, show_default=True, help="Runs fetched at once.")
@click.option("--threads", default=4, show_default=True, help="fasterq-dump / pigz threads per run (ncbi).")
@click.option("--tmp", default=None, type=click.Path(file_okay=False),
              help="fasterq-dump scratch, ideally fast local disk (default: <staging>/tmp).")
@click.option("--limit-rate", "rate_mb", default=None, type=float, help="Total download rate cap, MB/s.")
@click.option("--max-gb", default=None, type=float,
              help="Stop starting new files once this much has been queued in this run.")
@click.option("--library-strategy", "library_strategies", multiple=True,
              help="For unregistered projects: mirror runs with this ENA library_strategy (repeatable; "
                   "default the methylation strategies). Registered studies use their fastq_manifest.tsv.")
@click.option("--no-verify-remote", is_flag=True, help="Skip re-computing each MD5 on the data host.")
@click.option("--plan", "plan_only", is_flag=True, help="Show what is mirrored and what is left; download nothing.")
def mirror_fastq_cmd(accessions, dest, staging, source, jobs, threads, tmp, rate_mb, max_gb,
                     library_strategies, no_verify_remote, plan_only):
    """Prefetch backlog FASTQ to the data host (gannet), verified and resumable.

    By default each run is taken from NCBI as .sra (MD5-checked), converted with
    fasterq-dump and checked against ENA's read and base counts; ENA's FASTQ
    (MD5-checked) is the fallback. Safe to rerun, e.g. nightly from cron on
    raven: files already mirrored are skipped, partial downloads resume, and a
    second copy started while one is running exits at once.
    """
    from whale.mirror import (
        LocalStore,
        MirrorError,
        SSHStore,
        backlog_accessions,
        default_staging,
        load_data_host,
        mirror,
        plan_mirror,
    )

    repo = _repo()
    accs = list(accessions) or backlog_accessions(repo)
    try:
        store = LocalStore(Path(dest)) if dest else SSHStore(load_data_host(repo)["upload_target"])
        if plan_only:
            plans = plan_mirror(repo, accs, store, library_strategies=list(library_strategies) or None)
        else:
            res = mirror(repo, accs, store, staging=Path(staging) if staging else default_staging(repo),
                         jobs=jobs, rate_mb=rate_mb, max_gb=max_gb, verify_remote=not no_verify_remote,
                         source=source, threads=threads, tmp=Path(tmp) if tmp else None,
                         library_strategies=list(library_strategies) or None, log=click.echo)
            if res["locked"]:
                click.echo("Another mirror-fastq is running with this staging directory; nothing to do.")
                return
            plans = res["plans"]
    except MirrorError as exc:
        raise click.ClickException(str(exc)) from None

    just_done = set() if plan_only else set(res["mirrored"])
    rows = []
    for p in plans:
        todo = [f for f in p["todo"] if f"{p['accession']}/{f['file']}" not in just_done]
        if p.get("error"):
            rows.append({"accession": p["accession"], "files": "", "GB": "", "mirrored": "", "left GB": "",
                         "state": f"error: {p['error'][:60]}"})
            continue
        total = sum(int(f["bytes"] or 0) for f in p["files"]) / 1e9
        left = sum(int(f["bytes"] or 0) for f in todo) / 1e9
        done = len(p["files"]) - len(todo)
        rows.append({"accession": p["accession"], "files": len(p["files"]), "GB": f"{total:.1f}",
                     "mirrored": done, "left GB": f"{left:.1f}",
                     "state": "complete" if p["files"] and not todo else ("no files" if not p["files"] else "")})
    if plan_only:
        click.echo(f"Mirror at {store.label}\n")
    else:
        click.echo(f"\nMirrored {len(res['mirrored'])} file(s), {res['bytes'] / 1e9:.1f} GB this run "
                   f"({res['sources']['ncbi_sra']} converted from NCBI, {res['sources']['ena']} from ENA).\n")
    click.echo(_table(rows, ["accession", "files", "GB", "mirrored", "left GB", "state"]))
    for p in plans:
        for n in p.get("notes", []):
            click.echo(click.style(f"  {p['accession']}: {n}", fg="yellow"))
    if not plan_only and res["failed"]:
        click.echo(click.style(f"\n{len(res['failed'])} failed:", fg="red"))
        for what, why in res["failed"]:
            click.echo(click.style(f"  {what}: {why}", fg="red"))
        sys.exit(1)


@main.command("build-dashboard")
@click.option("--out", default=None, type=click.Path(dir_okay=False), help="Default: site/data.json.")
@click.option("--no-claims", is_flag=True, help="Do not read claims from the remote's study branches.")
@click.option("--no-mirror", is_flag=True, help="Do not read the gannet mirror ledgers.")
@click.option("--remote", default="origin", show_default=True)
def build_dashboard_cmd(out, no_claims, no_mirror, remote):
    """Write site/data.json for the GitHub Pages dashboard (site/index.html)."""
    import os
    import re

    from whale.dashboard import build, write
    from whale.fastq import load_ledger
    from whale.mirror import LEDGER_NAME, load_data_host

    repo = _repo()
    claims = None
    if not no_claims:
        try:
            claims = list_claims(repo, remote=remote)
        except ClaimError as exc:
            click.echo(click.style(f"claims unavailable: {exc}", fg="yellow"), err=True)

    mirrors = None
    if not no_mirror:
        try:
            base = load_data_host(repo).get("url_base")
        except FileNotFoundError:
            base = None
        if base:
            mirrors = {}
            for b in load_backlog(repo):
                url = f"{base.rstrip('/')}/raw/{b['accession']}/{LEDGER_NAME}"
                try:
                    mirrors[b["accession"]] = load_ledger(url)
                except (OSError, ValueError) as exc:
                    click.echo(click.style(f"mirror ledger unreadable for {b['accession']}: {exc}", fg="yellow"), err=True)
                    mirrors[b["accession"]] = None

    if os.environ.get("GITHUB_REPOSITORY"):
        repo_url = f"{os.environ.get('GITHUB_SERVER_URL', 'https://github.com')}/{os.environ['GITHUB_REPOSITORY']}"
    else:
        origin = subprocess.run(["git", "remote", "get-url", remote], cwd=repo, capture_output=True, text=True).stdout.strip()
        m = re.search(r"github\.com[:/](.+?)(\.git)?$", origin)
        repo_url = f"https://github.com/{m.group(1)}" if m else "https://github.com/RobertsLab/project-white-whale"

    data = build(repo, claims=claims, mirrors=mirrors, repo_url=repo_url)
    path = write(data, Path(out) if out else repo / "site" / "data.json")
    t = data["totals"]
    click.echo(f"Wrote {path}: {t['backlog']} studies, {t['registered']} registered, {t['in_progress']} in progress, "
               f"{t['done']} done; {t['mirrored_gb']} GB mirrored; {len(data['findings'])} finding(s)")


# --- claims ------------------------------------------------------------------

@main.command("start-study")
@click.option("--accession", required=True, help="BioProject accession, e.g. PRJNA682817.")
@click.option("--study-id", required=True, help="e.g. ROBERTS2022_PH_PLOIDY_WGBS (see registry/backlog.tsv).")
@click.option("--profile", required=True, help="Where it will run: desktop, raven or klone.")
@click.option("--part", default=None, help="Claim one part of a registered study, e.g. gill or 1of4.")
@click.option("--dir", "worktree_dir", default=None, type=click.Path(),
              help="Where to create the worktree (default: ../white-whale-<STUDY_ID>[-<PART>]).")
@click.option("--machine", default=None, help="Machine name recorded in the claim (default: host name).")
@click.option("--registered-by", default=None, help="Default: git user.name.")
@click.option("--base", default="main", show_default=True)
@click.option("--remote", default="origin", show_default=True)
@click.option("--no-push", is_flag=True, help="Create the branch and stubs locally only; nothing is claimed.")
@click.option("--no-pr", is_flag=True, help="Push the claim but do not open a draft PR.")
def start_study_cmd(accession, study_id, profile, part, worktree_dir, machine, registered_by,
                    base, remote, no_push, no_pr):
    """Check a study out to this machine: its own branch and worktree, claimed by pushing."""
    repo = _repo()
    try:
        res = start_study(accession, study_id, repo=repo, profile=profile, part=part,
                          worktree_dir=worktree_dir, registered_by=registered_by, machine=machine,
                          base=base, remote=remote, push=not no_push, open_pr=not no_pr)
    except ClaimError as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(f"Branch   : {res['branch']}")
    click.echo(f"Worktree : {res['worktree']}")
    for f in res["files"]:
        click.echo(f"  wrote  {f}")
    for link in res["linked"]:
        click.echo(f"  linked {link}")
    if res["pushed"]:
        click.echo(click.style("Claimed: branch pushed.", fg="green"))
        if res["pr_url"]:
            click.echo(f"Draft PR : {res['pr_url']}")
    else:
        click.echo(click.style("NOT claimed: --no-push given; other machines cannot see this study.",
                               fg="yellow"))
    click.echo(f"\nNext: cd {res['worktree']} and follow docs/adding_a_study.md.")


@main.command("set-status")
@click.option("--stage", required=True, help="A pipeline_stage.yaml term, e.g. fastq_staged.")
@click.option("--note", default=None)
@click.option("--run-name", "nextflow_run", default=None, help="Nextflow run name to record.")
@click.option("--profile", default=None, help="Record a move to another profile.")
@click.option("--machine", default=None, help="Default: host name.")
@click.option("--force", is_flag=True, help="Allow moving back a stage or reopening an excluded study.")
@click.option("--no-commit", is_flag=True)
@click.option("--no-push", is_flag=True)
@click.option("--remote", default="origin", show_default=True)
def set_status_cmd(stage, note, nextflow_run, profile, machine, force, no_commit, no_push, remote):
    """Record the stage the study on this branch has reached (commits and pushes by default)."""
    repo = _repo()
    try:
        new = set_status(stage, repo=repo, note=note, nextflow_run=nextflow_run, profile=profile,
                         machine=machine, force=force, commit=not no_commit,
                         push=not (no_push or no_commit), remote=remote)
    except ClaimError as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(f"{new['study_id']}{' part ' + new['part'] if new.get('part') else ''}: "
               f"{new['stage']} on {new['machine']} [{new['profile']}]")


@main.command("claims")
@click.option("--no-fetch", is_flag=True)
@click.option("--remote", default="origin", show_default=True)
def claims_cmd(no_fetch, remote):
    """Every study checked out to a machine, with its stage, read from the pushed branches."""
    repo = _repo()
    try:
        rows = list_claims(repo, remote=remote, fetch=not no_fetch)
    except ClaimError as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(_table(rows, ["study_id", "part", "accession", "stage", "profile", "machine",
                             "claimed_by", "claimed_on", "updated"]))


@main.command("release-study")
@click.option("--branch", required=True, help="study/<ACCESSION>/<STUDY_ID>[.part-<PART>]")
@click.option("--reason", required=True, help="Why; posted on the PR when it is closed.")
@click.option("--remote", default="origin", show_default=True)
@click.option("--yes", is_flag=True, help="Do not ask for confirmation.")
def release_study_cmd(branch, reason, remote, yes):
    """Give a claim up so another machine can take it: close the PR, delete the remote branch."""
    repo = _repo()
    if not yes:
        click.confirm(f"Delete {remote}/{branch} and close its PR? Local worktree is kept.", abort=True)
    try:
        release_study(branch, repo=repo, remote=remote, reason=reason)
    except ClaimError as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(f"Released {branch}. The local worktree and branch are untouched; remove them with\n"
               f"  git worktree list   # find it\n  git worktree remove <path> && git branch -D {branch}")


@main.command("check-study-scope")
@click.option("--branch", required=True)
@click.option("--base", default="origin/main", show_default=True)
def check_study_scope_cmd(branch, base):
    """Fail if a study branch changes files its study does not own (run by CI on study PRs)."""
    repo = _repo()
    parsed = parse_branch(branch)
    if parsed is None:
        click.echo(f"{branch} is not a study branch; nothing to check.")
        return
    _, study_id, part = parsed
    changed = subprocess.run(["git", "diff", "--name-only", f"{base}...HEAD"], cwd=repo,
                             capture_output=True, text=True, check=True).stdout.split()
    yaml_rel = f"registry/studies/{study_id}.yaml"
    source = base if part else "HEAD"
    proc = subprocess.run(["git", "show", f"{source}:{yaml_rel}"], cwd=repo, capture_output=True, text=True)
    study = yaml.safe_load(proc.stdout) if proc.returncode == 0 else None
    try:
        registered = registered_on(base, repo)
    except ClaimError:
        registered = {}
    errors, warnings = check_scope(branch, changed, study, registered)
    for w in warnings:
        click.echo(click.style(f"warning: {w}", fg="yellow"))
    for e in errors:
        click.echo(click.style(f"error: {e}", fg="red"))
    if errors:
        sys.exit(1)
    click.echo(f"ok: {len(changed)} changed file(s), all owned by {study_id}{' part ' + part if part else ''}.")


@main.command("profiles")
def profiles_cmd():
    """The compute profiles a study can be claimed on."""
    repo = _repo()
    rows = [{"profile": k, "executor": v["executor"], "label": v["label"], "suited for": v["suited_for"]}
            for k, v in load_vocab(repo, "compute_profile").items()]
    click.echo(_table(rows, ["profile", "executor", "label", "suited for"]))


if __name__ == "__main__":
    main()
