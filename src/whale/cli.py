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
