"""Check studies out to a machine: one branch per study (or study part).

Several computers (lab desktops, raven, klone) work through the backlog at
once. Each job is one study on its own branch,

    study/<ACCESSION>/<STUDY_ID>                 e.g. study/PRJNA682817/ROBERTS2022_PH_PLOIDY_WGBS
    study/<ACCESSION>/<STUDY_ID>.part-<PART>     e.g. study/PRJEB60400/POMS2023_ADAPTATION_WGBS.part-gill

opened as a draft pull request against main. Pushing the branch is the claim:
two machines racing for one study push different commits to the same branch
name and git rejects the second, so the loser stops before downloading
anything. Before pushing, `start_study` also refuses a study that is already
registered on main or already claimed by a pushed `study/` branch.

Each branch carries data/studies/<STUDY_ID>/status.yaml (or, for a part,
data/studies/<STUDY_ID>/parts/<PART>/status.yaml): which machine and profile
hold the study and which pipeline stage it has reached. `list_claims` reads
those files straight from the remote branches, so the ledger needs no server.

A part claim is for a study too large for one machine. The study is first
claimed whole and its registration merged; parts are then claimed against the
registered study and own only their own `parts/<PART>/` directory.

Adapted from AREE's intake/study_claims.py (github.com/sr320/AREE).
"""
from __future__ import annotations

import getpass
import re
import shutil
import socket
import subprocess
from datetime import date
from pathlib import Path

import yaml

from whale.paths import status_path
from whale.vocab import load_vocab, stage_order

BRANCH_PREFIX = "study/"
ACCESSION_RE = re.compile(r"^PRJ[NED][AB][0-9]+$")
# Same pattern as study_id in schemas/study.schema.json.
STUDY_ID_RE = re.compile(r"^[A-Z][A-Z0-9]*[0-9]{4}(_[A-Z0-9]+)+$")
PART_RE = re.compile(r"^[A-Za-z0-9-]+$")
_BRANCH_RE = re.compile(
    r"^study/(?P<accession>[^/]+)/(?P<study_id>[A-Za-z0-9_-]+)(?:\.part-(?P<part>[A-Za-z0-9-]+))?$"
)
TERMINAL_STAGES = {"excluded", "released"}
STUDY_PR_TEMPLATE = Path(".github/PULL_REQUEST_TEMPLATE/study.md")

# Gitignored symlinks pointing a checkout at this machine's data drives. A new
# worktree does not get them from git, so they are recreated in it.
DRIVE_LINK_GLOBS = ["data/raw", "data/methylomes", "data/reference/GCF_*"]


class ClaimError(RuntimeError):
    pass


def machine_name() -> str:
    return socket.gethostname().split(".")[0]


def default_user(repo: Path) -> str:
    proc = subprocess.run(["git", "config", "user.name"], cwd=repo, capture_output=True, text=True)
    return proc.stdout.strip() or getpass.getuser()


def branch_name(accession: str, study_id: str, part: str | None = None) -> str:
    if not ACCESSION_RE.match(accession):
        raise ClaimError(f"accession {accession!r} must be a BioProject accession (PRJNA…, PRJEB…, PRJDB…).")
    if not STUDY_ID_RE.match(study_id):
        raise ClaimError(f"study_id {study_id!r} must look like FIRSTAUTHOR<YEAR>_<TOPIC>[_<CHEMISTRY>] "
                         f"({STUDY_ID_RE.pattern}).")
    if part is not None and not PART_RE.match(part):
        raise ClaimError(f"part {part!r} must match {PART_RE.pattern}.")
    suffix = f".part-{part}" if part else ""
    return f"{BRANCH_PREFIX}{accession}/{study_id}{suffix}"


def parse_branch(name: str) -> tuple[str, str, str | None] | None:
    """(accession, study_id, part) for a study branch, or None for any other branch."""
    m = _BRANCH_RE.match(name)
    return (m["accession"], m["study_id"], m["part"]) if m else None


def accessions_in(study: dict) -> set[str]:
    """Every accession-like token in a study's `accessions` block (values may be free text)."""
    out = set()
    for value in (study.get("accessions") or {}).values():
        if value:
            out.add(str(value))
            out.update(t for t in re.split(r"[^A-Za-z0-9._]+", str(value)) if t)
    return out


def find_conflicts(accession: str, study_id: str, part: str | None,
                   registered: dict[str, set[str]], claimed_branches: list[str]) -> list[str]:
    """Reasons this claim cannot be made. Empty means it is free.

    `registered` maps each study_id on main to its accessions; `claimed_branches`
    are study branch names pushed (or created locally).
    """
    problems = []
    if part is None:
        if study_id in registered:
            problems.append(f"study_id {study_id} is already registered on main.")
        for sid, accs in sorted(registered.items()):
            if sid != study_id and accession in accs:
                problems.append(f"accession {accession} is already registered on main as {sid}.")
        for branch in sorted(claimed_branches):
            parsed = parse_branch(branch)
            if parsed and (parsed[0] == accession or parsed[1] == study_id):
                problems.append(f"already claimed by branch {branch}.")
        return problems

    if study_id not in registered:
        problems.append(f"a part can only be claimed for a study registered on main; {study_id} is not. "
                        "Claim the study whole first and merge its registration.")
    elif accession not in registered[study_id]:
        problems.append(f"{accession} is not an accession of {study_id} on main.")
    for branch in sorted(claimed_branches):
        parsed = parse_branch(branch)
        if not parsed or parsed[1] != study_id:
            continue
        if parsed[2] is None:
            problems.append(f"{study_id} is still claimed whole by branch {branch}; merge or release it first.")
        elif parsed[2] == part:
            problems.append(f"part {part} is already claimed by branch {branch}.")
    return problems


# --- files a claim writes ----------------------------------------------------

def _sub_once(text: str, pattern: str, repl: str, what: str) -> str:
    new, n = re.subn(pattern, repl, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise ClaimError(f"registry/studies/_TEMPLATE.yaml has no {what} line to fill; "
                         "update whale/claims.py to match the template.")
    return new


def render_study_yaml(template: str, *, study_id: str, accession: str, registered_by: str,
                      today: str) -> str:
    """The template with identity and provenance filled in; the rest is the curator's."""
    text = _sub_once(template, r"^study_id: .*$", f"study_id: {study_id}", "study_id")
    text = _sub_once(text, r"^(  bioproject:) null\s*$", rf'\1 "{accession}"', "accessions.bioproject")
    text = _sub_once(text, r"^(  registered_by:) .*$", rf'\1 "{registered_by}"', "registered_by")
    text = _sub_once(text, r"^(  date_registered:) .*$", rf'\1 "{today}"', "date_registered")
    return text


def config_name(study_id: str, part: str | None = None) -> str:
    return f"config/{study_id}.part-{part}.config" if part else f"config/{study_id}.config"


def render_config(study_id: str, accession: str, profile: str, part: str | None = None) -> str:
    part_line = f"\n        part         = '{part}'" if part else ""
    sheet = (f"data/studies/{study_id}/parts/{part}/samplesheet.tsv" if part
             else f"data/studies/{study_id}/samplesheet.tsv")
    return f"""/*
 * White Whale :: run configuration for {study_id}{f' part {part}' if part else ''} ({accession})
 *
 * Claimed with `whale start-study` for profile `{profile}`. Run from the
 * repository root of this study's worktree:
 *
 *   nextflow run workflows/methylation -c {config_name(study_id, part)} -profile {profile}{',docker' if profile == 'desktop' else ''}
 *
 * Machine-specific paths come from WHALE_RAW, WHALE_REFERENCE, WHALE_METHYLOMES
 * and WHALE_WORKDIR (docs/compute_profiles.md); never write one here.
 */

params {{
    study_id     = '{study_id}'{part_line}
    chemistry    = null   // assay_chemistry.yaml term; must match the study YAML
    samplesheet  = '{sheet}'
}}

// Study-specific resource overrides go in a `process {{ withName: ... }}` block.
// Klone runs: also pass `-c` with a run-only config last, as described in
// docs/compute_profiles.md, so workstation time limits do not override the profile.
"""


def new_status(*, study_id: str, accession: str, part: str | None, profile: str, machine: str,
               claimed_by: str, today: str) -> dict:
    return {
        "study_id": study_id,
        "accession": accession,
        "part": part,
        "stage": "claimed",
        "profile": profile,
        "machine": machine,
        "claimed_by": claimed_by,
        "claimed_on": today,
        "nextflow_run": None,
        "notes": None,
        "history": [{"stage": "claimed", "date": today, "machine": machine, "note": None}],
    }


def advance_status(status: dict, stage: str, *, repo: Path, machine: str, today: str,
                   note: str | None = None, nextflow_run: str | None = None,
                   profile: str | None = None, force: bool = False) -> dict:
    """A copy of `status` moved to `stage`. Stages only move forward unless `force`."""
    order = stage_order(repo)
    if stage not in order:
        raise ClaimError(f"unknown stage {stage!r}; one of: {', '.join(sorted(order, key=order.get))}")
    current = status["stage"]
    if current in TERMINAL_STAGES and not force:
        raise ClaimError(f"study is {current}; pass --force to reopen it.")
    if stage not in TERMINAL_STAGES and order[stage] < order[current] and not force:
        raise ClaimError(f"{stage} comes before the current stage {current}; pass --force to go back.")
    if profile is not None and profile not in load_vocab(repo, "compute_profile"):
        raise ClaimError(f"unknown profile {profile!r}")
    new = dict(status)
    new["stage"] = stage
    new["machine"] = machine
    if profile:
        new["profile"] = profile
    if nextflow_run:
        new["nextflow_run"] = nextflow_run
    new["history"] = [*status["history"], {"stage": stage, "date": today, "machine": machine, "note": note}]
    return new


def dump_status(status: dict) -> str:
    header = ("# Processing status of this study on its claim branch. Update with\n"
              "# `whale set-status --stage <stage>`; read by `whale claims`.\n")
    return header + yaml.safe_dump(status, sort_keys=False)


# --- git / GitHub ------------------------------------------------------------

def _run(args: list[str], cwd) -> str:
    proc = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ClaimError(f"`{' '.join(args)}` failed:\n{proc.stderr.strip() or proc.stdout.strip()}")
    return proc.stdout


def registered_on(ref: str, repo) -> dict[str, set[str]]:
    """study_id -> accessions for every study YAML at `ref` (e.g. origin/main)."""
    names = _run(["git", "ls-tree", "--name-only", ref, "registry/studies/"], repo).split()
    out = {}
    for name in names:
        if not name.endswith(".yaml") or Path(name).name.startswith("_"):
            continue
        study = yaml.safe_load(_run(["git", "show", f"{ref}:{name}"], repo)) or {}
        out[study.get("study_id") or Path(name).stem] = accessions_in(study)
    return out


def remote_study_branches(remote: str, repo) -> list[str]:
    out = _run(["git", "ls-remote", "--heads", remote, f"refs/heads/{BRANCH_PREFIX}*"], repo)
    return [line.split("refs/heads/", 1)[1] for line in out.splitlines() if "refs/heads/" in line]


def local_study_branches(repo) -> list[str]:
    """Study branches on this machine, including unpushed (`--no-push`) ones."""
    out = _run(["git", "branch", "--list", f"{BRANCH_PREFIX}*", "--format=%(refname:short)"], repo)
    return [b.strip() for b in out.splitlines() if b.strip()]


def current_branch(repo) -> str:
    return _run(["git", "rev-parse", "--abbrev-ref", "HEAD"], repo).strip()


def link_drives(repo: Path, worktree: Path) -> list[str]:
    """Recreate the repo's data-drive symlinks in a new worktree. Returns the paths linked."""
    linked = []
    for pattern in DRIVE_LINK_GLOBS:
        for src in sorted(repo.glob(pattern)):
            if not src.is_symlink():
                continue
            rel = src.relative_to(repo)
            dest = worktree / rel
            if dest.exists() or dest.is_symlink():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.symlink_to(src.readlink())
            linked.append(str(rel))
    return linked


def start_study(accession: str, study_id: str, *, repo, profile: str, part: str | None = None,
                worktree_dir=None, registered_by: str | None = None, machine: str | None = None,
                base: str = "main", remote: str = "origin", push: bool = True,
                open_pr: bool = True) -> dict:
    """Create the study branch in its own worktree, commit the stubs, and claim it.

    Returns {branch, worktree, files, linked, pushed, pr_url}.
    """
    repo = Path(repo)
    if profile not in load_vocab(repo, "compute_profile"):
        raise ClaimError(f"unknown profile {profile!r}; see registry/controlled_vocabularies/compute_profile.yaml")
    branch = branch_name(accession, study_id, part)
    suffix = f"-{part}" if part else ""
    worktree = Path(worktree_dir) if worktree_dir else repo.parent / f"white-whale-{study_id}{suffix}"
    if worktree.exists():
        raise ClaimError(f"{worktree} already exists; pass --dir to put the worktree elsewhere.")
    registered_by = registered_by or default_user(repo)
    machine = machine or machine_name()

    _run(["git", "fetch", "--quiet", remote, base], repo)
    base_ref = f"{remote}/{base}"
    claimed = sorted(set(remote_study_branches(remote, repo)) | set(local_study_branches(repo)))
    problems = find_conflicts(accession, study_id, part, registered_on(base_ref, repo), claimed)
    if problems:
        raise ClaimError("Cannot claim this study:\n  - " + "\n  - ".join(problems))

    _run(["git", "worktree", "add", "--quiet", "-b", branch, str(worktree), base_ref], repo)

    today = date.today().isoformat()
    status = new_status(study_id=study_id, accession=accession, part=part, profile=profile,
                        machine=machine, claimed_by=registered_by, today=today)
    files = {
        str(status_path(Path("."), study_id, part)): dump_status(status),
        config_name(study_id, part): render_config(study_id, accession, profile, part),
    }
    if part is None:
        template = (worktree / "registry/studies/_TEMPLATE.yaml").read_text()
        files[f"registry/studies/{study_id}.yaml"] = render_study_yaml(
            template, study_id=study_id, accession=accession, registered_by=registered_by, today=today)
    for rel, text in files.items():
        (worktree / rel).parent.mkdir(parents=True, exist_ok=True)
        (worktree / rel).write_text(text)
    _run(["git", "add", *files], worktree)
    what = f"{study_id} part {part}" if part else study_id
    _run(["git", "commit", "--quiet", "-m", f"Claim {what} ({accession}) on {machine} [{profile}]"], worktree)

    result = {"branch": branch, "worktree": str(worktree), "files": sorted(files),
              "linked": link_drives(repo, worktree), "pushed": False, "pr_url": None}
    if not push:
        return result
    try:
        _run(["git", "push", "--quiet", "-u", remote, branch], worktree)
    except ClaimError as exc:
        raise ClaimError(
            f"{exc}\n\nThe push was rejected, most likely because another machine claimed "
            f"{branch} first. Nothing was downloaded. Remove the local claim with:\n"
            f"  git worktree remove {worktree} && git branch -D {branch}"
        ) from None
    result["pushed"] = True

    if open_pr and shutil.which("gh"):
        body = (f"Claimed by `whale start-study` on **{machine}** (profile `{profile}`), {today}. "
                f"Draft until the study reaches `summarized`.\n\n"
                + (worktree / STUDY_PR_TEMPLATE).read_text())
        title = f"[study] {study_id}{f' part {part}' if part else ''} ({accession})"
        result["pr_url"] = _run(
            ["gh", "pr", "create", "--draft", "--base", base, "--head", branch,
             "--title", title, "--body", body], worktree).strip()
    return result


def set_status(stage: str, *, repo, note: str | None = None, nextflow_run: str | None = None,
               profile: str | None = None, machine: str | None = None, force: bool = False,
               commit: bool = True, push: bool = True, remote: str = "origin") -> dict:
    """Move the study on the current branch to `stage`; commit and push so the ledger sees it."""
    repo = Path(repo)
    branch = current_branch(repo)
    parsed = parse_branch(branch)
    if parsed is None:
        raise ClaimError(f"{branch!r} is not a study branch; run this inside the study's worktree.")
    _, study_id, part = parsed
    path = status_path(repo, study_id, part)
    if not path.exists():
        raise ClaimError(f"{path.relative_to(repo)} is missing on this branch.")
    status = yaml.safe_load(path.read_text())
    new = advance_status(status, stage, repo=repo, machine=machine or machine_name(),
                         today=date.today().isoformat(), note=note, nextflow_run=nextflow_run,
                         profile=profile, force=force)
    path.write_text(dump_status(new))
    if commit:
        rel = str(path.relative_to(repo))
        _run(["git", "add", rel], repo)
        _run(["git", "commit", "--quiet", "-m", f"{study_id}{f' part {part}' if part else ''}: {stage}"], repo)
        if push:
            _run(["git", "push", "--quiet", remote, branch], repo)
    return new


def list_claims(repo, *, remote: str = "origin", fetch: bool = True) -> list[dict]:
    """One row per pushed study branch, from its status file."""
    repo = Path(repo)
    if fetch:
        _run(["git", "fetch", "--quiet", "--prune", remote,
              f"+refs/heads/{BRANCH_PREFIX}*:refs/remotes/{remote}/{BRANCH_PREFIX}*"], repo)
    rows = []
    for branch in sorted(remote_study_branches(remote, repo)):
        parsed = parse_branch(branch)
        if parsed is None:
            continue
        accession, study_id, part = parsed
        rel = status_path(Path("."), study_id, part)
        proc = subprocess.run(["git", "show", f"{remote}/{branch}:{rel}"], cwd=repo,
                              capture_output=True, text=True)
        status = (yaml.safe_load(proc.stdout) or {}) if proc.returncode == 0 else {}
        last = (status.get("history") or [{}])[-1]
        rows.append({
            "study_id": study_id,
            "part": part or "",
            "accession": accession,
            "stage": status.get("stage", "unknown"),
            "profile": status.get("profile", ""),
            "machine": status.get("machine", ""),
            "claimed_by": status.get("claimed_by", ""),
            "claimed_on": status.get("claimed_on", ""),
            "updated": last.get("date", ""),
            "branch": branch,
        })
    return rows


def release_study(branch: str, *, repo, remote: str = "origin", reason: str,
                  close_pr: bool = True) -> None:
    """Give up a claim: close its PR with the reason and delete the remote branch.

    The local worktree and branch are left alone (they may hold downloaded data);
    remove them yourself once nothing in them is needed.
    """
    if parse_branch(branch) is None:
        raise ClaimError(f"{branch!r} is not a study branch.")
    repo = Path(repo)
    if close_pr and shutil.which("gh"):
        subprocess.run(["gh", "pr", "close", branch, "--comment", f"Released: {reason}"],
                       cwd=repo, capture_output=True, text=True)
    _run(["git", "push", "--quiet", remote, "--delete", branch], repo)
