"""Read the state of every study, without changing anything.

Sources, all read-only: registry/backlog.tsv and registry/studies/ on main;
every pushed study branch (its status.yaml, config, study YAML, qc_summary.tsv,
study_summary.yaml, run_manifest.json); whether the branch is merged into
main; and, for a branch checked out in a worktree on this machine, its latest
Nextflow run (whale.autopilot.nflog). The only network call is an optional
`git fetch` of main and the study branches, which changes nothing but this
checkout's remote-tracking refs.
"""
from __future__ import annotations

import csv
import io
import json
import subprocess
from collections.abc import Callable
from pathlib import Path

import yaml

from whale.autopilot.model import StudyState
from whale.autopilot.nflog import pid_alive, read_run
from whale.autopilot.playbook import config_string
from whale.claims import BRANCH_PREFIX, config_name, parse_branch


class Git:
    def __init__(self, repo: Path):
        self.repo = Path(repo)

    def run(self, *args: str, check: bool = False) -> subprocess.CompletedProcess:
        return subprocess.run(["git", *args], cwd=self.repo, capture_output=True, text=True, check=check)

    def show(self, ref: str, path: str) -> str | None:
        proc = self.run("show", f"{ref}:{path}")
        return proc.stdout if proc.returncode == 0 else None

    def ls(self, ref: str, path: str) -> list[str]:
        proc = self.run("ls-tree", "--name-only", ref, path.rstrip("/") + "/")
        return proc.stdout.split() if proc.returncode == 0 else []

    def has_ref(self, ref: str) -> bool:
        return self.run("rev-parse", "--verify", "--quiet", ref).returncode == 0

    def is_ancestor(self, a: str, b: str) -> bool:
        return self.run("merge-base", "--is-ancestor", a, b).returncode == 0

    def remote_branches(self, remote: str) -> list[str]:
        proc = self.run("for-each-ref", "--format=%(refname)", f"refs/remotes/{remote}/{BRANCH_PREFIX}")
        prefix = f"refs/remotes/{remote}/"
        return sorted(line[len(prefix):] for line in proc.stdout.splitlines() if line.startswith(prefix))

    def worktrees(self) -> dict[str, str]:
        """branch -> worktree path, for every worktree of this repository."""
        out, path = {}, None
        for line in self.run("worktree", "list", "--porcelain").stdout.splitlines():
            if line.startswith("worktree "):
                path = line.split(" ", 1)[1]
            elif line.startswith("branch refs/heads/") and path:
                out[line[len("branch refs/heads/"):]] = path
        return out


def fetch(repo: Path, remote: str = "origin", base: str = "main") -> str | None:
    """Update remote-tracking refs for main and the study branches. Returns an error, or None."""
    proc = Git(repo).run("fetch", "--quiet", "--prune", remote, f"+refs/heads/{base}:refs/remotes/{remote}/{base}",
                         f"+refs/heads/{BRANCH_PREFIX}*:refs/remotes/{remote}/{BRANCH_PREFIX}*")
    return None if proc.returncode == 0 else (proc.stderr.strip() or "git fetch failed")


def _yaml(text: str | None) -> dict | None:
    if text is None:
        return None
    try:
        doc = yaml.safe_load(text)
    except yaml.YAMLError:
        return None
    return doc if isinstance(doc, dict) else None


def _tsv(text: str | None) -> list[dict] | None:
    if text is None:
        return None
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def _json(text: str | None) -> dict | None:
    if text is None:
        return None
    try:
        return json.loads(text)
    except ValueError:
        return None


def _data_dir(study_id: str, part: str | None) -> str:
    base = f"data/studies/{study_id}"
    return f"{base}/parts/{part}" if part else base


def _fill_from_ref(git: Git, ref: str, s: StudyState) -> None:
    sid, part = s.study_id, s.part
    d = _data_dir(sid, part)
    s.status = _yaml(git.show(ref, f"{d}/status.yaml"))
    s.study_yaml = _yaml(git.show(ref, f"registry/studies/{sid}.yaml"))
    s.config_text = git.show(ref, config_name(sid, part))
    run_cfgs = [p for p in git.ls(ref, d) if p.endswith(".run.config")]
    s.run_config_text = git.show(ref, run_cfgs[0]) if run_cfgs else None
    s.qc = _tsv(git.show(ref, f"{d}/qc_summary.tsv"))
    s.summary = _yaml(git.show(ref, f"{d}/study_summary.yaml"))
    s.run_manifest = _json(git.show(ref, f"{d}/run_manifest.json"))


def newest_masks(doc: dict | None) -> dict[str, str]:
    """assembly -> the newest mask id in data/reference/snp_masks.yaml (masks are appended, never replaced)."""
    out: dict[str, str] = {}
    for m in (doc or {}).get("masks") or []:
        if m.get("genome_assembly") and m.get("mask_id"):
            out[m["genome_assembly"]] = m["mask_id"]
    return out


def attach_masks(states, masks: dict[str, str], default_assembly: str | None) -> None:
    """Each study's genome assembly (its config, else base.config's) and the SNP mask for it."""
    for s in states:
        s.genome_assembly = config_string(s.config_text, "genome_assembly") or default_assembly
        s.snp_mask_id = masks.get(s.genome_assembly) if s.genome_assembly else None


def gather(repo: Path, *, remote: str = "origin", base: str = "main",
           alive: Callable[[int], bool] = pid_alive) -> tuple[list[StudyState], dict]:
    """Every study's state, plus context: {main_ref, main_assembly, workflows}."""
    git = Git(repo)
    main_ref = f"{remote}/{base}" if git.has_ref(f"{remote}/{base}") else base

    backlog = _tsv(git.show(main_ref, "registry/backlog.tsv")) or []
    registered = {}
    for name in git.ls(main_ref, "registry/studies"):
        if name.endswith(".yaml") and not Path(name).name.startswith("_"):
            doc = _yaml(git.show(main_ref, name)) or {}
            registered[doc.get("study_id") or Path(name).stem] = doc
    reg_by_acc = {(d.get("accessions") or {}).get("bioproject"): sid for sid, d in registered.items()}

    worktrees = git.worktrees()
    states: dict[str, StudyState] = {}
    for branch in git.remote_branches(remote):
        parsed = parse_branch(branch)
        if parsed is None:
            continue
        acc, sid, part = parsed
        key = f"{sid}.part-{part}" if part else sid
        s = StudyState(key=key, accession=acc, study_id=sid, part=part, branch=branch,
                       registered_on_main=sid in registered)
        ref = f"{remote}/{branch}"
        _fill_from_ref(git, ref, s)
        s.merged = git.is_ancestor(ref, main_ref)
        wt = worktrees.get(branch)
        if wt and Path(wt).is_dir():
            s.worktree = wt
            s.run = read_run(Path(wt), alive=alive)
        states[key] = s

    # Registered on main with no branch left (merged and the branch deleted, or registered for parts).
    for sid, doc in registered.items():
        if sid in states or any(s.study_id == sid for s in states.values()):
            continue
        acc = (doc.get("accessions") or {}).get("bioproject") or ""
        s = StudyState(key=sid, accession=acc, study_id=sid, registered_on_main=True)
        _fill_from_ref(git, main_ref, s)
        s.merged = s.summary is not None or (s.stage in ("summarized", "merged"))
        states[sid] = s

    claimed_accs = {s.accession for s in states.values()}
    for row in backlog:
        acc = row.get("accession")
        held = [s for s in states.values() if s.accession == acc]
        for s in held:
            s.backlog = row
        if held or acc in claimed_accs:
            continue
        suggested = row.get("suggested_study_id") or ""
        named = suggested if suggested and not suggested.startswith("(") else None
        key = named or acc
        states[key] = StudyState(key=key, accession=acc, study_id=named, backlog=row,
                                 registered_on_main=acc in reg_by_acc)

    assemblies = [((s.summary or {}).get("genome_assembly")) for s in states.values() if s.merged]
    assemblies = [a for a in assemblies if a]
    main_assembly = max(set(assemblies), key=assemblies.count) if assemblies else None
    workflows = set(git.ls(main_ref, "workflows"))
    masks = newest_masks(_yaml(git.show(main_ref, "data/reference/snp_masks.yaml")))
    default_assembly = config_string(git.show(main_ref, "config/base.config"), "genome_assembly")
    attach_masks(states.values(), masks, default_assembly)
    context = {"main_ref": main_ref, "main_assembly": main_assembly,
               "workflows": {Path(w).name for w in workflows}, "snp_masks": masks}

    order = {row.get("accession"): i for i, row in enumerate(backlog)}
    ordered = sorted(states.values(), key=lambda s: (order.get(s.accession, len(order)), s.key))
    return ordered, context
