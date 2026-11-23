"""Per-machine limits (config/autopilot.yaml) and the slot bookkeeping of one tick."""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path

import yaml

CONFIG = Path("config/autopilot.yaml")


@dataclass
class Machine:
    name: str
    profile: str
    max_concurrent_runs: int = 1
    disk_floor_gb: float = 0
    disk_path: str | None = None
    launch: str = "local"
    gannet_key: bool = False
    gh: bool = True
    env: dict = field(default_factory=dict)
    extra_configs: list = field(default_factory=list)
    limits_config: str | None = None
    slurm_accounts: list = field(default_factory=list)
    head_job: str | None = None
    ssh_host: str | None = None       # sbatch machines: where a tick elsewhere reads `squeue`
    checkout_root: str | None = None  # holds the study worktrees, white-whale-<STUDY_ID>[-<PART>]
    path_prepend: list = field(default_factory=list)  # put first on PATH in a job (whale, nextflow, java)
    setup: str | None = None          # shell run first in a job, e.g. "unset TMPDIR"
    log_dir: str | None = None        # sbatch machines: where a job's --output goes
    snp_mask_root: str | None = None  # C/T SNP mask folders, <root>/<MASK_ID>; may contain {assembly}; unset: none

    def exports(self, extra: dict | None = None) -> str:
        return " ".join(f"{k}={v}" for k, v in {**self.env, **(extra or {})}.items())

    def snp_mask(self, mask_id: str | None, assembly: str | None = None) -> str | None:
        """The mask folder to pass as WHALE_SNP_MASK, or None when this machine has no masks."""
        if not mask_id or not self.snp_mask_root:
            return None
        root = self.snp_mask_root
        if "{assembly}" in root:
            if not assembly:
                return None
            root = root.replace("{assembly}", assembly)
        return f"{root.rstrip('/')}/{mask_id}"

    def checkout(self, study_id: str, part: str | None = None) -> str | None:
        """The study's worktree on this machine (as `whale start-study` names it), if checkout_root is set."""
        if not self.checkout_root:
            return None
        return f"{self.checkout_root.rstrip('/')}/white-whale-{study_id}{f'-{part}' if part else ''}"

    def sbatch(self, cmd: str, *, study_id: str, part: str | None = None, kind: str = "head",
               env: dict | None = None) -> str:
        """`sbatch ... --wrap` that runs cmd in the study's worktree with this machine's PATH and env.

        The job is named whale-<kind>-<STUDY_ID>[-<PART>] and runs with --chdir in the worktree,
        so `whale autopilot` reading the queue (slurm.assign) attributes it to the study.
        """
        name = f"whale-{kind}-{study_id}{f'-{part}' if part else ''}"
        opts = [self.head_job or "", f"--job-name={name}"]
        where = self.checkout(study_id, part)
        if where:
            opts.append(f"--chdir={where}")
        if self.log_dir:
            opts.append(f"--output={self.log_dir.rstrip('/')}/{name}-%j.log")
        steps = []
        if self.path_prepend:
            steps.append("export PATH=" + ":".join(self.path_prepend) + ":$PATH")
        if self.setup:
            steps.append(self.setup)
        if self.env or env:
            steps.append(f"export {self.exports(env)}")
        steps.append(cmd)
        return "sbatch " + " ".join(o for o in opts if o) + " --wrap '" + "; ".join(steps) + "'"


def load_machines(repo: Path, path: Path | None = None) -> dict[str, Machine]:
    p = path or (Path(repo) / CONFIG)
    if not p.exists():
        return {}
    doc = yaml.safe_load(p.read_text()) or {}
    out = {}
    for name, cfg in (doc.get("machines") or {}).items():
        known = {k: v for k, v in (cfg or {}).items() if k in Machine.__dataclass_fields__}
        out[name] = Machine(name=name, **known)
    return out


def find_machine(machines: dict[str, Machine], name: str | None) -> Machine | None:
    if not name:
        return None
    for m in machines.values():
        if m.name.lower() == name.lower():
            return m
    return None


def machines_for_profile(machines: dict[str, Machine], profile: str | None) -> list[Machine]:
    return [m for m in machines.values() if m.profile == profile]


def free_gb(machine: Machine, here: str) -> float | None:
    """Free space on the machine's disk_path, or None when the tick is not running there."""
    if not machine.disk_path or machine.name.lower() != (here or "").lower():
        return None
    try:
        return shutil.disk_usage(machine.disk_path).free / 1e9
    except OSError:
        return None


class Slots:
    """Concurrent-run slots per machine, filled as one tick's decisions are made."""

    def __init__(self, machines: dict[str, Machine], used: dict[str, int] | None = None):
        self.machines = machines
        self.used = {k.lower(): v for k, v in (used or {}).items()}

    def count(self, name: str) -> int:
        return self.used.get(name.lower(), 0)

    def has_room(self, name: str) -> bool:
        m = find_machine(self.machines, name)
        return m is None or self.count(name) < m.max_concurrent_runs

    def take(self, name: str) -> None:
        self.used[name.lower()] = self.count(name) + 1
