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

    def exports(self) -> str:
        return " ".join(f"{k}={v}" for k, v in self.env.items())


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
