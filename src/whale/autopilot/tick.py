"""One autopilot tick: gather, decide one transition per study, report.

Step A has only a dry run: `tick(..., dry_run=True)` returns the decisions and
changes nothing. A tick without dry_run raises, and an AUTOPILOT_PAUSED file
in the repository root makes every tick a no-op.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from whale.autopilot import slurm
from whale.autopilot.decide import ACTIVE_STATES, classify, decide
from whale.autopilot.gather import fetch, gather
from whale.autopilot.machines import Slots, load_machines
from whale.autopilot.model import Decision, StudyState
from whale.autopilot.nflog import pid_alive
from whale.claims import machine_name

PAUSE_FILE = "AUTOPILOT_PAUSED"


class AutopilotError(RuntimeError):
    pass


@dataclass
class TickResult:
    here: str
    paused: str | None = None             # contents of AUTOPILOT_PAUSED, if present
    decisions: list[Decision] = field(default_factory=list)
    states: list[StudyState] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def paused(repo: Path) -> str | None:
    p = Path(repo) / PAUSE_FILE
    if not p.exists():
        return None
    return p.read_text().strip() or "(no reason given)"


def _counts(states: list[StudyState]) -> tuple[dict[str, int], dict[str, int]]:
    """Per machine: studies running now, and studies active (claimed until summarized)."""
    runs: dict[str, int] = {}
    active: dict[str, int] = {}
    for s in states:
        state, m = classify(s), (s.machine or "")
        if not m:
            continue
        if state in ACTIVE_STATES:
            active[m] = active.get(m, 0) + 1
        live = s.run is not None and s.run.status == "live"
        if state == "running" and (live or s.run is None):
            runs[m] = runs.get(m, 0) + 1
    return runs, active


def decide_all(states: list[StudyState], ctx: dict, machines: dict, here: str) -> list[Decision]:
    runs, active = _counts(states)
    run_slots, active_slots = Slots(machines, runs), Slots(machines, active)
    claimed: set[str] = set()
    return [decide(s, ctx, machines, here=here, runs=run_slots, active=active_slots, claimed_this_tick=claimed)
            for s in states]


def attach_queues(states: list[StudyState], machines: dict, here: str,
                  runner: slurm.Runner | None = None) -> list[str]:
    """Read each sbatch machine's queue once and give every study on it its jobs. Returns warnings."""
    warnings = []
    for m in machines.values():
        if m.launch != "sbatch":
            continue
        on_it = [s for s in states if (s.machine or "").lower() == m.name.lower() and s.study_id]
        if not on_it:
            continue
        q = slurm.read_queue(m.name, m.ssh_host, here=here, **({"runner": runner} if runner else {}))
        if not q.ok:
            warnings.append(f"could not read the {m.name} queue from {here}: {q.error}")
        by_study = slurm.assign([s.study_id for s in on_it], q)
        for s in on_it:
            s.slurm_checked, s.slurm_error = q.ok, q.error
            s.slurm_jobs = by_study.get(s.study_id, [])
    return warnings


def tick(repo: Path, *, dry_run: bool, remote: str = "origin", base: str = "main", do_fetch: bool = True,
         here: str | None = None, machines_path: Path | None = None,
         alive: Callable[[int], bool] = pid_alive, queue_runner: slurm.Runner | None = None) -> TickResult:
    if not dry_run:
        raise AutopilotError("step A of the autopilot only has a dry run; pass --dry-run "
                             "(docs/autopilot.md: steps B-F add live transitions).")
    repo = Path(repo)
    result = TickResult(here=here or machine_name())
    result.paused = paused(repo)
    if result.paused is not None:
        return result
    if do_fetch and (err := fetch(repo, remote, base)):
        result.warnings.append(f"could not fetch {remote}; using refs as last fetched: {err}")
    machines = load_machines(repo, machines_path)
    if not machines:
        result.warnings.append("config/autopilot.yaml has no machines; launch decisions ignore limits")
    states, ctx = gather(repo, remote=remote, base=base, alive=alive)
    result.warnings += attach_queues(states, machines, result.here, queue_runner)
    result.states = states
    result.decisions = decide_all(states, ctx, machines, result.here)
    return result
