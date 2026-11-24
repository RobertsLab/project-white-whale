"""Read a Slurm machine's queue, so a tick elsewhere can see runs there.

One read-only call per sbatch machine per tick: `squeue -u $USER`, run locally
when the tick is on that machine, otherwise over SSH with BatchMode (it never
prompts, so a missing Duo/ControlMaster session fails fast instead of hanging).
Jobs are attached to a study when the job name or the last component of its
working directory contains the study id (a head job runs in the study
checkout, white-whale-<STUDY_ID>), or when a `whale-head-<PREFIX>` job's prefix
names exactly one study on that machine (the hand-written DECICOMP2025 head
job is `whale-head-DECICOMP2025`). Nextflow task jobs (`nf-PROCESS_(tag)`) run
in the work directory, so they are not attributed to a study.
"""
from __future__ import annotations

import shlex
import subprocess
from collections.abc import Callable
from dataclasses import dataclass, field

FORMAT = "%i|%j|%P|%a|%T|%M|%Z"
HEAD_PREFIX = "whale-head-"


@dataclass
class SlurmJob:
    job_id: str
    name: str
    partition: str
    account: str
    state: str                        # PENDING, RUNNING, ...
    elapsed: str
    workdir: str

    @property
    def is_head(self) -> bool:
        """Anything but a Nextflow task job (nf-...), e.g. a staging or head job of the study."""
        return not self.name.startswith("nf-")

    def describe(self) -> str:
        return f"job {self.job_id} {self.name} {self.state} {self.elapsed} on {self.partition} ({self.account})"


@dataclass
class Queue:
    machine: str
    jobs: list[SlurmJob] = field(default_factory=list)
    error: str | None = None          # why the queue could not be read; jobs is then empty

    @property
    def ok(self) -> bool:
        return self.error is None


Runner = Callable[[list[str]], subprocess.CompletedProcess]


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60)


def parse(text: str) -> list[SlurmJob]:
    jobs = []
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) != 7 or not parts[0].strip():
            continue                  # login banners, conda warnings
        jobs.append(SlurmJob(*(p.strip() for p in parts)))
    return jobs


def read_queue(machine_name: str, ssh_host: str | None, *, here: str, runner: Runner = _run) -> Queue:
    squeue = ["squeue", "-u", "$USER", "-h", "-o", FORMAT]
    if machine_name.lower() == (here or "").lower():
        cmd = ["bash", "-lc", " ".join(shlex.quote(a) if a != "$USER" else a for a in squeue)]
    elif ssh_host:
        remote = " ".join(shlex.quote(a) if a != "$USER" else a for a in squeue)
        cmd = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", ssh_host, remote]
    else:
        return Queue(machine_name, error=f"no ssh_host for {machine_name} in config/autopilot.yaml")
    try:
        proc = runner(cmd)
    except (OSError, subprocess.TimeoutExpired) as e:
        return Queue(machine_name, error=f"{cmd[0]} failed: {e}")
    if proc.returncode != 0:
        err = (proc.stderr or "").strip().splitlines()
        return Queue(machine_name, error=(err[-1] if err else f"exit {proc.returncode}"))
    return Queue(machine_name, jobs=parse(proc.stdout))


def _prefix_matches(prefix: str, study_id: str) -> bool:
    return study_id == prefix or study_id.startswith(prefix + "_") or study_id.startswith(prefix + ".")


def assign(study_ids: list[str], queue: Queue) -> dict[str, list[SlurmJob]]:
    """study id -> its jobs. By the full study id in the job name or the last component of its
    working directory (the study checkout); otherwise a `whale-head-<PREFIX>` job whose prefix
    names exactly one of study_ids (OUC2024 would be ambiguous, so it matches nothing)."""
    out: dict[str, list[SlurmJob]] = {}
    for j in queue.jobs:
        leaf = j.workdir.rstrip("/").split("/")[-1]
        hits = [sid for sid in study_ids if sid in j.name or sid in leaf]
        if not hits and j.name.startswith(HEAD_PREFIX):
            prefix = j.name[len(HEAD_PREFIX):]
            hits = [sid for sid in study_ids if prefix and _prefix_matches(prefix, sid)]
        if len(hits) == 1:
            out.setdefault(hits[0], []).append(j)
    return out
