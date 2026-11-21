"""Read a study's latest Nextflow run from its local checkout.

Only what this machine can read: the `.nextflow.log` in the study's worktree,
`.nextflow/history`, `.nextflow.pid`, the `.exitcode` / `.command.err` of a
failed task's work directory when it is on a local disk, and any `*.out` /
`*.log` beside the log (nohup output, sbatch output, fetch-fastq logs). A run
on another machine is simply not seen; the decision then says so.
"""
from __future__ import annotations

import contextlib
import os
import re
from collections.abc import Callable
from pathlib import Path

from whale.autopilot.model import CachedTask, FailedTask, RunEvidence

_LAUNCH = re.compile(r"nextflow\.cli\.Launcher - \$> (.+)$", re.M)
_SESSION = re.compile(r"Session UUID: (\S+)")
_RUN_NAME = re.compile(r"Run name: (\S+)")
_WORKDIR = re.compile(r"Work-dir: (\S+)")
_CACHED = re.compile(r"\[([0-9a-f]{2}/[0-9a-f]{6})\] Cached process > (\S+)(?: \((.+?)\))?\s*$", re.M)
_ERROR_PROC = re.compile(r"Error executing process > '(\S+)(?: \((.+?)\))?'")
_EXIT_STATUS = re.compile(r"Command exit status:\s*\n\s*(-?\d+)")
_EXIT_INLINE = re.compile(r"terminated with an error exit status \((-?\d+)\)")
_TASK_WORKDIR = re.compile(r"^Work dir:\s*\n\s*(\S+)", re.M)
_COMMAND_ERR = re.compile(r"Command error:\s*\n(.*?)(?:\n\s*\nWork dir:|\Z)", re.S)
_ABORTED = re.compile(r"Session aborted|Execution aborted|ERROR ~|Error executing process")
_COMPLETE = re.compile(r"Execution complete -- Goodbye")
_STATS = re.compile(r"WorkflowStats\[succeededCount=(\d+); failedCount=(\d+); ignoredCount=(\d+)")

EXTRA_LOG_GLOBS = ["*.out", "nextflow*.log", "fetch-fastq*.log", "slurm-*.out"]
MAX_EXTRA_BYTES = 200_000


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def parse_log(text: str, source: str = ".nextflow.log") -> RunEvidence:
    """A RunEvidence from one .nextflow.log; status is ok, failed or unknown (maybe live)."""
    cmd = _LAUNCH.search(text)
    command = cmd.group(1).strip() if cmd else None
    ev = RunEvidence(source=source, status="unknown", command=command, text=text)
    if (m := _SESSION.search(text)):
        ev.session = m.group(1)
    if (m := _RUN_NAME.search(text)):
        ev.run_name = m.group(1)
    if (m := _WORKDIR.search(text)):
        ev.workdir = m.group(1)
    ev.resumed = bool(command and re.search(r"(^|\s)-resume(\s|$)", command))
    ev.cached = [CachedTask(h, p, t or None) for h, p, t in _CACHED.findall(text)]

    if (m := _ERROR_PROC.search(text)):
        block = text[m.start():]
        task = FailedTask(process=m.group(1), tag=m.group(2) or None)
        if (e := _EXIT_STATUS.search(block)) or (e := _EXIT_INLINE.search(block)):
            task.exit_status = int(e.group(1))
        if (w := _TASK_WORKDIR.search(block)):
            task.workdir = w.group(1)
        if (c := _COMMAND_ERR.search(block)):
            task.stderr = c.group(1).strip()[-4000:]
        ev.failed.append(task)
    elif (m := _EXIT_INLINE.search(text)) and _ABORTED.search(text):
        proc = re.search(r"Process `(\S+)(?: \((.+?)\))?` terminated", text)
        ev.failed.append(FailedTask(process=proc.group(1) if proc else "?",
                                    tag=(proc.group(2) if proc else None),
                                    exit_status=int(m.group(1))))

    complete = bool(_COMPLETE.search(text))
    stats = _STATS.findall(text)
    failed_count = int(stats[-1][1]) - int(stats[-1][2]) if stats else 0
    if ev.failed or (complete and _ABORTED.search(text)):
        ev.status = "failed"
    elif complete and failed_count <= 0:
        ev.status = "ok"
    elif complete:
        ev.status = "failed"
    return ev


def read_history(path: Path) -> list[dict]:
    """Rows of .nextflow/history: date, duration, run_name, status (OK, ERR, -), session, command."""
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="replace").splitlines():
        cols = line.split("\t")
        if len(cols) >= 7:
            rows.append({"date": cols[0], "duration": cols[1], "run_name": cols[2],
                         "status": cols[3], "session": cols[5], "command": cols[6]})
    return rows


def _task_files(task: FailedTask) -> None:
    """Fill a failed task's exit status and stderr from its work dir, when it is readable here."""
    if not task.workdir:
        return
    wd = Path(task.workdir)
    try:
        code = (wd / ".exitcode").read_text().strip()
        if task.exit_status is None and code.lstrip("-").isdigit():
            task.exit_status = int(code)
    except OSError:
        pass
    if not task.stderr:
        with contextlib.suppress(OSError):
            task.stderr = (wd / ".command.err").read_text(errors="replace")[-4000:]


def _extra_logs(checkout: Path, main_log: Path) -> str:
    """The newest few other logs in the checkout, for patterns a Nextflow log does not hold."""
    seen = []
    for pattern in EXTRA_LOG_GLOBS:
        seen += [p for p in checkout.glob(pattern) if p.is_file() and p != main_log]
    seen = sorted(set(seen), key=lambda p: p.stat().st_mtime, reverse=True)[:4]
    out = []
    for p in seen:
        data = p.read_bytes()[-MAX_EXTRA_BYTES:]
        out.append(f"==> {p.name} <==\n{data.decode(errors='replace')}")
    return "\n".join(out)


def read_run(checkout: Path, *, alive: Callable[[int], bool] = pid_alive) -> RunEvidence | None:
    """The latest run in a study checkout, or None if it has never run Nextflow here."""
    checkout = Path(checkout)
    log = checkout / ".nextflow.log"
    extra = _extra_logs(checkout, log) if checkout.is_dir() else ""
    if not log.exists():
        if extra:
            return RunEvidence(source=str(checkout), status="unknown", extra_text=extra)
        return None
    ev = parse_log(log.read_text(errors="replace"), source=str(log))
    ev.extra_text = extra
    for task in ev.failed:
        _task_files(task)

    history = read_history(checkout / ".nextflow" / "history")
    if history:
        last = history[-1]
        if ev.session:
            ev.previous_runs = sum(1 for r in history[:-1] if r["session"] == ev.session)
        if last["run_name"] == ev.run_name and last["status"] in ("OK", "ERR") and ev.status == "unknown":
            ev.status = "ok" if last["status"] == "OK" else "failed"

    if ev.status == "unknown":
        pid_file = checkout / ".nextflow.pid"
        try:
            pid = int(pid_file.read_text().strip())
        except (OSError, ValueError):
            pid = None
        if pid is not None and alive(pid):
            ev.status = "live"
        elif pid is not None:
            # The head process is gone without writing an end: killed, power cut, reboot.
            ev.status = "stopped"
    return ev
