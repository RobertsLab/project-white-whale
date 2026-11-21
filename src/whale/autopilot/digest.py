"""Render a tick: one line per study, and a Markdown digest for people."""
from __future__ import annotations

from datetime import date

from whale.autopilot.model import Decision
from whale.autopilot.tick import TickResult

SECTIONS = [
    ("advance", "Would advance"),
    ("running", "Running"),
    ("blocked", "Blocked: needs a person"),
    ("waiting", "Waiting"),
    ("done", "Done"),
]


def _arrow(d: Decision) -> str:
    return f"{d.state} → {d.next_state}" if d.next_state and d.next_state != d.state else d.state


def lines(result: TickResult) -> list[str]:
    """One line per study: key, transition, kind, machine, reason."""
    if not result.decisions:
        return []
    w_key = max(len(d.key) for d in result.decisions)
    w_arrow = max(len(_arrow(d)) for d in result.decisions)
    out = []
    for d in result.decisions:
        tag = f"[{d.playbook}] " if d.playbook else ""
        out.append(f"{d.key.ljust(w_key)}  {_arrow(d).ljust(w_arrow)}  {d.kind:<7}  "
                   f"{(d.machine or '-'):<6}  {tag}{d.reason}")
    return out


def markdown(result: TickResult, today: str | None = None) -> str:
    today = today or date.today().isoformat()
    head = [f"# Autopilot digest, {today} (dry run on {result.here})", ""]
    if result.paused is not None:
        return "\n".join(head + [f"Paused: `AUTOPILOT_PAUSED` says: {result.paused}", ""])
    by_kind = {k: [d for d in result.decisions if d.kind == k] for k, _ in SECTIONS}
    counts = ", ".join(f"{len(by_kind[k])} {label.split(':')[0].lower()}" for k, label in SECTIONS if by_kind[k])
    out = head + [f"{len(result.decisions)} studies: {counts}.", ""]
    for w in result.warnings:
        out += [f"> warning: {w}", ""]
    for kind, label in SECTIONS:
        ds = by_kind[kind]
        if not ds:
            continue
        out += [f"## {label} ({len(ds)})", ""]
        for d in ds:
            where = f" on {d.machine}" if d.machine else ""
            if kind == "blocked":
                out.append(f"- **{d.key}** ({_arrow(d)}{where}): {d.question}")
                out.append(f"  - why: {d.reason}")
            elif kind == "done":
                out.append(f"- **{d.key}**: {d.reason}")
            else:
                out.append(f"- **{d.key}** ({_arrow(d)}{where}): {d.reason}")
                if kind == "advance" and d.playbook:
                    out.append(f"  - remedy: {d.action}")
                elif kind == "advance":
                    out.append(f"  - would run: `{d.action}`")
                elif kind == "waiting":
                    out.append(f"  - waiting for: {d.action}")
            if d.playbook:
                out.append(f"  - playbook: `{d.playbook}`")
            if kind != "done":
                out += [f"  - evidence: {e}" for e in d.evidence[:6]]
            out += [f"  - proposed (not applied): {p}" for p in d.proposals]
        out.append("")
    return "\n".join(out)
