"""The per-study state machine: one transition per study per tick.

    free → claimed → registered → fastq_staged → running → summarized → reviewed → merged
                                                     ↘ blocked (a person must answer one question)
    excluded (terminal)

`classify` maps what was read to a state; `decide` picks the single next
transition and says why. Deterministic and pure: everything it needs is in
the StudyState, the context from gather, the machine limits and the slot
counters, so every transition is unit-tested with hand-built states.
"""
from __future__ import annotations

import re

from whale.autopilot import playbook
from whale.autopilot.library import library_findings, readiness
from whale.autopilot.machines import (
    Machine,
    Slots,
    find_machine,
    free_gb,
    machines_for_profile,
)
from whale.autopilot.model import Decision, StudyState

RUN_STAGES = {"aligned", "called", "qc_passed"}
ACTIVE_STATES = {"claimed", "registered", "fastq_staged", "running"}
ENRICHMENT = {"medip", "mbd"}

# Backlog notes that mean a person must decide before a claim (PLAN.md §10 and the notes column).
NOTE_GATES = [
    (re.compile(r"claim in parts", re.I), "the backlog says to claim it in parts",
     "How should {sid} be split into parts (by tissue, by experiment, or 1ofN), and on which machines?"),
    (re.compile(r"taxon filter|run list first", re.I), "the backlog says it needs a run list or taxon filter",
     "Which runs of {acc} are C. gigas (a run list or taxon filter for fetch-samplesheet)?"),
]


def classify(s: StudyState) -> str:
    stage = s.stage
    if stage == "excluded":
        return "excluded"
    if s.merged or stage == "merged":
        return "merged"
    if s.branch is None or stage == "released":
        return "free"
    if stage in (None, "claimed"):
        return "claimed"
    if stage == "registered":
        return "registered"
    if stage == "fastq_staged":
        # Any Nextflow run of it here (even one we cannot tell is alive) means it is past launching.
        if s.run is not None and s.run.session:
            return "running"
        # A head job in its machine's Slurm queue (read this tick) means it was launched there.
        if s.head_jobs:
            return "running"
        return "fastq_staged"
    if stage in RUN_STAGES:
        return "running"
    if stage == "summarized":
        return "summarized"
    return "blocked"


def _sid(s: StudyState) -> str:
    return s.study_id or s.key


def _fastq_gb(s: StudyState) -> float | None:
    for v in ((s.study_yaml or {}).get("data_availability") or {}).get("fastq_size_gb"), \
             (s.backlog or {}).get("fastq_gb"):
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return None


def _from_match(s: StudyState, state: str, m: playbook.Match, machine: Machine | None,
                next_state: str | None) -> Decision:
    return Decision(key=s.key, state=state, next_state=next_state, kind="advance",
                    action=m.remedy, reason=f"known failure: {m.title}", machine=machine.name if machine else s.machine,
                    playbook=m.entry, evidence=m.evidence)


def _staging_command(s: StudyState, machine: Machine | None) -> str:
    cmd = "whale fetch-fastq --mark-staged"
    exports = f"export {machine.exports()}; " if machine is not None and machine.env else ""
    if machine is not None and machine.launch == "sbatch":
        return f"sbatch {machine.head_job} --wrap '{exports}{cmd}'"
    return f"{exports}{cmd}"


def _launch_command(s: StudyState, machine: Machine | None) -> str:
    cmd = playbook._relaunch(s, machine, resume=False)
    if machine is not None and machine.launch == "local":
        return f"tmux new -d -s {_sid(s)} \"{cmd}\""
    return cmd


def decide(s: StudyState, ctx: dict, machines: dict[str, Machine], *, here: str,
           runs: Slots, active: Slots, claimed_this_tick: set[str]) -> Decision:
    state = classify(s)
    machine = find_machine(machines, s.machine)
    sid = _sid(s)

    if state == "excluded":
        note = ((s.status or {}).get("history") or [{}])[-1].get("note") or (s.status or {}).get("notes")
        return Decision(s.key, state, None, "done", "nothing", f"excluded: {note or 'no reason recorded'}")

    if state == "merged":
        findings = [f for f in library_findings(s) if not f.ok]
        return Decision(s.key, state, None, "done", "nothing", "merged to main",
                        proposals=[f"{f.detail}; proposed: {f.proposal}" for f in findings if f.proposal],
                        evidence=[e for f in findings for e in f.evidence])

    if state == "free":
        return _decide_free(s, ctx, machines, active, claimed_this_tick)

    if state == "claimed":
        notes = (s.backlog or {}).get("notes") or ""
        how = "import the AREE registration, then " if "AREE" in notes else ""
        return Decision(
            s.key, state, "registered", "blocked",
            f"register: {how}fill registry/studies/{sid}.yaml, `whale fetch-samplesheet --bioproject "
            f"{s.accession} --show-attributes`, then with --condition-attribute; `whale validate-study --strict`",
            "registration needs someone to read the deposit (the register agent, step C)",
            question=f"Which sample attributes define the conditions of {sid} ({s.accession})?",
            machine=s.machine)

    if state in ("registered", "fastq_staged"):
        chem = playbook.chemistry_null(s, machine)
        if chem:
            if (s.study_yaml or {}).get("chemistry") is None:
                return Decision(s.key, state, None, "blocked", chem.remedy, chem.title,
                                question=f"What assay chemistry is {sid} (wgbs, emseq, pbat, rrbs)?",
                                machine=s.machine, playbook=chem.entry, evidence=chem.evidence)
            return _from_match(s, state, chem, machine, state)

    if state == "registered":
        # A staging log on this machine (fetch-fastq output) may show a retryable failure.
        if s.run is not None:
            m = playbook.transient_network(s, machine)
            if m:
                return _from_match(s, state, m, machine, "fastq_staged")
        need, have = _fastq_gb(s), free_gb(machine, here) if machine else None
        if machine and need is not None and have is not None and have < need * 1.02 + machine.disk_floor_gb:
            return Decision(s.key, state, None, "waiting", "free disk space",
                            f"{machine.name} has {have:.0f} GB free; staging {need:.0f} GB must leave "
                            f"{machine.disk_floor_gb:.0f} GB", machine=machine.name)
        where = "" if have is not None or not machine else f" (free space on {machine.name} not checked from {here})"
        return Decision(s.key, state, "fastq_staged", "advance", _staging_command(s, machine),
                        f"registered; stage {need or '?'} GB of FASTQ on {s.machine}{where}", machine=s.machine)

    if state == "fastq_staged":
        name = s.machine or "?"
        if not runs.has_room(name):
            m = find_machine(machines, name)
            return Decision(s.key, state, None, "waiting", "a free run slot",
                            f"{name} is running {runs.count(name)} of {m.max_concurrent_runs if m else '?'} studies",
                            machine=name)
        if machine is not None and machine.launch == "sbatch":
            if s.slurm_checked:
                runs.take(name)
                return Decision(s.key, state, "running", "advance", _launch_command(s, machine),
                                f"FASTQ staged and {name}'s queue has no job for it", machine=name)
            # Never propose a launch that might duplicate a head job we could not see.
            why = s.slurm_error or "its queue was not read"
            return Decision(s.key, state, None, "waiting", f"read the {name} queue",
                            f"FASTQ staged, but {name}'s queue could not be read from {here} ({why}); "
                            f"launch only once `squeue -u $USER` there shows no head job for it",
                            machine=name, evidence=[f"would run: {_launch_command(s, machine)}"])
        runs.take(name)
        if s.worktree is None and (machine is None or machine.name.lower() != here.lower()):
            return Decision(s.key, state, "running", "advance", _launch_command(s, machine),
                            f"FASTQ staged; no run visible from {here}: launch once no Nextflow is running it "
                            f"on {name}", machine=name)
        return Decision(s.key, state, "running", "advance", _launch_command(s, machine),
                        "FASTQ staged and nothing has run yet", machine=name)

    if state == "running":
        return _decide_running(s, machine)

    if state == "summarized":
        findings = readiness(s, ctx.get("main_assembly"))
        bad = [f for f in findings if not f.ok]
        evidence = [f"{f.check}: {f.detail}" for f in findings]
        if not bad:
            gh = "" if (machine is None or machine.gh) else " (mark the PR ready through the REST API: no gh here)"
            return Decision(s.key, state, "reviewed", "advance",
                            f"mark reviewed; merging the study PR stays with a person (gate 1){gh}",
                            "summarized and every readiness check passes", machine=s.machine, evidence=evidence)
        first = bad[0]
        if first.check == "pass_rate":
            q = f"{sid}: {first.detail}. Accept, re-run, or exclude the failing libraries?"
        elif first.check == "snp_mask":
            q = f"Re-run {sid} with the SNP mask (-resume) before it is merged?"
        elif first.proposal:
            q = f"{sid}: {first.detail}. Apply the proposed change and re-run?"
        else:
            q = f"{sid}: {first.detail}. Merge anyway?"
        return Decision(s.key, state, None, "blocked", "review the findings below", "readiness check failed",
                        question=q, machine=s.machine, evidence=evidence,
                        proposals=[f"{f.check}: {f.proposal}" for f in bad if f.proposal])

    return Decision(s.key, "blocked", None, "blocked", "read status.yaml",
                    f"unrecognized stage {s.stage!r}",
                    question=f"What stage is {sid} really at (status.yaml says {s.stage!r})?", machine=s.machine)


def _decide_free(s: StudyState, ctx: dict, machines: dict[str, Machine], active: Slots,
                 claimed_this_tick: set[str]) -> Decision:
    row = s.backlog or {}
    chem, notes = row.get("chemistry") or "", row.get("notes") or ""
    sid = s.study_id
    if chem in ENRICHMENT and "enrichment" not in ctx.get("workflows", set()):
        return Decision(s.key, "free", None, "waiting", "the enrichment workflow on main",
                        f"{chem} needs the enrichment workflow, which is not on main yet")
    for pattern, why, question in NOTE_GATES:
        if pattern.search(notes):
            return Decision(s.key, "free", None, "blocked", "decide how to claim it", why,
                            question=question.format(sid=sid or s.accession, acc=s.accession))
    if s.registered_on_main:
        return Decision(s.key, "free", None, "blocked", "claim its parts",
                        "registered on main but not claimed",
                        question=f"Which parts of {sid or s.accession} should be claimed, and where?")
    if sid is None:
        return Decision(s.key, "free", "claimed", "blocked", "name the study, then claim it",
                        "the backlog has no study_id for it yet (the register agent names it in step C)",
                        question=f"What study_id should {s.accession} get (FIRSTAUTHOR<YEAR>_<TOPIC>_<CHEMISTRY>)?")
    profile = row.get("suggested_profile") or "desktop"
    candidates = machines_for_profile(machines, profile)
    for m in candidates:
        if m.name.lower() in claimed_this_tick or not active.has_room(m.name):
            continue
        claimed_this_tick.add(m.name.lower())
        active.take(m.name)
        extra = " --library-strategy WGS at fetch-samplesheet" if "--library-strategy WGS" in notes else ""
        return Decision(s.key, "free", "claimed", "advance",
                        f"whale start-study --accession {s.accession} --study-id {sid} --profile {profile}",
                        f"next in the backlog for {profile}; {m.name} has room"
                        + (f" (note: pass{extra})" if extra else ""), machine=m.name)
    if not candidates:
        return Decision(s.key, "free", None, "waiting", f"a machine for profile {profile}",
                        f"config/autopilot.yaml has no machine running profile {profile}")
    why = []
    for m in candidates:
        if m.name.lower() in claimed_this_tick:
            why.append(f"{m.name} already takes a claim this tick")
        else:
            why.append(f"{m.name} has {active.count(m.name)} of {m.max_concurrent_runs} studies in progress")
    return Decision(s.key, "free", None, "waiting", f"a free {profile} slot", "; ".join(why))


def _decide_running(s: StudyState, machine: Machine | None) -> Decision:
    run, sid = s.run, _sid(s)
    name = s.machine
    matches = playbook.match(s, machine)

    if run is None:
        gannet = playbook.gannet_upload(s, machine)
        if gannet:
            return _from_match(s, "running", gannet, machine, "summarized")
        passed = (s.summary or {}).get("passed")
        if s.stage == "called" and passed == 0:
            return Decision(s.key, "running", None, "blocked", "decide on the study",
                            "no library passed the QC gate",
                            question=f"No library of {sid} passed QC: exclude it, or re-run with other settings?",
                            machine=name, evidence=[f"stage called, 0 of {(s.summary or {}).get('libraries')} pass"])
        if s.head_jobs:
            head = s.head_jobs[0]
            return Decision(s.key, "running", None, "running", "wait",
                            f"head job {head.job_id} is {head.state} on {name} ({head.elapsed})",
                            machine=name, evidence=[j.describe() for j in s.slurm_jobs])
        where = f"; {name}'s queue has no head job for it" if s.slurm_checked else ""
        return Decision(s.key, "running", None, "running", "wait",
                        f"stage {s.stage} on {name}; no run log readable from here{where}", machine=name)

    hazards = [m for m in matches if m.hazard]
    if run.status == "live":
        if hazards:
            return _from_match(s, "running", hazards[0], machine, "running")
        done = len(run.cached)
        return Decision(s.key, "running", None, "running", "wait",
                        f"run {run.run_name} is live on {name}" + (f" ({done} tasks cached)" if done else ""),
                        machine=name, evidence=[run.source])

    if run.status == "ok":
        if hazards:
            return _from_match(s, "running", hazards[0], machine, "running")
        gannet = playbook.gannet_upload(s, machine)
        if gannet:
            return _from_match(s, "running", gannet, machine, "summarized")
        results = f"$WHALE_METHYLOMES/{sid}" + (f"/parts/{s.part}" if s.part else "")
        return Decision(s.key, "running", "summarized", "advance", f"whale record-run --results {results}",
                        f"run {run.run_name} finished but the study is at {s.stage}", machine=name,
                        evidence=[run.source])

    if run.status in ("failed", "stopped"):
        if matches:
            return _from_match(s, "running", matches[0], machine, "running")
        task = run.failed[0] if run.failed else None
        what = f"{task.process} ({task.tag}) failed with exit {task.exit_status}" if task else "the run failed"
        return Decision(s.key, "running", None, "blocked", "needs diagnosis", f"no playbook entry matches: {what}",
                        question=f"Why did {what} in {sid}, and what should change before a relaunch?",
                        machine=name, evidence=playbook.diagnosis_evidence(s))

    return Decision(s.key, "running", None, "running", "wait",
                    f"run {run.run_name or '?'} has not finished and has no pid file to check", machine=name,
                    evidence=[run.source])
