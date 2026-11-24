"""Known failures and their remedies: deterministic, no model.

Each entry looks at what was read about a study (its status, config and the
latest run this machine can see) and, when it recognizes the failure, returns
the evidence it matched and a concrete remedy. The tick proposes the remedy;
in step A nothing is applied. A failed run no entry matches is reported as
"needs diagnosis" with its evidence (step C hands that to a model).

To add an entry: a function `(study, machine) -> Match | None`, registered in
ENTRIES, with a fixture test in tests/test_autopilot_playbook.py.
"""
from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from whale.autopilot.machines import Machine
from whale.autopilot.model import StudyState
from whale.claims import config_chemistry, config_name


@dataclass
class Match:
    entry: str
    title: str
    remedy: str
    evidence: list[str] = field(default_factory=list)
    # A hazard is found on a run that did not fail (or is still running) but whose
    # results would be wrong; it must be dealt with before the run is recorded.
    hazard: bool = False


def _text(study: StudyState) -> str:
    run = study.run
    return f"{run.text}\n{run.extra_text}" if run else ""


def _relaunch(study: StudyState, machine: Machine | None, *, resume: bool = True) -> str:
    """The command a person would use to relaunch the study's run on its machine."""
    sid, part = study.study_id or study.key, study.part
    cfgs = [config_name(sid, part)]
    if study.run_config_text is not None and machine is not None:
        cfgs.append(f"data/studies/{sid}/{machine.name.lower()}.run.config")
    if machine is not None:
        cfgs += machine.extra_configs
    profile = (machine.profile if machine else study.profile) or "desktop"
    prof = f"{profile},docker" if profile == "desktop" else profile
    session = study.run.session if (resume and study.run and study.run.session) else ""
    cmd = "nextflow run workflows/methylation " + " ".join(f"-c {c}" for c in cfgs) + f" -profile {prof}"
    if resume:
        cmd += f" -resume {session}".rstrip()
    # The SNP mask is per assembly, and each machine keeps its masks in its own place.
    mask = machine.snp_mask(study.snp_mask_id, study.genome_assembly) if machine is not None else None
    extra = {"WHALE_SNP_MASK": mask} if mask else {}
    if machine is not None and machine.launch == "sbatch":
        return machine.sbatch(cmd, study_id=sid, part=part, kind="head", env=extra)
    exports = f"export {machine.exports(extra)}; " if machine is not None and (machine.env or extra) else ""
    return f"{exports}{cmd}"


# --- config helpers ----------------------------------------------------------

def _strip_comments(text: str) -> str:
    return re.sub(r"//[^\n]*", "", re.sub(r"/\*.*?\*/", "", text, flags=re.S))


def config_list(text: str | None, name: str) -> list[str]:
    """The quoted strings of a `name = [ ... ]` list in a Nextflow config."""
    if not text:
        return []
    m = re.search(rf"\b{name}\s*=\s*\[(.*?)\]", _strip_comments(text), re.S)
    return re.findall(r"""['"]([^'"]+)['"]""", m.group(1)) if m else []


def config_string(text: str | None, name: str) -> str | None:
    if not text:
        return None
    m = re.search(rf"""\b{name}\s*=\s*(?:'([^']*)'|"([^"]*)")""", _strip_comments(text))
    return (m.group(1) if m.group(1) is not None else m.group(2)) if m else None


# --- entries -----------------------------------------------------------------

SWAP_PROCESSES = ("TRIM_GALORE", "BISMARK_ALIGN")


def swap_mates_resume(study: StudyState, machine: Machine | None) -> Match | None:
    """swap_mates set, then -resume: Nextflow reuses trim/align tasks, since file order is not hashed."""
    run = study.run
    swapped = set(config_list(study.config_text, "swap_mates"))
    if not run or not run.resumed or not swapped:
        return None
    stale = [c for c in run.cached if c.process in SWAP_PROCESSES and c.tag in swapped]
    if not stale:
        return None
    wd = run.workdir or "$WHALE_WORKDIR"
    aside = f"{wd.rstrip('/')}-unswapped"
    moves = "; ".join(f"mv {wd}/{c.hash}* {aside}/" for c in stale[:6])
    more = f" (and {len(stale) - 6} more)" if len(stale) > 6 else ""
    return Match(
        "swap_mates_resume", "swap_mates + -resume reused unswapped trim/align tasks",
        f"cancel the run before record-run; move the {len(stale)} cached TRIM_GALORE/BISMARK_ALIGN "
        f"work dirs of swapped libraries aside (mkdir -p {aside}; {moves}{more}), then relaunch: "
        f"{_relaunch(study, machine)}",
        [f"{run.source}: resumed run {run.run_name} cached {c.process} ({c.tag}) [{c.hash}], "
         f"a library in swap_mates" for c in stale[:4]],
        hazard=run.status in ("live", "ok"),
    )


ENV_VARS = ("WHALE_RAW", "WHALE_REFERENCE", "WHALE_METHYLOMES", "WHALE_WORKDIR")


def resume_env_cache_miss(study: StudyState, machine: Machine | None) -> Match | None:
    """-resume that cached nothing: WHALE_* were not exported, so inputs hash differently."""
    run = study.run
    if not run or not run.resumed or run.cached:
        return None
    evidence = [f"{run.source}: run {run.run_name} was launched with -resume but cached no task"]
    want = (machine.env if machine else {}).get("WHALE_WORKDIR")
    wd = run.workdir or ""
    checkout = study.worktree or ""
    suspicious = False
    if wd and checkout and wd.startswith(checkout.rstrip("/") + "/"):
        evidence.append(f"Work-dir {wd} is inside the checkout: WHALE_WORKDIR was not set")
        suspicious = True
    elif wd and want and "$" not in want and wd.rstrip("/") != want.rstrip("/"):
        evidence.append(f"Work-dir {wd} is not this machine's WHALE_WORKDIR ({want})")
        suspicious = True
    if run.previous_runs:
        evidence.append(f"{run.previous_runs} earlier run(s) of session {run.session} in .nextflow/history")
        suspicious = True
    if not suspicious:
        return None
    return Match(
        "resume_env_cache_miss", "-resume missed the cache (WHALE_* not exported)",
        f"stop the run if it is still going; export {', '.join(ENV_VARS)} and relaunch with the "
        f"session id: {_relaunch(study, machine)}; check for 'Cached process' lines before leaving it",
        evidence, hazard=run.status == "live",
    )


_SLURM_ACCOUNT = re.compile(r"Invalid account or account/partition combination|Invalid account|"
                            r"Invalid partition name|QOSMax\w+Limit|AssocGrp\w+Limit")
_PREEMPTED = re.compile(r"PREEMPTED|terminated for an unknown reason|Job was preempted|DUE TO PREEMPTION")
_SCONTROL = re.compile(r"scontrol\s+update\b[^\n]*(?:Partition|Account)=", re.I)


def klone_partition(study: StudyState, machine: Machine | None) -> Match | None:
    """klone: wrong account for a partition, repeated preemption, or tasks moved with scontrol."""
    if (machine.profile if machine else study.profile) != "klone":
        return None
    text = _text(study)
    evidence = []
    if (m := _SLURM_ACCOUNT.search(text)):
        evidence.append(f"Slurm refused a job: {m.group(0)!r}")
    preempted = len(_PREEMPTED.findall(text))
    if preempted >= 3:
        evidence.append(f"{preempted} preemption / external-termination messages in the log")
    if (m := _SCONTROL.search(text)):
        evidence.append(f"a task was moved with scontrol ({m.group(0).strip()!r}); Nextflow resubmits moved tasks")
    if not evidence:
        return None
    table = ", ".join(f"{a['partition']}={a['account']}" for a in (machine.slurm_accounts if machine else []))
    sid = study.study_id or study.key
    return Match(
        "klone_partition", "klone partition/account problem",
        f"put a per-process override in data/studies/{sid}/klone.run.config that sets the account with "
        f"the partition ({table or 'see docs/compute_profiles.md'}), capped with maxForks to the lab's "
        f"share; never move Nextflow tasks with scontrol; relaunch the head job: {_relaunch(study, machine)}",
        evidence,
    )


def chemistry_null(study: StudyState, machine: Machine | None) -> Match | None:
    """start-study writes `chemistry = null`; the run refuses it only after staging."""
    if study.config_text is None:
        return None
    have = config_chemistry(study.config_text)
    want = (study.study_yaml or {}).get("chemistry")
    if have is not None and (want is None or have == want):
        return None
    name = config_name(study.study_id or study.key, study.part)
    if have is None:
        ev = [f"{name} sets chemistry = null"]
    else:
        ev = [f"{name} sets chemistry = {have!r} but the study YAML says {want!r}"]
    if want is None:
        return Match("chemistry_null", "config chemistry not set",
                     "fill in `chemistry` in the study YAML first, then set the same term in the config",
                     ev + ["the study YAML has no chemistry either"])
    return Match("chemistry_null", "config chemistry not set",
                 f"set chemistry = '{want}' in {name} (to match the study YAML) and commit it on the study branch",
                 ev)


_NO_GH = re.compile(r"gh: (?:command )?not found|No such file or directory: 'gh'")


def raven_no_gh(study: StudyState, machine: Machine | None) -> Match | None:
    """raven has no `gh`; a PR step run there silently does nothing."""
    on_raven = (machine and not machine.gh) or (study.machine or "").lower() == "raven"
    if not on_raven or not _NO_GH.search(_text(study)):
        return None
    return Match(
        "raven_no_gh", "GitHub CLI missing on raven",
        "do the GitHub step through the REST/GraphQL API with the repo token instead of gh "
        "(mark ready: GraphQL markPullRequestReadyForReview on the study PR)",
        [f"{study.run.source if study.run else 'log'}: gh not found on {study.machine or 'raven'}"],
    )


_UPLOAD = re.compile(r"methylomes not on the data host[^\n]*|Permission denied \(publickey[^\n]*|"
                     r"Host key verification failed[^\n]*")


def gannet_upload(study: StudyState, machine: Machine | None) -> Match | None:
    """record-run could not upload the methylomes (no SSH key for gannet on that machine)."""
    on_host = ((study.summary or {}).get("methylomes") or {}).get("on_data_host")
    m = _UPLOAD.search(_text(study))
    stuck = study.stage == "qc_passed" and on_host is False
    if not (stuck or (m and study.stage in ("aligned", "called", "qc_passed"))):
        return None
    evidence = []
    if stuck:
        evidence.append("stage qc_passed and study_summary.yaml methylomes.on_data_host is false")
    if m:
        evidence.append(f"log: {m.group(0).strip()[:160]}")
    sid = study.study_id or study.key
    results = f"$WHALE_METHYLOMES/{sid}" + (f"/parts/{study.part}" if study.part else "")
    return Match(
        "gannet_upload", "methylome upload to gannet failed",
        f"rerun `whale record-run --results {results}` from the study worktree on a machine with a gannet "
        "SSH key (raven), with the results readable there; it only moves stages forward",
        evidence,
    )


_NETWORK = re.compile(
    r"Connection reset by peer|Connection timed out|timed out|Temporary failure in name resolution|"
    r"HTTP Error 5\d\d|50[234] (?:Bad Gateway|Service Unavailable|Gateway Time-?out)|IncompleteRead|"
    r"RemoteDisconnected|curl: \((?:6|7|18|28|35|52|56)\)|Network is unreachable|ConnectionError")


def transient_network(study: StudyState, machine: Machine | None) -> Match | None:
    """ENA or network hiccup while staging; fetch-fastq resumes and verifies, so retry."""
    m = _NETWORK.search(_text(study))
    if not m:
        return None
    if study.run and study.run.failed and study.run.failed[0].exit_status not in (None, 1, 6, 7, 28, 56):
        return None
    attempts = sum(1 for line in _text(study).splitlines() if _NETWORK.search(line))
    wait = min(10 * 2 ** min(attempts - 1, 5), 360)
    return Match(
        "transient_network", "transient ENA/network error",
        f"retry `whale fetch-fastq --mark-staged` after ~{wait} min (it resumes and MD5-checks); "
        f"prefer the gannet mirror (--source mirror) and, if the accession is not mirrored yet, "
        f"mirror it on raven first (`whale mirror-fastq --accession {study.accession}`)",
        [f"{m.group(0)!r} ({attempts} network error line(s))"],
    )


_OOM = re.compile(r"exit status \(137\)|\bKilled\b|OutOfMemoryError|Cannot allocate memory|oom-kill")


def oom_killed(study: StudyState, machine: Machine | None) -> Match | None:
    """A task killed for memory (exit 137): on Hete, two big tasks in the 12 GB Colima VM."""
    run = study.run
    if not run or run.status != "failed" or not run.failed:
        return None
    task = run.failed[0]
    if task.exit_status != 137 and not _OOM.search(task.stderr or ""):
        return None
    where = (machine.limits_config if machine and machine.limits_config
             else f"data/studies/{study.study_id or study.key}/{(machine.name if machine else 'machine').lower()}.run.config")
    evidence = [f"{task.process} ({task.tag}) exit {task.exit_status}" + (f" in {task.workdir}" if task.workdir else "")]
    if _already_capped(where, task.process):
        evidence.append(f"{where} already runs {task.process} one at a time")
        remedy = (f"{task.process} is already capped at maxForks = 1 in {where}; give it more memory there "
                  f"(the machine's container VM permitting) or move the study to a bigger machine "
                  f"(`whale set-status --profile raven --force`), then relaunch: {_relaunch(study, machine)}")
    else:
        remedy = (f"in {where} cap {task.process} with `withName: '{task.process}' {{ maxForks = 1 }}` (and raise "
                  f"its memory if the machine has room), then relaunch: {_relaunch(study, machine)}")
    return Match("oom_killed", "task killed for memory (exit 137)", remedy, evidence)


def _already_capped(config_path: str, process: str) -> bool:
    """Whether a readable config already gives `process` maxForks = 1."""
    try:
        text = _strip_comments(Path(config_path).read_text())
    except OSError:
        return False
    m = re.search(rf"withName:\s*['\"]?{re.escape(process)}['\"]?\s*\{{([^}}]*)\}}", text)
    return bool(m and re.search(r"maxForks\s*=\s*1\b", m.group(1)))


def head_stopped(study: StudyState, machine: Machine | None) -> Match | None:
    """The head process is gone and the log has no end or error: a power cut, reboot or kill."""
    run = study.run
    if not run or run.status != "stopped" or run.failed:
        return None
    return Match(
        "head_stopped", "run stopped without finishing (head process gone)",
        f"relaunch with -resume and this machine's WHALE_* exported: {_relaunch(study, machine)}",
        [f"{run.source}: run {run.run_name} has no end and its .nextflow.pid process is gone"],
    )


# Most specific first; the first match is the remedy.
ENTRIES: list[Callable[[StudyState, Machine | None], Match | None]] = [
    swap_mates_resume,
    resume_env_cache_miss,
    gannet_upload,
    raven_no_gh,
    klone_partition,
    chemistry_null,
    oom_killed,
    transient_network,
    head_stopped,
]


def match(study: StudyState, machine: Machine | None) -> list[Match]:
    return [m for entry in ENTRIES if (m := entry(study, machine))]


def diagnosis_evidence(study: StudyState) -> list[str]:
    """What a person (or step C's model) needs to diagnose a failure no entry matched."""
    run = study.run
    if run is None:
        return ["no run log readable on this machine"]
    out = [f"log: {run.source} (run {run.run_name or '?'}, session {run.session or '?'})"]
    for t in run.failed[:3]:
        out.append(f"failed: {t.process} ({t.tag}) exit {t.exit_status}" + (f", work dir {t.workdir}" if t.workdir else ""))
        tail = [ln for ln in (t.stderr or "").splitlines() if ln.strip()][-3:]
        out += [f"  stderr: {ln.strip()[:200]}" for ln in tail]
    if not run.failed:
        tail = [ln for ln in run.text.splitlines() if re.search(r"ERROR|WARN|Exception", ln)][-3:]
        out += [f"  {ln.strip()[:200]}" for ln in tail]
    return out
