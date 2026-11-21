"""What the autopilot knows about a study, and what it decides to do with it."""
from __future__ import annotations

from dataclasses import dataclass, field

# The autopilot's own states. They are coarser than pipeline_stage.yaml:
# aligned/called/qc_passed are written in one go by `whale record-run`, so a
# study at one of them is still "running" until it reaches summarized.
STATES = ["free", "claimed", "registered", "fastq_staged", "running", "summarized",
          "reviewed", "merged", "blocked", "excluded"]

# Forward order of the main line (blocked and excluded sit outside it).
MAIN_LINE = ["free", "claimed", "registered", "fastq_staged", "running", "summarized",
             "reviewed", "merged"]

# How the digest groups decisions.
KINDS = ["advance", "running", "blocked", "waiting", "done"]


@dataclass
class CachedTask:
    hash: str                         # work dir prefix, e.g. 40/7189a7
    process: str
    tag: str | None = None


@dataclass
class FailedTask:
    process: str
    tag: str | None = None
    exit_status: int | None = None
    workdir: str | None = None
    stderr: str = ""


@dataclass
class RunEvidence:
    """The latest Nextflow run of a study that this machine can read."""
    source: str                       # the log file read
    status: str                       # live | ok | failed | stopped | unknown
    run_name: str | None = None
    session: str | None = None
    command: str | None = None
    workdir: str | None = None
    resumed: bool = False
    cached: list[CachedTask] = field(default_factory=list)
    failed: list[FailedTask] = field(default_factory=list)
    previous_runs: int = 0            # earlier runs of the same session in .nextflow/history
    text: str = ""                    # the log text, for playbook patterns
    extra_text: str = ""              # other logs beside it (fetch-fastq, sbatch .out files)


@dataclass
class StudyState:
    """Everything read about one study (or study part). Built by gather; pure data."""
    key: str                          # STUDY_ID or STUDY_ID.part-PART; the backlog accession if unnamed
    accession: str
    study_id: str | None = None
    part: str | None = None
    backlog: dict | None = None       # registry/backlog.tsv row
    branch: str | None = None         # study/<ACC>/<ID> if claimed
    merged: bool = False              # the branch (or its data) is on main
    registered_on_main: bool = False
    status: dict | None = None        # status.yaml (branch, else main)
    study_yaml: dict | None = None
    config_text: str | None = None    # config/<ID>.config on the branch
    run_config_text: str | None = None  # data/studies/<ID>/<machine>.run.config, if any
    qc: list[dict] | None = None      # qc_summary.tsv
    summary: dict | None = None       # study_summary.yaml
    run_manifest: dict | None = None
    run: RunEvidence | None = None
    worktree: str | None = None       # local checkout of the branch, if on this machine

    @property
    def stage(self) -> str | None:
        return (self.status or {}).get("stage")

    @property
    def machine(self) -> str | None:
        return (self.status or {}).get("machine")

    @property
    def profile(self) -> str | None:
        return (self.status or {}).get("profile") or (self.backlog or {}).get("suggested_profile")


@dataclass
class Decision:
    """The one transition a tick would make for a study, and why."""
    key: str
    state: str                        # one of STATES
    next_state: str | None            # None when nothing moves
    kind: str                         # one of KINDS
    action: str                       # what would be done (a command, or what a person must do)
    reason: str
    question: str | None = None       # for blocked: the one thing a person must answer
    machine: str | None = None
    playbook: str | None = None       # id of the matched playbook entry
    evidence: list[str] = field(default_factory=list)
    proposals: list[str] = field(default_factory=list)  # config changes proposed, never applied

    def __post_init__(self):
        if self.state not in STATES:
            raise ValueError(f"unknown state {self.state!r}")
        if self.next_state is not None and self.next_state not in STATES:
            raise ValueError(f"unknown next state {self.next_state!r}")
        if self.kind not in KINDS:
            raise ValueError(f"unknown kind {self.kind!r}")
        if self.kind == "blocked" and not self.question:
            raise ValueError(f"{self.key}: a blocked decision needs a question")
