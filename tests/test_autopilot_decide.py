"""Every transition of the autopilot state machine, on hand-built study states."""
import pytest
from autopilot_fakes import CONFIG, nf_log, qc_rows, run_of, state, summary

from whale.autopilot.decide import classify, decide
from whale.autopilot.machines import Machine, Slots
from whale.autopilot.model import Decision, StudyState
from whale.autopilot.slurm import SlurmJob

MACHINES = {
    "Hete": Machine(name="Hete", profile="desktop", max_concurrent_runs=1, disk_floor_gb=100,
                    env={"WHALE_WORKDIR": "/w/work"}),
    "raven": Machine(name="raven", profile="raven", max_concurrent_runs=2, gh=False, gannet_key=True),
    "klone": Machine(name="klone", profile="klone", max_concurrent_runs=2, launch="sbatch",
                     head_job="--account=coenv --partition=compute"),
}
CTX = {"main_assembly": "GCF_963853765.1", "workflows": {"methylation", "snp_mask"}}


def run(s, *, here="Hete", runs=None, active=None, claimed=None, ctx=CTX):
    return decide(s, ctx, MACHINES, here=here, runs=Slots(MACHINES, runs or {}),
                  active=Slots(MACHINES, active or {}), claimed_this_tick=set() if claimed is None else claimed)


def free(key="TESTER2099_FREE_WGBS", *, profile="raven", chemistry="wgbs", notes="", named=True):
    row = {"accession": "PRJNA999909", "suggested_study_id": key if named else "(name on claim)",
           "chemistry": chemistry, "suggested_profile": profile, "notes": notes, "fastq_gb": "50"}
    return StudyState(key=key if named else "PRJNA999909", accession="PRJNA999909",
                      study_id=key if named else None, backlog=row)


# --- free ----------------------------------------------------------------------

def test_free_to_claimed():
    d = run(free())
    assert (d.state, d.next_state, d.kind, d.machine) == ("free", "claimed", "advance", "raven")
    assert d.action == "whale start-study --accession PRJNA999909 --study-id TESTER2099_FREE_WGBS --profile raven"


def test_free_one_claim_per_machine_per_tick():
    claimed = set()
    assert run(free("A2099_ONE_WGBS"), claimed=claimed).kind == "advance"
    d = run(free("B2099_TWO_WGBS"), claimed=claimed)
    assert d.kind == "waiting" and "already takes a claim" in d.reason


def test_free_waits_when_machine_full():
    d = run(free(profile="desktop"), active={"Hete": 1})
    assert d.kind == "waiting" and "1 of 1" in d.reason


def test_free_blocked_without_a_name():
    d = run(free(named=False))
    assert d.kind == "blocked" and "study_id" in d.question and d.next_state == "claimed"


@pytest.mark.parametrize("notes,word", [("Largest set; claim in parts once read.", "parts"),
                                        ("needs a taxon filter or run list first", "runs")])
def test_free_blocked_by_backlog_note(notes, word):
    d = run(free(notes=notes))
    assert d.kind == "blocked" and word in d.question


def test_free_enrichment_waits_for_workflow():
    assert run(free(chemistry="medip")).kind == "waiting"
    assert run(free(chemistry="medip"), ctx={**CTX, "workflows": {"methylation", "enrichment"}}).kind == "advance"


def test_free_registered_on_main_unclaimed():
    s = free()
    s.registered_on_main = True
    assert run(s).kind == "blocked"


def test_free_unknown_profile_waits():
    d = run(free(profile="cloud"))
    assert d.kind == "waiting" and "no machine" in d.reason


# --- claimed / registered ------------------------------------------------------

def test_claimed_to_registered_needs_a_person():
    d = run(state("claimed"))
    assert (d.state, d.next_state, d.kind) == ("claimed", "registered", "blocked")
    assert "conditions" in d.question and "fetch-samplesheet" in d.action


def test_registered_to_fastq_staged():
    s = state("registered", machine="klone", profile="klone")
    d = run(s)
    assert (d.state, d.next_state, d.kind) == ("registered", "fastq_staged", "advance")
    assert d.action.startswith("sbatch --account=coenv") and "whale fetch-fastq --mark-staged" in d.action
    assert "--job-name=whale-fastq-TESTER2099_X_WGBS" in d.action


def test_registered_chemistry_null_is_fixed_first():
    d = run(state("registered", config=CONFIG.replace("'wgbs'", "null")))
    assert d.playbook == "chemistry_null" and d.kind == "advance" and "chemistry = 'wgbs'" in d.action


def test_registered_without_any_chemistry_blocks():
    d = run(state("registered", config=CONFIG.replace("'wgbs'", "null"), chemistry=None))
    assert d.kind == "blocked" and "chemistry" in d.question


def test_registered_waits_for_disk(monkeypatch):
    monkeypatch.setattr("whale.autopilot.decide.free_gb", lambda m, here: 120.0)
    s = state("registered", backlog={"fastq_gb": "50"})
    d = run(s)
    assert d.kind == "waiting" and "120 GB free" in d.reason
    monkeypatch.setattr("whale.autopilot.decide.free_gb", lambda m, here: 500.0)
    assert run(s).kind == "advance"


def test_registered_staging_network_error_retries():
    s = state("registered", run=run_of("", extra_text="curl: (28) Connection timed out"))
    s.run.status = "unknown"
    d = run(s)
    assert d.playbook == "transient_network" and d.next_state == "fastq_staged"


# --- fastq_staged --------------------------------------------------------------

def test_fastq_staged_to_running_here():
    s = state("fastq_staged", worktree="/checkout")
    d = run(s)
    assert (d.state, d.next_state, d.kind) == ("fastq_staged", "running", "advance")
    # Always -resume: harmless on a first launch, and keeps the work of an earlier failed run.
    assert d.action.startswith("tmux new -d -s TESTER2099_X_WGBS") and d.action.endswith(' -resume"')
    assert d.evidence == []


def test_launch_resumes_the_session_seen_on_this_machine():
    s = state("fastq_staged", worktree="/checkout",
              run=run_of(nf_log(command="nextflow run x", complete=False), session="f97ff558-1da0"))
    s.run.status = "unknown"
    s.status["stage"] = "fastq_staged"
    from whale.autopilot.decide import _launch_command
    assert "-resume f97ff558-1da0" in _launch_command(s, MACHINES["Hete"])


def test_launch_warns_about_swap_mates():
    s = state("fastq_staged", worktree="/checkout",
              config=CONFIG + "\nparams { swap_mates = ['S1', 'S2'] }\n")
    d = run(s)
    assert "-resume" in d.action and "swap_mates is set (S1, S2)" in d.evidence[0]


def _job(name="whale-head-TESTER2099", state="RUNNING", workdir="/gscratch/x/white-whale-TESTER2099_X_WGBS"):
    return SlurmJob("43445189", name, "compute", "compute-coenv", state, "10:17:55", workdir)


def test_fastq_staged_on_klone_queue_unread_waits():
    # The queue could not be read: never propose a launch that might duplicate a head job.
    s = state("fastq_staged", machine="klone", profile="klone", slurm_error="Permission denied (publickey)")
    d = run(s, runs={})
    assert (d.kind, d.next_state) == ("waiting", None)
    assert "Permission denied" in d.reason and "squeue" in d.reason
    assert d.evidence[0].startswith("would run: sbatch")


def test_fastq_staged_on_klone_empty_queue_launches():
    s = state("fastq_staged", machine="klone", profile="klone", slurm_checked=True)
    d = run(s)
    assert (d.kind, d.next_state) == ("advance", "running") and d.action.startswith("sbatch")
    assert "queue has no job" in d.reason
    # An empty queue may mean the last run failed (DECICOMP2025, disk quota, 2026-11-21): never start over.
    assert d.action.endswith(" -resume'")


def test_fastq_staged_with_head_job_is_running():
    s = state("fastq_staged", machine="klone", profile="klone", slurm_checked=True, slurm_jobs=[_job()])
    assert classify(s) == "running"
    d = run(s)
    assert (d.state, d.kind, d.action) == ("running", "running", "wait")
    assert "head job 43445189 is RUNNING on klone" in d.reason


def test_klone_head_job_takes_a_run_slot():
    from whale.autopilot.tick import decide_all
    running = state("fastq_staged", key="TESTER2099_A_WGBS", machine="klone", profile="klone",
                    slurm_checked=True, slurm_jobs=[_job()])
    other = [state("fastq_staged", key=f"TESTER2099_{x}_WGBS", machine="klone", profile="klone",
                   slurm_checked=True) for x in "BC"]
    kinds = [d.kind for d in decide_all([running, *other], CTX, MACHINES, "Hete")]
    assert kinds == ["running", "advance", "waiting"]   # klone allows 2 concurrent runs


def test_fastq_staged_waits_for_slot():
    d = run(state("fastq_staged", machine="raven", profile="raven"), runs={"raven": 2})
    assert d.kind == "waiting" and "2 of 2" in d.reason


def test_launch_takes_a_slot():
    runs = Slots(MACHINES, {})
    for expected in ("advance", "waiting"):
        d = decide(state("fastq_staged", worktree="/c"), CTX, MACHINES, here="Hete", runs=runs,
                   active=Slots(MACHINES), claimed_this_tick=set())
        assert d.kind == expected


def test_fastq_staged_chemistry_checked_before_launch():
    assert run(state("fastq_staged", config=CONFIG.replace("'wgbs'", "null"))).playbook == "chemistry_null"


# --- running -------------------------------------------------------------------

def test_running_live_waits():
    s = state("fastq_staged", run=run_of(nf_log(complete=False), status="live"))
    assert classify(s) == "running"
    d = run(s)
    assert (d.state, d.kind, d.next_state) == ("running", "running", None)


def test_running_live_with_hazard():
    cfg = CONFIG.replace("}", "    swap_mates = ['S1']\n}")
    log = nf_log(command="nextflow run x -resume", cached=[("40/7189a7", "TRIM_GALORE", "S1")], complete=False)
    d = run(state("fastq_staged", config=cfg, run=run_of(log, status="live")))
    assert d.playbook == "swap_mates_resume" and d.kind == "advance"


def test_running_ok_to_summarized():
    d = run(state("fastq_staged", run=run_of(nf_log())))
    assert (d.state, d.next_state, d.kind) == ("running", "summarized", "advance")
    assert d.action == "whale record-run --results $WHALE_METHYLOMES/TESTER2099_X_WGBS"


def test_running_failed_known():
    d = run(state("fastq_staged", run=run_of(nf_log(failure={"process": "BISMARK_DEDUPLICATE", "exit": 137}))))
    assert d.playbook == "oom_killed" and d.next_state == "running" and d.kind == "advance"


def test_running_failed_unknown_needs_diagnosis():
    d = run(state("fastq_staged", run=run_of(nf_log(failure={"process": "QC_GATE", "exit": 2}))))
    assert d.kind == "blocked" and d.playbook is None
    assert "QC_GATE" in d.question and any("exit 2" in e for e in d.evidence)


def test_running_stopped_relaunches():
    d = run(state("fastq_staged", run=run_of(nf_log(complete=False), status="stopped")))
    assert d.playbook == "head_stopped" and "-resume" in d.action


def test_running_upload_failed():
    d = run(state("qc_passed", machine="klone", summary=summary(on_host=False)))
    assert d.playbook == "gannet_upload" and d.next_state == "summarized"


def test_running_nothing_passed_blocks():
    d = run(state("called", summary=summary(passed=0)))
    assert d.kind == "blocked" and "exclude" in d.question


def test_running_remote_unseen():
    d = run(state("aligned", machine="klone", summary=summary(on_host=True)))
    assert d.kind == "running" and "no run log" in d.reason


def test_running_unknown_log():
    d = run(state("fastq_staged", run=run_of(nf_log(complete=False), status="unknown")))
    # An unfinished log with no pid file to check: never relaunched on a guess.
    assert (d.state, d.kind) == ("running", "running") and "no pid file" in d.reason


# --- summarized / reviewed / merged / excluded --------------------------------

def test_summarized_to_reviewed():
    d = run(state("summarized", machine="raven", profile="raven", qc=qc_rows(), summary=summary()))
    assert (d.state, d.next_state, d.kind) == ("summarized", "reviewed", "advance")
    assert "gate 1" in d.action and "REST" in d.action   # raven has no gh


def test_summarized_low_pass_rate_blocks():
    d = run(state("summarized", qc=qc_rows(48, passed=13), summary=summary(48, 13)))
    assert d.kind == "blocked" and "13 of 48" in d.question


def test_summarized_unmasked_blocks():
    d = run(state("summarized", qc=qc_rows(), summary=summary(snp_mask=False)))
    assert d.kind == "blocked" and "SNP mask" in d.question


def test_summarized_library_problem_proposes_change():
    d = run(state("summarized", qc=qc_rows(noncpg=3.0), summary=summary()))
    assert d.kind == "blocked" and "proposed change" in d.question
    assert any("--clip_r2 18" in p for p in d.proposals)


def test_merged_done_with_proposals():
    s = state("summarized", merged=True, qc=qc_rows(mapping=30.0), summary=summary())
    d = run(s)
    assert (d.state, d.kind) == ("merged", "done")
    assert any("--non_directional" in p for p in d.proposals)


def test_excluded_done():
    d = run(state("excluded", note="RRBS with MspI fragments < 40 bp"))
    assert (d.state, d.kind) == ("excluded", "done") and "MspI" in d.reason


def test_unknown_stage_blocks():
    d = run(state("teleported"))
    assert d.state == "blocked" and d.kind == "blocked"


# --- the model -----------------------------------------------------------------

def test_decision_validation():
    with pytest.raises(ValueError):
        Decision("X", "nowhere", None, "done", "", "")
    with pytest.raises(ValueError):
        Decision("X", "claimed", "registered", "blocked", "", "")   # blocked without a question
    with pytest.raises(ValueError):
        Decision("X", "claimed", None, "sleeping", "", "")


def test_classify_main_line():
    assert classify(free()) == "free"
    assert classify(state("claimed")) == "claimed"
    assert classify(state(None)) == "claimed"
    assert classify(state("registered")) == "registered"
    assert classify(state("fastq_staged")) == "fastq_staged"
    for st in ("aligned", "called", "qc_passed"):
        assert classify(state(st)) == "running"
    assert classify(state("summarized")) == "summarized"
    assert classify(state("summarized", merged=True)) == "merged"
    assert classify(state("merged")) == "merged"
    assert classify(state("released")) == "free"
    assert classify(state("excluded")) == "excluded"
