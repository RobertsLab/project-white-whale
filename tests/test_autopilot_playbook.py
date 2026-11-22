"""Every playbook entry and library check, on fake logs and reports."""
from pathlib import Path

from autopilot_fakes import CONFIG, SESSION, nf_log, qc_rows, run_of, state, summary

from whale.autopilot import playbook
from whale.autopilot.library import check_nondirectional, check_read2_tail, readiness
from whale.autopilot.machines import Machine
from whale.autopilot.nflog import parse_log, read_run

HETE = Machine(name="Hete", profile="desktop", env={"WHALE_RAW": "/w/raw", "WHALE_REFERENCE": "/w/ref",
                                                     "WHALE_METHYLOMES": "/w/meth", "WHALE_WORKDIR": "/w/work"},
               extra_configs=["/w/hete-limits.config"], gannet_key=False)
RAVEN = Machine(name="raven", profile="raven", gh=False, gannet_key=True)
KLONE = Machine(name="klone", profile="klone", launch="sbatch", head_job="--account=coenv --partition=compute",
                slurm_accounts=[{"partition": "compute", "account": "coenv"},
                                {"partition": "cpu-g2-mem2x", "account": "srlab"}])
RESUME = f"nextflow run workflows/methylation -c config/X.config -profile desktop,docker -resume {SESSION}"


def ids(s, m):
    return [x.entry for x in playbook.match(s, m)]


# --- log parsing ---------------------------------------------------------------

def test_parse_failed_log():
    ev = parse_log(nf_log(failure={"process": "BISMARK_DEDUPLICATE", "tag": "SRR1", "exit": 137,
                                   "stderr": "Killed"}))
    assert ev.status == "failed" and ev.session == SESSION and ev.run_name == "maniac_pasteur"
    t = ev.failed[0]
    assert (t.process, t.tag, t.exit_status) == ("BISMARK_DEDUPLICATE", "SRR1", 137)
    assert t.workdir.endswith("d5/5c1956a18a48988ccb799de07be7ee") and "Killed" in t.stderr


def test_parse_ok_and_resumed_log():
    ev = parse_log(nf_log(command=RESUME, cached=[("40/7189a7", "TRIM_GALORE", "S1")]))
    assert ev.status == "ok" and ev.resumed
    assert ev.cached[0].hash == "40/7189a7" and ev.cached[0].tag == "S1"


def test_read_run_live_stopped_and_history(tmp_path):
    (tmp_path / ".nextflow.log").write_text(nf_log(complete=False))
    (tmp_path / ".nextflow.pid").write_text("4242\n")
    assert read_run(tmp_path, alive=lambda pid: True).status == "live"
    assert read_run(tmp_path, alive=lambda pid: False).status == "stopped"
    (tmp_path / ".nextflow").mkdir()
    (tmp_path / ".nextflow/history").write_text(
        f"2026-10-25 09:02:22\t4m\tearlier\tERR\tabc\t{SESSION}\tnextflow run x\n"
        f"2026-11-18 21:21:13\t9h\tmaniac_pasteur\tERR\tabc\t{SESSION}\tnextflow run x -resume\n")
    ev = read_run(tmp_path, alive=lambda pid: True)
    assert ev.status == "failed" and ev.previous_runs == 1
    assert read_run(tmp_path / "missing") is None


def test_read_run_task_files(tmp_path):
    wd = tmp_path / "work" / "ab" / "cdef"
    wd.mkdir(parents=True)
    (wd / ".exitcode").write_text("137")
    (wd / ".command.err").write_text("deduplicate_bismark: Killed\n")
    log = nf_log(failure={"process": "P", "tag": "S1", "exit": 137, "workdir": str(wd), "stderr": ""})
    log = log.replace("Command exit status:\n  137\n", "").replace("Command error:\n  \n", "")
    (tmp_path / ".nextflow.log").write_text(log)
    t = read_run(tmp_path).failed[0]
    assert t.exit_status == 137 and "Killed" in t.stderr


# --- entries -------------------------------------------------------------------

def test_swap_mates_resume():
    cfg = CONFIG.replace("}", "    swap_mates = ['S1', // read 2 first\n        'S2']\n}")
    log = nf_log(command=RESUME, cached=[("40/7189a7", "TRIM_GALORE", "S1"), ("41/000000", "TRIM_GALORE", "S9"),
                                         ("42/aaaaaa", "BISMARK_ALIGN", "S2")])
    s = state("fastq_staged", config=cfg, run=run_of(log))
    m = playbook.swap_mates_resume(s, HETE)
    assert m and m.hazard            # a successful run with stale tasks must not be recorded
    assert "mv /Volumes/ugreen-2TB/whale/work/40/7189a7*" in m.remedy and "42/aaaaaa" in m.remedy
    assert "41/000000" not in m.remedy
    assert ids(s, HETE)[0] == "swap_mates_resume"
    # Not resumed, or nothing swapped: no match.
    assert playbook.swap_mates_resume(state("fastq_staged", config=cfg, run=run_of(log.replace(f" -resume {SESSION}", ""))), HETE) is None
    assert playbook.swap_mates_resume(state("fastq_staged", run=run_of(log)), HETE) is None


def test_resume_env_cache_miss():
    log = nf_log(command=RESUME, workdir="/checkout/work", failure={"process": "FASTQC", "exit": 1})
    s = state("fastq_staged", run=run_of(log), worktree="/checkout")
    m = playbook.resume_env_cache_miss(s, HETE)
    assert m and "WHALE_WORKDIR=/w/work" in m.remedy and f"-resume {SESSION}" in m.remedy
    assert any("inside the checkout" in e for e in m.evidence)
    # Work dir is this machine's: no evidence of a miss.
    ok = state("fastq_staged", run=run_of(log.replace("/checkout/work", "/w/work")), worktree="/checkout")
    assert playbook.resume_env_cache_miss(ok, HETE) is None
    # A resume that cached tasks is fine.
    cached = run_of(nf_log(command=RESUME, cached=[("40/7189a7", "FASTQC", "S1")]))
    assert playbook.resume_env_cache_miss(state("fastq_staged", run=cached), HETE) is None


def test_klone_partition_account_and_scontrol():
    log = nf_log(failure={"process": "BISMARK_ALIGN", "exit": 1},
                 body="sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified")
    s = state("fastq_staged", machine="klone", profile="klone", run=run_of(log))
    m = playbook.klone_partition(s, KLONE)
    assert m and "compute=coenv" in m.remedy and "cpu-g2-mem2x=srlab" in m.remedy and "never move" in m.remedy
    moved = run_of(nf_log(body="scontrol update JobId=4242 Partition=cpu-g2"))
    assert any("scontrol" in e for e in playbook.klone_partition(state("fastq_staged", run=moved), KLONE).evidence)
    preempt = run_of(nf_log(body="\n".join(["task terminated for an unknown reason -- Likely it has been terminated by the external system"] * 3)))
    assert playbook.klone_partition(state("fastq_staged", run=preempt), KLONE)
    assert playbook.klone_partition(s, HETE) is None   # not klone


def test_chemistry_null():
    null_cfg = CONFIG.replace("'wgbs'", "null")
    m = playbook.chemistry_null(state("registered", config=null_cfg), HETE)
    assert m and "chemistry = 'wgbs'" in m.remedy
    m = playbook.chemistry_null(state("registered", config=CONFIG.replace("'wgbs'", "'emseq'")), HETE)
    assert m and "emseq" in m.evidence[0]
    assert playbook.chemistry_null(state("registered"), HETE) is None
    m = playbook.chemistry_null(state("registered", config=null_cfg, chemistry=None), HETE)
    assert m and "study YAML first" in m.remedy


def test_raven_no_gh():
    s = state("summarized", machine="raven", profile="raven", run=run_of(nf_log(body="bash: gh: command not found")))
    m = playbook.raven_no_gh(s, RAVEN)
    assert m and "REST" in m.remedy
    assert playbook.raven_no_gh(state("summarized", run=run_of(nf_log())), RAVEN) is None


def test_gannet_upload():
    s = state("qc_passed", machine="klone", summary=summary(on_host=False))
    m = playbook.gannet_upload(s, KLONE)
    assert m and "record-run --results $WHALE_METHYLOMES/TESTER2099_X_WGBS" in m.remedy and "raven" in m.remedy
    log = run_of(nf_log(body="methylomes not on the data host: Permission denied (publickey)"))
    assert playbook.gannet_upload(state("called", run=log), HETE)
    assert playbook.gannet_upload(state("summarized", summary=summary()), HETE) is None


def test_transient_network():
    s = state("registered", run=run_of("", extra_text="curl: (56) Recv failure: Connection reset by peer\n" * 2))
    s.run.status = "unknown"
    m = playbook.transient_network(s, HETE)
    assert m and "--source mirror" in m.remedy and "mirror-fastq --accession PRJNA999901" in m.remedy
    assert "~20 min" in m.remedy      # backoff grows with the number of errors
    assert playbook.transient_network(state("registered", run=run_of(nf_log())), HETE) is None


def test_oom_killed(tmp_path):
    log = nf_log(failure={"process": "BISMARK_DEDUPLICATE", "tag": "SRR1", "exit": 137})
    limits = tmp_path / "limits.config"
    limits.write_text("process {\n  withName: BISMARK_ALIGN {\n    maxForks = 1\n  }\n}\n")
    hete = Machine(**{**HETE.__dict__, "limits_config": str(limits)})
    m = playbook.oom_killed(state("fastq_staged", run=run_of(log)), hete)
    assert m and "withName: 'BISMARK_DEDUPLICATE' { maxForks = 1 }" in m.remedy and str(limits) in m.remedy
    assert "-c /w/hete-limits.config" in m.remedy and f"-resume {SESSION}" in m.remedy
    limits.write_text("process {\n  withName: 'BISMARK_DEDUPLICATE' { maxForks = 1 }\n}\n")
    m = playbook.oom_killed(state("fastq_staged", run=run_of(log)), hete)
    assert "already capped" in m.remedy
    other = run_of(nf_log(failure={"process": "P", "exit": 2}))
    assert playbook.oom_killed(state("fastq_staged", run=other), HETE) is None


def test_head_stopped():
    s = state("fastq_staged", run=run_of(nf_log(complete=False), status="stopped"))
    m = playbook.head_stopped(s, HETE)
    assert m and "-resume" in m.remedy and "export WHALE_RAW=/w/raw" in m.remedy


def test_no_match_gives_diagnosis_evidence():
    s = state("fastq_staged", run=run_of(nf_log(failure={"process": "QC_GATE", "exit": 2,
                                                         "stderr": "KeyError: 'cpgs_5x'"})))
    assert playbook.match(s, HETE) == []
    ev = playbook.diagnosis_evidence(s)
    assert any("QC_GATE" in e and "exit 2" in e for e in ev) and any("KeyError" in e for e in ev)


def test_relaunch_on_klone_is_sbatch():
    s = state("fastq_staged", machine="klone", profile="klone", run=run_of(nf_log()), run_config_text="params {}")
    cmd = playbook._relaunch(s, KLONE)
    assert cmd.startswith("sbatch --account=coenv --partition=compute --job-name=whale-head-TESTER2099_X_WGBS --wrap")
    assert "-c data/studies/TESTER2099_X_WGBS/klone.run.config" in cmd and "-profile klone" in cmd


def test_klone_head_job_runs_in_the_worktree_with_its_tools():
    # Like the hand-written run-DECICOMP2025.sbatch: worktree, conda env on PATH, no TMPDIR, a log.
    from whale.autopilot.machines import load_machines
    from whale.autopilot.slurm import Queue, SlurmJob, assign
    klone = load_machines(Path(__file__).resolve().parents[1])["klone"]
    s = state("fastq_staged", key="TESTER2099_X_WGBS", machine="klone", profile="klone", run_config_text="params {}")
    cmd = playbook._relaunch(s, klone)
    assert "--job-name=whale-head-TESTER2099_X_WGBS" in cmd
    assert "--chdir=/gscratch/srlab/$USER/github/white-whale-TESTER2099_X_WGBS" in cmd
    assert "--output=/gscratch/scrubbed/$USER/whale/logs/whale-head-TESTER2099_X_WGBS-%j.log" in cmd
    wrap = cmd.split("--wrap '", 1)[1]
    steps = wrap.rstrip("'").split("; ")
    assert steps[0] == "export PATH=/mmfs1/gscratch/srlab/$USER/miniforge3/envs/whale/bin:$PATH"
    assert steps[1] == "unset TMPDIR" and "WHALE_WORKDIR=/gscratch/scrubbed/$USER/whale/work" in steps[2]
    assert steps[3].startswith("nextflow run workflows/methylation") and steps[3].endswith("-resume")
    # The queue reader attributes the job to the study by its name and working directory.
    job = SlurmJob("1", "whale-head-TESTER2099_X_WGBS", "compute", "compute-coenv", "RUNNING", "0:01",
                   "/mmfs1/gscratch/srlab/sr320/github/white-whale-TESTER2099_X_WGBS")
    assert list(assign(["TESTER2099_X_WGBS", "TESTER2099_Y_WGBS"], Queue("klone", [job]))) == ["TESTER2099_X_WGBS"]


def test_klone_part_worktree_and_job_name():
    m = Machine(name="klone", profile="klone", launch="sbatch", checkout_root="/g/github/")
    assert m.checkout("X_WGBS", "gill") == "/g/github/white-whale-X_WGBS-gill"
    assert "--job-name=whale-fastq-X_WGBS-gill --chdir=/g/github/white-whale-X_WGBS-gill" in \
        m.sbatch("whale fetch-fastq --mark-staged", study_id="X_WGBS", part="gill", kind="fastq")
    assert Machine(name="raven", profile="raven").checkout("X_WGBS") is None


MASK = "GCF_963853765.1.snpmask.v1"


def test_launch_exports_the_snp_mask_of_the_study_assembly():
    from whale.autopilot.machines import load_machines
    machines = load_machines(Path(__file__).resolve().parents[1])
    for name, folder in [("klone", f"/gscratch/srlab/$USER/whale/snp_mask/{MASK}"),
                         ("raven", f"/home/shared/8TB_HDD_03/$USER/whale/reference/GCF_963853765.1/snp_mask/{MASK}")]:
        s = state("fastq_staged", machine=name, profile=name, genome_assembly="GCF_963853765.1", snp_mask_id=MASK)
        assert f"WHALE_SNP_MASK={folder}" in playbook._relaunch(s, machines[name]), name
    # Hete has no masks: nothing exported (a missing folder would fail the run's checkIfExists).
    s = state("fastq_staged", genome_assembly="GCF_963853765.1", snp_mask_id=MASK)
    assert "WHALE_SNP_MASK" not in playbook._relaunch(s, machines["Hete"])


def test_snp_mask_root_with_assembly_placeholder():
    m = Machine(name="raven", profile="raven", snp_mask_root="/r/{assembly}/snp_mask/")
    assert m.snp_mask(MASK, "GCF_963853765.1") == f"/r/GCF_963853765.1/snp_mask/{MASK}"
    assert m.snp_mask(MASK, None) is None and m.snp_mask(None, "GCF_963853765.1") is None


def test_newest_mask_per_assembly_and_study_assembly():
    from whale.autopilot.gather import attach_masks, newest_masks
    masks = newest_masks({"masks": [{"mask_id": "A.v1", "genome_assembly": "A"},
                                    {"mask_id": "B.v1", "genome_assembly": "B"},
                                    {"mask_id": "A.v2", "genome_assembly": "A"}]})
    assert masks == {"A": "A.v2", "B": "B.v1"}
    default = state("fastq_staged", key="TESTER2099_D_WGBS")
    own = state("fastq_staged", key="TESTER2099_O_WGBS", config=CONFIG + "\nparams { genome_assembly = 'B' }\n")
    none = state("fastq_staged", key="TESTER2099_N_WGBS", config=CONFIG + "\nparams { genome_assembly = 'C' }\n")
    attach_masks([default, own, none], masks, "A")
    assert (default.genome_assembly, default.snp_mask_id) == ("A", "A.v2")
    assert (own.genome_assembly, own.snp_mask_id) == ("B", "B.v1")
    assert (none.genome_assembly, none.snp_mask_id) == ("C", None)


# --- library type and readiness ------------------------------------------------

def test_nondirectional_detected_and_proposed():
    s = state("summarized", qc=qc_rows(mapping=29.0))
    f = check_nondirectional(s)
    assert not f.ok and "--non_directional" in f.proposal
    assert check_nondirectional(state("summarized", qc=qc_rows(mapping=55.0))).ok
    done = state("summarized", qc=qc_rows(mapping=29.0),
                 run_manifest={"params": {"bismark_align_args": "--non_directional"}})
    assert check_nondirectional(done).ok
    # Config fixed since the recorded run: propose a re-run, not the same change again.
    fixed = state("merged", qc=qc_rows(mapping=29.0), run_manifest={"params": {"bismark_align_args": ""}},
                  config=CONFIG.replace("}", "    bismark_align_args = '--non_directional'\n}"))
    f = check_nondirectional(fixed)
    assert not f.ok and "re-run" in f.proposal


def test_read2_tail_detected_and_proposed():
    f = check_read2_tail(state("summarized", qc=qc_rows(noncpg=2.9)))
    assert not f.ok and "--clip_r2 18" in f.proposal and "98%" in f.proposal
    assert check_read2_tail(state("summarized", qc=qc_rows(noncpg=0.3))).ok
    # Poor conversion in some libraries (spread out) is not a tail.
    rows = qc_rows(10, noncpg=0.5) + qc_rows(10, noncpg=6.0)
    assert check_read2_tail(state("summarized", qc=rows)).ok
    clipped = state("summarized", qc=qc_rows(noncpg=2.9), run_manifest={"params": {"trim_args": "--clip_r2 18"}})
    assert check_read2_tail(clipped).ok


def test_readiness():
    good = readiness(state("summarized", qc=qc_rows(), summary=summary()), "GCF_963853765.1")
    assert all(f.ok for f in good)
    bad = {f.check: f for f in readiness(state("summarized", qc=qc_rows(48, passed=13),
                                               summary=summary(48, 13, snp_mask=False, assembly="GCF_X")),
                                         "GCF_963853765.1")}
    assert not bad["pass_rate"].ok and "13 of 48" in bad["pass_rate"].detail and "conversion" in bad["pass_rate"].detail
    assert not bad["snp_mask"].ok and not bad["assembly"].ok
