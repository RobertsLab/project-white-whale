"""Reading a Slurm queue for the autopilot (src/whale/autopilot/slurm.py). No network: the runner is faked."""
import subprocess

from autopilot_fakes import state

from whale.autopilot.machines import Machine
from whale.autopilot.slurm import Queue, assign, parse, read_queue
from whale.autopilot.tick import attach_queues

# Shape of real `squeue -u $USER -h -o "%i|%j|%P|%a|%T|%M|%Z"` output on klone (2026-11-21),
# with the conda warning klone's login shell prints first.
SQUEUE = """WARNING conda.cli.main_config:_set_key(451): Key auto_activate_base is an alias of auto_activate
41261110|$HOME/bin/hyakalloc_snapshot.sh|ckpt|ckpt-srlab|PENDING|0:00|/mmfs1/home/sr320
43450651|nf-BISMARK_ALIGN_(ERR15996566)|compute|compute-coenv|PENDING|0:00|/mmfs1/gscratch/scrubbed/sr320/whale/work/be/ed5fac
43445189|whale-head-DECICOMP2025|compute|compute-coenv|RUNNING|10:17:55|/mmfs1/gscratch/srlab/sr320/github/white-whale-DECICOMP2025_POMS_AGE_EMSEQ
43451043|nf-BISMARK_ALIGN_(ERR15996567)|cpu-g2-mem2x|cpu-g2-mem2x-srlab|RUNNING|1:42:38|/mmfs1/gscratch/scrubbed/sr320/whale/work/49/4c7d64
"""


def _runner(stdout="", returncode=0, stderr="", calls=None):
    def run(cmd):
        if calls is not None:
            calls.append(cmd)
        return subprocess.CompletedProcess(cmd, returncode, stdout, stderr)
    return run


def test_parse_skips_banners():
    jobs = parse(SQUEUE)
    assert [j.job_id for j in jobs] == ["41261110", "43450651", "43445189", "43451043"]
    head = jobs[2]
    assert (head.name, head.account, head.state, head.is_head) == \
        ("whale-head-DECICOMP2025", "compute-coenv", "RUNNING", True)
    assert not jobs[1].is_head


def test_read_queue_over_ssh_never_prompts():
    calls = []
    q = read_queue("klone", "klone.hyak.uw.edu", here="Hete", runner=_runner(SQUEUE, calls=calls))
    assert q.ok and len(q.jobs) == 4
    cmd = calls[0]
    assert cmd[:3] == ["ssh", "-o", "BatchMode=yes"] and "klone.hyak.uw.edu" in cmd
    assert cmd[-1].startswith("squeue -u $USER")       # $USER expands on klone, not here


def test_read_queue_locally_on_the_machine():
    calls = []
    read_queue("klone", "klone.hyak.uw.edu", here="klone", runner=_runner(SQUEUE, calls=calls))
    assert calls[0][0] == "bash" and "squeue -u $USER" in calls[0][-1]


def test_read_queue_reports_why_it_failed():
    q = read_queue("klone", "klone.hyak.uw.edu", here="Hete",
                   runner=_runner(returncode=255, stderr="sr320@klone.hyak.uw.edu: Permission denied (keyboard-interactive)."))
    assert not q.ok and "Permission denied" in q.error and q.jobs == []
    assert "no ssh_host" in read_queue("klone", None, here="Hete").error


def test_assign_by_workdir_name_and_unique_head_prefix():
    q = Queue("klone", parse(SQUEUE))
    got = assign(["DECICOMP2025_POMS_AGE_EMSEQ", "PESTO2022_METHYLSEQ_WGBS"], q)
    assert [j.job_id for j in got["DECICOMP2025_POMS_AGE_EMSEQ"]] == ["43445189"]
    assert "PESTO2022_METHYLSEQ_WGBS" not in got
    # A prefix shared by two studies is ambiguous unless the workdir names one of them.
    q2 = Queue("klone", parse("1|whale-head-OUC2024|compute|compute-coenv|RUNNING|1:00|/home/sr320\n"))
    assert assign(["OUC2024_SEX_REVERSAL_WGBS", "OUC2024_GROWTH_STRAINS_WGBS"], q2) == {}
    assert list(assign(["OUC2024_SEX_REVERSAL_WGBS"], q2)) == ["OUC2024_SEX_REVERSAL_WGBS"]


def test_attach_queues_reads_each_sbatch_machine_once():
    machines = {"klone": Machine(name="klone", profile="klone", launch="sbatch", ssh_host="klone.hyak.uw.edu"),
                "Hete": Machine(name="Hete", profile="desktop")}
    deci = state("fastq_staged", key="DECICOMP2025_POMS_AGE_EMSEQ", machine="klone", profile="klone")
    other = state("fastq_staged", key="TESTER2099_Y_WGBS", machine="klone", profile="klone")
    local = state("fastq_staged", key="TESTER2099_Z_WGBS")
    calls = []
    warnings = attach_queues([deci, other, local], machines, "Hete", _runner(SQUEUE, calls=calls))
    assert warnings == [] and len(calls) == 1
    assert [j.job_id for j in deci.head_jobs] == ["43445189"] and other.slurm_checked and other.slurm_jobs == []
    assert not local.slurm_checked
    warnings = attach_queues([deci], machines, "Hete", _runner(returncode=255, stderr="Connection timed out"))
    assert "could not read the klone queue from Hete: Connection timed out" in warnings[0]
    assert (deci.slurm_checked, deci.slurm_jobs, deci.slurm_error) == (False, [], "Connection timed out")
