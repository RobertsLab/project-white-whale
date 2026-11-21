"""Fake Nextflow logs and study states for the autopilot tests.

The log lines copy the shape of real .nextflow.log lines (Nextflow 26.04) from
the ROBERTSLAB2013 run on Hete, trimmed to what the parser reads.
"""
from __future__ import annotations

from whale.autopilot.model import StudyState
from whale.autopilot.nflog import parse_log

SESSION = "f97ff558-1da0-4f7a-a026-2f37f8668acc"


def nf_log(*, command="nextflow -bg run workflows/methylation -c config/X.config -profile desktop,docker",
           run_name="maniac_pasteur", workdir="/Volumes/ugreen-2TB/whale/work", cached=(),
           submitted=(), failure=None, complete=True, body="") -> str:
    """A .nextflow.log. cached: (hash, process, tag); failure: dict(process, tag, exit, stderr, workdir)."""
    lines = [
        f"Nov-18 21:21:12.641 [main] DEBUG nextflow.cli.Launcher - $> {command}",
        f"Nov-18 21:21:13.686 [main] DEBUG nextflow.Session - Session UUID: {SESSION}",
        f"Nov-18 21:21:13.686 [main] DEBUG nextflow.Session - Run name: {run_name}",
        f"Nov-18 21:21:13.744 [main] DEBUG nextflow.Session - Work-dir: {workdir} [Mac OS X]",
    ]
    for h, proc, tag in cached:
        lines.append(f"Nov-18 21:21:15.205 [Actor Thread 5] INFO  nextflow.processor.TaskProcessor - "
                     f"[{h}] Cached process > {proc} ({tag})")
    for h, proc, tag in submitted:
        lines.append(f"Nov-18 21:22:15.205 [Task submitter] INFO  nextflow.Session - [{h}] Submitted process > {proc} ({tag})")
    if body:
        lines.append(body)
    if failure:
        p, t, code = failure["process"], failure.get("tag", "S1"), failure.get("exit", 1)
        lines += [
            f"Nov-19 06:37:14.724 [TaskFinalizer-5] ERROR nextflow.processor.TaskProcessor - Error executing process > '{p} ({t})'",
            "",
            "Caused by:",
            f"  Process `{p} ({t})` terminated with an error exit status ({code})",
            "",
            "Command exit status:",
            f"  {code}",
            "",
            "Command error:",
            f"  {failure.get('stderr', 'something went wrong')}",
            "",
            "Work dir:",
            f"  {failure.get('workdir', workdir + '/d5/5c1956a18a48988ccb799de07be7ee')}",
            "",
            f"Nov-19 06:37:14.754 [TaskFinalizer-5] DEBUG nextflow.Session - Session aborted -- Cause: Process "
            f"`{p} ({t})` terminated with an error exit status ({code})",
        ]
        stats = "succeededCount=32; failedCount=1; ignoredCount=0"
    else:
        stats = "succeededCount=40; failedCount=0; ignoredCount=0"
    if complete:
        lines += [
            f"Nov-19 06:37:15.059 [main] DEBUG n.trace.WorkflowStatsObserver - Workflow completed > WorkflowStats[{stats}; cachedCount=45]",
            "Nov-19 06:37:15.319 [main] DEBUG nextflow.script.ScriptRunner - > Execution complete -- Goodbye",
        ]
    return "\n".join(lines) + "\n"


def run_of(text: str, status: str | None = None, **kw):
    ev = parse_log(text, source="/checkout/.nextflow.log")
    if status:
        ev.status = status
    for k, v in kw.items():
        setattr(ev, k, v)
    return ev


CONFIG = """params {
    study_id     = 'TESTER2099_X_WGBS'
    chemistry    = 'wgbs'   // assay_chemistry.yaml term
    samplesheet  = 'data/studies/TESTER2099_X_WGBS/samplesheet.tsv'
}
"""


def state(stage=None, *, key="TESTER2099_X_WGBS", accession="PRJNA999901", machine="Hete", profile="desktop",
          branch=True, chemistry="wgbs", config=CONFIG, **kw) -> StudyState:
    status = None
    if stage is not None:
        status = {"study_id": key, "stage": stage, "machine": machine, "profile": profile,
                  "history": [{"stage": stage, "note": kw.pop("note", None)}]}
    s = StudyState(key=key, accession=accession, study_id=key,
                   branch=f"study/{accession}/{key}" if branch else None,
                   status=status, study_yaml={"study_id": key, "chemistry": chemistry},
                   config_text=config, **kw)
    return s


def qc_rows(n=12, *, mapping=55.0, noncpg=0.3, passed=None):
    passed = n if passed is None else passed
    rows = []
    for i in range(n):
        rows.append({"sample_id": f"S{i}", "mapping_pct": str(mapping + (i % 3) * 0.5),
                     "non_cpg_meth_pct": str(noncpg + (i % 3) * 0.05),
                     "conversion_pct": str(100 - noncpg), "qc_pass": "true" if i < passed else "false"})
    return rows


def summary(n=12, passed=None, *, snp_mask=True, assembly="GCF_963853765.1", on_host=True):
    passed = n if passed is None else passed
    return {"libraries": n, "passed": passed, "genome_assembly": assembly,
            "snp_mask": {"mask_id": "m1"} if snp_mask else None,
            "failed": [{"sample_id": f"S{i}", "reasons": "conversion 93% < 98.0%"} for i in range(passed, n)],
            "methylomes": {"on_data_host": on_host}}
