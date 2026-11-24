"""Checks on a finished run's reports: library type and PR readiness.

Library-type problems found by auditing the first summarized studies
(docs/adding_a_study.md, "Check the library type"):

- non-directional libraries (Zymo Pico Methyl-Seq: read 1 ~9% C) aligned in
  directional mode map about half of what they should (ROBERTS2022: 30%);
- a synthetic tail at the start of read 2 reads as methylated C in every
  context: 2-4% non-CpG methylation in every library, which looks like poor
  conversion (WANG2023).

Both are detected from qc_summary.tsv (and FastQC / M-bias reports when they
are readable) and turned into a *proposed* config change. Nothing is applied.

The readiness check is what a person reads before calling a study PR ready:
pass rate, SNP mask, genome assembly, and the library checks above.
"""
from __future__ import annotations

import contextlib
import statistics
from dataclasses import dataclass, field

from whale.autopilot.model import StudyState
from whale.autopilot.playbook import config_string
from whale.claims import config_name

# Directional mode on a non-directional library maps about half; a normal
# directional WGBS library on this assembly maps 45-60%.
NONDIR_MAX_MEDIAN_MAPPING = 40.0
NONDIR_MAX_SPREAD = 15.0
# A synthetic read-2 tail: every library at 1.5-5% non-CpG methylation, tightly.
TAIL_NONCPG = (1.5, 5.0)
TAIL_MIN_FRACTION = 0.8
TAIL_MAX_SPREAD = 2.5
MIN_PASS_FRACTION = 0.75

BASE_RESOLUTION = {"wgbs", "emseq", "pbat", "rrbs", "bs_capture"}


@dataclass
class Finding:
    check: str
    ok: bool
    detail: str
    proposal: str | None = None
    evidence: list[str] = field(default_factory=list)


def _col(qc: list[dict], name: str) -> list[float]:
    out = []
    for r in qc:
        with contextlib.suppress(KeyError, TypeError, ValueError):
            out.append(float(r[name]))
    return out


def _params(study: StudyState, key: str) -> str:
    """A run parameter as the recorded run used it, else as the config sets it now."""
    run_params = (study.run_manifest or {}).get("params") or {}
    if key in run_params:
        return str(run_params.get(key) or "")
    return config_string(study.config_text, key) or ""


def check_nondirectional(study: StudyState) -> Finding | None:
    qc = study.qc or []
    chem = (study.study_yaml or {}).get("chemistry") or _params(study, "chemistry")
    if not qc or chem not in BASE_RESOLUTION or chem in ("pbat", "rrbs"):
        return None
    used = _params(study, "bismark_align_args")
    now = config_string(study.config_text, "bismark_align_args") or ""
    mapping = _col(qc, "mapping_pct")
    if not mapping:
        return None
    med, spread = statistics.median(mapping), max(mapping) - min(mapping)
    name = config_name(study.study_id or study.key, study.part)
    if "--non_directional" in used:
        return Finding("library_type", True, "aligned non-directionally")
    if "--non_directional" in now:
        return Finding("library_type", False,
                       f"{name} now sets --non_directional but the recorded run did not use it",
                       proposal="re-run with the current config (-resume from trimming) and record it again",
                       evidence=[f"run_manifest bismark_align_args = {used!r}"])
    if med < NONDIR_MAX_MEDIAN_MAPPING and spread <= NONDIR_MAX_SPREAD:
        return Finding(
            "library_type", False,
            f"every library maps low in directional mode (median {med:.1f}%, range {min(mapping):.1f}-"
            f"{max(mapping):.1f}%): looks non-directional (e.g. Zymo Pico)",
            proposal=f"confirm read 1 is ~9% C in FastQC per-base content, then add --non_directional to "
                     f"bismark_align_args in {name} and re-run",
            evidence=[f"median mapping_pct {med:.1f}% over {len(mapping)} libraries, directional mode"],
        )
    return Finding("library_type", True, f"median mapping {med:.1f}%")


def check_read2_tail(study: StudyState) -> Finding | None:
    qc = study.qc or []
    noncpg = _col(qc, "non_cpg_meth_pct")
    if not noncpg:
        return None
    used_trim = _params(study, "trim_args")
    now_trim = config_string(study.config_text, "trim_args") or ""
    lo, hi = TAIL_NONCPG
    inside = sum(lo <= v <= hi for v in noncpg)
    spread = max(noncpg) - min(noncpg)
    name = config_name(study.study_id or study.key, study.part)
    if "--clip_r2" in used_trim:
        return Finding("read2_tail", True, f"read 2 clipped ({used_trim})")
    if "--clip_r2" in now_trim:
        return Finding("read2_tail", False,
                       f"{name} now clips read 2 but the recorded run did not",
                       proposal="re-run with the current config and record it again",
                       evidence=[f"run_manifest trim_args = {used_trim!r}"])
    if inside >= TAIL_MIN_FRACTION * len(noncpg) and spread <= TAIL_MAX_SPREAD:
        return Finding(
            "read2_tail", False,
            f"{inside} of {len(noncpg)} libraries at {lo}-{hi}% non-CpG methylation (median "
            f"{statistics.median(noncpg):.2f}%): looks like a synthetic read-2 tail, not poor conversion",
            proposal=f"confirm in the read-2 M-bias plot (raised CHG/CHH at the 5' end), then set "
                     f"trim_args = '--clip_r2 18 --three_prime_clip_r1 18' in {name}, keep the 98% "
                     f"conversion gate, and re-run",
            evidence=[f"non_cpg_meth_pct range {min(noncpg):.2f}-{max(noncpg):.2f}%"],
        )
    return Finding("read2_tail", True, f"non-CpG median {statistics.median(noncpg):.2f}%")


def library_findings(study: StudyState) -> list[Finding]:
    return [f for f in (check_nondirectional(study), check_read2_tail(study)) if f]


def readiness(study: StudyState, main_assembly: str | None) -> list[Finding]:
    """The checks a person makes before calling a study PR ready."""
    out = []
    s = study.summary or {}
    n, passed = s.get("libraries") or len(study.qc or []), s.get("passed")
    if passed is None and study.qc:
        passed = sum(r.get("qc_pass") == "true" for r in study.qc)
    if n:
        frac = (passed or 0) / n
        reasons = sorted({(f.get("reasons") or "").split(" ")[0] for f in s.get("failed") or []} - {""})
        out.append(Finding("pass_rate", frac >= MIN_PASS_FRACTION,
                           f"{passed} of {n} libraries pass the QC gate"
                           + (f" (failures: {', '.join(reasons)})" if reasons and frac < 1 else "")))
    if "snp_mask" in s or study.summary is not None:
        out.append(Finding("snp_mask", bool(s.get("snp_mask")),
                           "SNP-masked" if s.get("snp_mask") else
                           "no SNP mask: the atlas leaves the study out",
                           proposal=None if s.get("snp_mask") else "re-run with WHALE_SNP_MASK set and -resume"))
    asm = s.get("genome_assembly")
    if asm and main_assembly:
        out.append(Finding("assembly", asm == main_assembly,
                           f"assembly {asm}" + ("" if asm == main_assembly else f", main uses {main_assembly}")))
    out += library_findings(study)
    return out
