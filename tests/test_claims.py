import pytest
import yaml
from conftest import git

from whale.claims import (
    ClaimError,
    advance_status,
    branch_name,
    find_conflicts,
    list_claims,
    new_status,
    parse_branch,
    release_study,
    render_study_yaml,
    set_status,
    start_study,
)
from whale.scope import check_scope
from whale.validate import validate_status, validate_study

# Fictional studies: the claim tests seed their repos from the real registry/,
# so a real study here would break them as soon as it was registered.
ACC, SID = "PRJNA999901", "TESTER2099_CLAIMS_WGBS"
PART_ACC, PART_SID = "PRJEB999902", "TESTER2099_PARTS_WGBS"


# --- naming -------------------------------------------------------------------

def test_branch_round_trip():
    assert branch_name(ACC, SID) == f"study/{ACC}/{SID}"
    assert parse_branch(branch_name(ACC, SID)) == (ACC, SID, None)
    b = branch_name(PART_ACC, PART_SID, "gill")
    assert b == f"study/{PART_ACC}/{PART_SID}.part-gill"
    assert parse_branch(b) == (PART_ACC, PART_SID, "gill")
    assert parse_branch("main") is None


@pytest.mark.parametrize("acc,sid,part", [("10.1038/x", SID, None), ("SRR123", SID, None),
                                          (ACC, "bad id", None), (ACC, SID, "a/b")])
def test_branch_name_rejects(acc, sid, part):
    with pytest.raises(ClaimError):
        branch_name(acc, sid, part)


def test_conflicts_whole_study():
    assert find_conflicts(ACC, SID, None, {}, []) == []
    assert find_conflicts(ACC, SID, None, {SID: {ACC}}, [])
    assert find_conflicts(ACC, SID, None, {"OTHER2020_X": {ACC}}, [])
    assert find_conflicts(ACC, "NEW2020_X", None, {}, [f"study/{ACC}/{SID}"])


def test_conflicts_parts():
    reg = {SID: {ACC}}
    assert find_conflicts(ACC, SID, "gill", {}, [])                       # not registered
    assert find_conflicts(ACC, SID, "gill", reg, [f"study/{ACC}/{SID}"])  # still claimed whole
    assert find_conflicts(ACC, SID, "gill", reg, [f"study/{ACC}/{SID}.part-gill"])
    assert find_conflicts(ACC, SID, "gill", reg, [f"study/{ACC}/{SID}.part-mantle"]) == []


def test_rendered_yaml_validates(repo_root):
    template = (repo_root / "registry/studies/_TEMPLATE.yaml").read_text()
    text = render_study_yaml(template, study_id=SID, accession=ACC, registered_by="t", today="2026-10-23")
    study = yaml.safe_load(text)
    assert study["study_id"] == SID and study["accessions"]["bioproject"] == ACC
    assert validate_study(study, repo_root)[0] == []


def test_status_moves_forward_only(repo_root):
    s = new_status(study_id=SID, accession=ACC, part=None, profile="raven", machine="raven",
                   claimed_by="t", today="2026-10-23")
    assert validate_status(s, repo_root) == []
    s2 = advance_status(s, "aligned", repo=repo_root, machine="raven", today="2026-10-24")
    assert s2["stage"] == "aligned" and len(s2["history"]) == 2 and s["stage"] == "claimed"
    with pytest.raises(ClaimError):
        advance_status(s2, "fastq_staged", repo=repo_root, machine="raven", today="2026-10-24")
    back = advance_status(s2, "fastq_staged", repo=repo_root, machine="raven", today="2026-10-24", force=True)
    assert back["stage"] == "fastq_staged"
    ex = advance_status(s2, "excluded", repo=repo_root, machine="raven", today="2026-10-24")
    with pytest.raises(ClaimError):
        advance_status(ex, "called", repo=repo_root, machine="raven", today="2026-10-25")
    moved = advance_status(s2, "called", repo=repo_root, machine="n3000", today="2026-10-25", profile="klone")
    assert moved["profile"] == "klone" and moved["machine"] == "n3000"
    assert validate_status(moved, repo_root) == []


# --- against a real (local) remote -------------------------------------------

def test_claim_lifecycle(claim_env):
    clone = claim_env["clone"]
    wt = claim_env["tmp"] / "wt"
    res = start_study(ACC, SID, repo=clone, profile="raven", worktree_dir=wt,
                      registered_by="Tester", machine="raven", open_pr=False)
    assert res["pushed"]
    assert (wt / f"registry/studies/{SID}.yaml").exists()
    assert (wt / f"config/{SID}.config").exists()
    status = yaml.safe_load((wt / f"data/studies/{SID}/status.yaml").read_text())
    assert status["stage"] == "claimed" and status["profile"] == "raven"

    # A second machine sees the claim and cannot take the study.
    other = claim_env["tmp"] / "machine-b"
    git("clone", "--quiet", str(claim_env["remote"]), str(other), cwd=claim_env["tmp"])
    with pytest.raises(ClaimError, match="already claimed"):
        start_study(ACC, SID, repo=other, profile="klone", worktree_dir=claim_env["tmp"] / "wt-b",
                    registered_by="Other", machine="klone", open_pr=False)

    # Progress is pushed and visible from the other machine.
    set_status("fastq_staged", repo=wt, machine="raven", note="140 GB, MD5 ok")
    rows = list_claims(other)
    assert [(r["study_id"], r["stage"], r["machine"]) for r in rows] == [(SID, "fastq_staged", "raven")]

    # The branch passes its own scope check.
    git("fetch", "--quiet", "origin", cwd=wt)
    changed = git("diff", "--name-only", "origin/main...HEAD", cwd=wt).split()
    study = yaml.safe_load((wt / f"registry/studies/{SID}.yaml").read_text())
    assert check_scope(res["branch"], changed, study, {}) == ([], [])

    # Releasing frees it for the other machine.
    release_study(res["branch"], repo=clone, reason="moving to klone", close_pr=False)
    assert list_claims(other) == []
    res_b = start_study(ACC, SID, repo=other, profile="klone", worktree_dir=claim_env["tmp"] / "wt-b",
                        registered_by="Other", machine="klone", open_pr=False)
    assert res_b["pushed"]


def test_no_push_claims_nothing(claim_env):
    clone = claim_env["clone"]
    res = start_study(ACC, SID, repo=clone, profile="desktop", worktree_dir=claim_env["tmp"] / "wt",
                      registered_by="T", machine="mac", push=False)
    assert not res["pushed"]
    assert list_claims(clone) == []
    # ...but this machine will not claim it twice.
    with pytest.raises(ClaimError, match="already claimed"):
        start_study(ACC, SID, repo=clone, profile="desktop", worktree_dir=claim_env["tmp"] / "wt2",
                    registered_by="T", machine="mac", push=False)


def test_unknown_profile(claim_env):
    with pytest.raises(ClaimError, match="unknown profile"):
        start_study(ACC, SID, repo=claim_env["clone"], profile="laptop",
                    worktree_dir=claim_env["tmp"] / "wt", registered_by="T", open_pr=False)


def test_part_claim_after_registration(claim_env, filled_study):
    clone = claim_env["clone"]
    acc, sid = PART_ACC, PART_SID
    with pytest.raises(ClaimError, match="registered on main"):
        start_study(acc, sid, repo=clone, profile="klone", part="gill",
                    worktree_dir=claim_env["tmp"] / "wt-gill", registered_by="T", open_pr=False)

    study = dict(filled_study, study_id=sid)
    study["accessions"] = dict(filled_study["accessions"], bioproject=acc)
    (clone / f"registry/studies/{sid}.yaml").write_text(yaml.safe_dump(study, sort_keys=False))
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", "register", cwd=clone)
    git("push", "--quiet", "origin", "main", cwd=clone)

    res = start_study(acc, sid, repo=clone, profile="klone", part="gill",
                      worktree_dir=claim_env["tmp"] / "wt-gill", registered_by="T",
                      machine="klone", open_pr=False)
    assert res["branch"].endswith(".part-gill")
    assert sorted(res["files"]) == [f"config/{sid}.part-gill.config",
                                    f"data/studies/{sid}/parts/gill/status.yaml"]
    start_study(acc, sid, repo=clone, profile="raven", part="mantle",
                worktree_dir=claim_env["tmp"] / "wt-mantle", registered_by="T",
                machine="raven", open_pr=False)
    rows = {r["part"]: r for r in list_claims(clone)}
    assert set(rows) == {"gill", "mantle"} and rows["mantle"]["profile"] == "raven"
