from whale.scope import check_scope

ACC, SID = "PRJNA682817", "ROBERTS2022_PH_PLOIDY_WGBS"
BRANCH = f"study/{ACC}/{SID}"
STUDY = {"study_id": SID, "accessions": {"bioproject": ACC}}


def test_owned_files_pass():
    changed = [f"registry/studies/{SID}.yaml", f"config/{SID}.config",
               f"data/studies/{SID}/status.yaml", f"data/studies/{SID}/samplesheet.tsv",
               f"tests/test_{SID.lower()}.py"]
    assert check_scope(BRANCH, changed, STUDY, {}) == ([], [])


def test_shared_file_fails():
    errors, _ = check_scope(BRANCH, ["src/whale/cli.py", f"config/{SID}.config"], STUDY, {})
    assert errors and "src/whale/cli.py" in errors[0]


def test_whole_study_branch_may_not_touch_parts():
    errors, _ = check_scope(BRANCH, [f"data/studies/{SID}/parts/gill/status.yaml"], STUDY, {})
    assert errors


def test_part_branch_owns_only_its_part():
    branch = f"{BRANCH}.part-gill"
    ok = [f"data/studies/{SID}/parts/gill/status.yaml", f"config/{SID}.part-gill.config"]
    assert check_scope(branch, ok, STUDY, {SID: {ACC}}) == ([], [])
    errors, _ = check_scope(branch, [f"registry/studies/{SID}.yaml"], STUDY, {SID: {ACC}})
    assert errors
    errors, _ = check_scope(branch, [f"data/studies/{SID}/parts/mantle/status.yaml"], STUDY, {SID: {ACC}})
    assert errors


def test_yaml_identity_checks():
    assert check_scope(BRANCH, [], None, {})[0]
    assert check_scope(BRANCH, [], {**STUDY, "study_id": "OTHER2020_X"}, {})[0]
    assert check_scope(BRANCH, [], {"study_id": SID, "accessions": {"bioproject": "PRJNA1"}}, {})[0]


def test_shared_accession_warns():
    _, warnings = check_scope(BRANCH, [], STUDY, {"OTHER2020_X": {ACC}})
    assert warnings


def test_not_a_study_branch():
    assert check_scope("main", [], STUDY, {})[0]
