import copy

import yaml

from whale.validate import validate_status, validate_study, validate_study_file


def test_template_has_no_errors(repo_root):
    res = validate_study_file(repo_root / "registry/studies/_TEMPLATE.yaml", repo_root)
    assert res.errors == []
    assert any("chemistry is not filled in" in w for w in res.warnings)


def test_filled_study_is_clean(repo_root, filled_study):
    errors, warnings = validate_study(filled_study, repo_root)
    assert errors == []
    assert warnings == []


def test_every_registered_study_is_valid(repo_root):
    for path in sorted((repo_root / "registry/studies").glob("*.yaml")):
        assert validate_study_file(path, repo_root).errors == [], path


def _errors(repo_root, study, **changes):
    s = copy.deepcopy(study)
    for dotted, value in changes.items():
        node = s
        *parents, last = dotted.split("__")
        for p in parents:
            node = node[p]
        node[last] = value
    return validate_study(s, repo_root)[0]


def test_unknown_terms_are_errors(repo_root, filled_study):
    assert any("chemistry" in e for e in _errors(repo_root, filled_study, chemistry="nanopore"))
    assert any("tissues" in e for e in _errors(repo_root, filled_study, samples__tissues=["liver"]))
    assert any("quality_flags" in e for e in _errors(repo_root, filled_study, quality_flags=["pooled_libraries", "made_up"]))


def test_pooling_must_be_flagged(repo_root, filled_study):
    errs = _errors(repo_root, filled_study, quality_flags=["shared_biosample"])
    assert any("pooled_libraries" in e for e in errs)


def test_enrichment_and_rrbs_must_be_flagged(repo_root, filled_study):
    assert any("enrichment_assay" in e for e in _errors(repo_root, filled_study, chemistry="medip"))
    assert any("rrbs_subset" in e for e in _errors(repo_root, filled_study, chemistry="rrbs"))


def test_goal_2_needs_rnaseq(repo_root, filled_study):
    errs = _errors(repo_root, filled_study, paired_rnaseq__relationship="none")
    assert any("goals includes 2" in e for e in errs)


def test_reference_level_must_be_a_level(repo_root, filled_study):
    s = copy.deepcopy(filled_study)
    s["factors"][0]["reference_level"] = "Mid"
    assert any("reference_level" in e for e in validate_study(s, repo_root)[0])


def test_study_id_pattern_and_filename(repo_root, filled_study, tmp_path):
    assert _errors(repo_root, filled_study, study_id="wang2020")
    path = tmp_path / "OTHER2020_X.yaml"
    path.write_text(yaml.safe_dump(filled_study))
    res = validate_study_file(path, repo_root)
    assert any("does not match the file name" in e for e in res.errors)


def test_config_chemistry_must_match_the_yaml(filled_study, tmp_path):
    from whale.claims import config_chemistry, render_config
    from whale.validate import config_chemistry_errors

    sid = filled_study["study_id"]
    rendered = render_config(sid, "PRJNA562805", "raven")
    assert config_chemistry(rendered) is None
    assert config_chemistry("params {\n    chemistry = \"emseq\" // x\n}") == "emseq"

    (tmp_path / "config").mkdir()
    (tmp_path / "config" / f"{sid}.config").write_text(rendered)
    (tmp_path / "config" / f"{sid}.part-gill.config").write_text(rendered.replace("null", "'wgbs'", 1))
    errs = config_chemistry_errors(filled_study, tmp_path)
    assert len(errs) == 1 and f"config/{sid}.config sets chemistry = None" in errs[0]
    # A YAML still without chemistry has nothing to compare against yet.
    assert config_chemistry_errors({**filled_study, "chemistry": None}, tmp_path) == []


def test_malformed_yaml(repo_root, tmp_path):
    path = tmp_path / "BAD2020_X.yaml"
    path.write_text("study_id: [unclosed\n")
    assert validate_study_file(path, repo_root).errors


def test_status_validation(repo_root):
    status = {
        "study_id": "WANG2020_TIDAL_HEAT_WGBS", "accession": "PRJNA562805", "part": None,
        "stage": "aligned", "profile": "raven", "machine": "raven", "claimed_by": "t",
        "claimed_on": "2026-10-23", "nextflow_run": None, "notes": None,
        "history": [{"stage": "claimed", "date": "2026-10-23", "machine": "raven", "note": None},
                    {"stage": "aligned", "date": "2026-10-24", "machine": "raven", "note": None}],
    }
    assert validate_status(status, repo_root) == []
    assert validate_status({**status, "profile": "laptop"}, repo_root)
    assert validate_status({**status, "stage": "called"}, repo_root)


def test_cli_reports_unparseable_status_yaml(repo_root, tmp_path, monkeypatch):
    from click.testing import CliRunner

    from whale.cli import main

    monkeypatch.setenv("WHALE_REPO", str(repo_root))
    bad = tmp_path / "status.yaml"
    bad.write_text("stage: [claimed\n")
    res = CliRunner().invoke(main, ["validate-status", str(bad)])
    assert res.exit_code == 1
    assert "FAIL" in res.output and "could not parse YAML" in res.output
