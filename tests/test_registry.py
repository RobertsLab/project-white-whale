import csv
import re

import yaml

from whale.claims import STUDY_ID_RE
from whale.registry import build_index, load_backlog
from whale.vocab import load_vocab


def test_backlog_rows_use_known_terms(repo_root):
    rows = load_backlog(repo_root)
    assert rows
    chem = load_vocab(repo_root, "assay_chemistry")
    profiles = load_vocab(repo_root, "compute_profile")
    accessions = [r["accession"] for r in rows]
    assert len(accessions) == len(set(accessions))
    for r in rows:
        assert re.match(r"^PRJ[NED][AB][0-9]+$", r["accession"]), r
        sid = r["suggested_study_id"]
        assert sid == "(name on claim)" or STUDY_ID_RE.match(sid), r
        assert r["chemistry"] in chem, r
        assert r["suggested_profile"] in profiles, r


def test_vocab_files_have_terms(repo_root):
    for path in (repo_root / "registry/controlled_vocabularies").glob("*.yaml"):
        doc = yaml.safe_load(path.read_text())
        assert doc["terms"], path


def test_stage_order_is_increasing(repo_root):
    stages = load_vocab(repo_root, "pipeline_stage")
    orders = [v["order"] for v in stages.values()]
    assert orders == sorted(orders)
    assert len(set(orders)) == len(orders)


def test_build_index(repo_root, tmp_path, filled_study, monkeypatch):
    out = tmp_path / "index.csv"
    build_index(repo_root, out)
    with out.open() as fh:
        header = next(csv.reader(fh))
    assert header[0] == "study_id"
