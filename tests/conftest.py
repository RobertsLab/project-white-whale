import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]

# What a claim needs from the repository; copied into throwaway git repos.
SEED_PATHS = ["registry", "schemas", "config", ".github/PULL_REQUEST_TEMPLATE"]


@pytest.fixture
def repo_root() -> Path:
    return REPO


@pytest.fixture
def template_study() -> dict:
    return yaml.safe_load((REPO / "registry/studies/_TEMPLATE.yaml").read_text())


@pytest.fixture
def filled_study(template_study) -> dict:
    s = template_study
    s["study_id"] = "WANG2020_TIDAL_HEAT_WGBS"
    s["accessions"]["bioproject"] = "PRJNA562805"
    s["citation"] = "Wang X et al. (2020) Heredity 126:10-22"
    s["chemistry"] = "wgbs"
    s["library_layout"] = "paired"
    s["samples"].update(n_runs=12, n_libraries=12, sample_unit="pool", pool_size=5,
                        tissues=["gill"], life_stages=["adult"])
    s["paired_rnaseq"].update(relationship="same_animals", bioproject="PRJNA407831")
    s["factors"] = [
        {"factor_id": "origin", "factor_original": "intertidal vs subtidal F1",
         "factor_standardized": "genotype_or_population", "levels": ["Sub", "In"], "reference_level": "Sub"},
        {"factor_id": "heat", "factor_original": "35 C for 6 h vs 0 h",
         "factor_standardized": "temperature", "levels": ["C", "H"], "reference_level": "C"},
    ]
    s["goals"] = [1, 2]
    s["data_availability"].update(raw_data_available=True, fastq_size_gb=140.7)
    s["quality_flags"] = ["pooled_libraries", "shared_biosample", "confounded_design"]
    s["limitations"] = "One population per tidal class."
    s["provenance"].update(registered_by="tester", date_registered="2026-10-23")
    return s


def git(*args, cwd):
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True).stdout


@pytest.fixture
def claim_env(tmp_path, monkeypatch):
    """A bare 'origin' and a clone of it seeded with the registry, as a second machine would have."""
    for k, v in {"GIT_AUTHOR_NAME": "Tester", "GIT_AUTHOR_EMAIL": "t@example.org",
                 "GIT_COMMITTER_NAME": "Tester", "GIT_COMMITTER_EMAIL": "t@example.org"}.items():
        monkeypatch.setenv(k, v)
    monkeypatch.delenv("WHALE_REPO", raising=False)
    remote = tmp_path / "origin.git"
    git("init", "--quiet", "--bare", "-b", "main", str(remote), cwd=tmp_path)
    clone = tmp_path / "machine-a"
    git("clone", "--quiet", str(remote), str(clone), cwd=tmp_path)
    git("checkout", "--quiet", "-b", "main", cwd=clone)
    for rel in SEED_PATHS:
        src, dest = REPO / rel, clone / rel
        if src.is_dir():
            shutil.copytree(src, dest)
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
    git("add", "-A", cwd=clone)
    git("commit", "--quiet", "-m", "seed", cwd=clone)
    git("push", "--quiet", "-u", "origin", "main", cwd=clone)
    return {"remote": remote, "clone": clone, "tmp": tmp_path}
