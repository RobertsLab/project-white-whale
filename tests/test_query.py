"""whale.query: the read-only tools behind `whale mcp`, on a small synthetic repository."""
from __future__ import annotations

import asyncio
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from whale.query import tools
from whale.query.store import Store

ASSEMBLY = "GCF_TEST.1"
MASK = "GCF_TEST.1.snpmask.v1"
GENES = {
    "1001": ("dnmt1", "DNA (cytosine-5)-methyltransferase 1"),
    "1002": ("LOC1002", "heat shock 70 kDa protein"),
    "1003": ("LOC1003", "heat shock 70 kDa protein 12A"),
    "1004": ("LOC1004", "uncharacterized LOC1004"),
}
# % methylation per gene and library; cpgs per gene (under 5 = not measured)
ALPHA = {"libs": ["A1", "A2", "A3", "A4"], "heat": ["C", "C", "H", "H"],
         "pct": {"1001": [80, 82, 78, 80], "1002": [2, 4, 2, 4], "1003": [50, 50, 70, 70], "1004": [90, 90, 92, 88]},
         "cpgs": {"1003": 10}}
BETA = {"libs": ["B1", "B2"], "heat": None,
        "pct": {"1001": [76, 78], "1002": [6, 6], "1003": [40, 40], "1004": [10, 12]},
        "cpgs": {"1003": 2}}


def _gz(path: Path, rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def _tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))


def _study(repo: Path, sid: str, spec: dict, *, chemistry: str, tissue: str, mask: str | None, run: str):
    libs = spec["libs"]
    factors = [{"factor_id": "heat", "factor_original": "heat", "factor_standardized": "temperature",
                "levels": ["C", "H"], "reference_level": "C"}] if spec["heat"] else []
    (repo / "registry" / "studies").mkdir(parents=True, exist_ok=True)
    (repo / "registry" / "studies" / f"{sid}.yaml").write_text(yaml.safe_dump({
        "study_id": sid, "accessions": {"bioproject": "PRJ" + sid[:4]}, "citation": f"{sid} et al.",
        "chemistry": chemistry, "samples": {"n_runs": len(libs), "n_libraries": len(libs),
                                            "sample_unit": "individual", "tissues": [tissue],
                                            "life_stages": ["adult"]},
        "factors": factors, "paired_rnaseq": {"relationship": "none"}, "quality_flags": ["sex_unknown"],
        "limitations": "Small.", "goals": [1]}))
    d = repo / "data" / "studies" / sid
    for region in ("gene_body", "promoter"):
        meth, calls, cpgs = [["gene_id", *libs]], [["gene_id", *libs]], [["gene_id", *libs]]
        for g in GENES:
            pct = spec["pct"][g] if region == "gene_body" else [p / 2 for p in spec["pct"][g]]
            calls.append([g] + [1000] * len(libs))
            meth.append([g] + [int(10 * p) for p in pct])
            cpgs.append([g] + [spec["cpgs"].get(g, 20)] * len(libs))
        for metric, rows in (("meth", meth), ("calls", calls), ("cpgs", cpgs)):
            _gz(d / "genes" / f"{region}.{metric}.tsv.gz", rows)
    _tsv(d / "samplesheet.tsv", [{"sample_id": lib, "sample_accession": f"SAM{lib}",
                                  "condition": spec["heat"][i] if spec["heat"] else "all",
                                  "heat": spec["heat"][i] if spec["heat"] else "", "tissue": tissue}
                                 for i, lib in enumerate(libs)])
    _tsv(d / "qc_summary.tsv", [{"sample_id": lib, "condition": "x", "mapping_pct": 60, "conversion_pct": 99.5,
                                 "cpg_meth_pct": 12.5, "cpgs_5x": 5000000, "mean_depth": 10, "qc_pass": "true",
                                 "qc_reasons": ""} for lib in libs])
    summary = {"study_id": sid, "part": None, "generated_on": "2026-11-01", "nextflow_run": run,
               "genome_assembly": ASSEMBLY, "chemistry": chemistry, "libraries": len(libs) + 1,
               "passed": len(libs), "failed": [{"sample_id": "X", "reasons": "conversion"}],
               "by_factor": [{"factor_id": "heat", "column": "heat", "levels": []}] if spec["heat"] else []}
    if mask:
        summary["snp_mask"] = {"mask_id": mask, "cpgs_masked": 10}
    (d / "study_summary.yaml").write_text(yaml.safe_dump(summary))
    (d / "status.yaml").write_text(yaml.safe_dump({"study_id": sid, "stage": "summarized"}))
    return d


@pytest.fixture(scope="module")
def store(tmp_path_factory) -> Store:
    repo = tmp_path_factory.mktemp("whale_repo")
    (repo / "registry" / "controlled_vocabularies").mkdir(parents=True)
    _gz(repo / "data" / "reference" / "gene_info" / f"{ASSEMBLY}.tsv.gz",
        [["gene_id", "symbol", "description", "biotype", "chrom", "start", "end", "strand"]]
        + [[g, s, desc, "protein_coding", "chr1", 100 * i, 100 * i + 50, "+"]
           for i, (g, (s, desc)) in enumerate(GENES.items())])
    a = _study(repo, "ALPHA2020_GILL_WGBS", ALPHA, chemistry="wgbs", tissue="gill", mask=MASK, run="brave_turing")
    _study(repo, "BETA2021_MANTLE_EMSEQ", BETA, chemistry="emseq", tissue="mantle", mask=MASK, run="calm_hopper")
    _study(repo, "GAMMA2022_NOMASK_WGBS", BETA | {"libs": ["G1", "G2"]}, chemistry="wgbs", tissue="gill",
           mask=None, run="old_run")

    # Optional outputs, for ALPHA only (formats of PRs #60, #62, #63).
    _tsv(a / "sample_curation.tsv", [{"sample_id": lib, "animal_id": f"an{lib}", "tissue": "gill",
                                      "life_stage": "adult", "age": "", "sex": "unknown", "ploidy": "diploid",
                                      "sample_unit": "individual", "pool_size": "", "population": "P",
                                      "treatment": "t", "confidence": "high", "sources": "s", "notes": "n"}
                                     for lib in ALPHA["libs"]])
    dd = a / "differential"
    _tsv(dd / "tests.tsv", [
        {"job_id": "main", "analysis": "all_units", "units": 4, "outliers": "", "region": "gene_body",
         "genes_tested": 4, "coef": "heatH", "n_sig": 1, "up": 1, "down": 0, "large": 1, "lambda": 1.5, "note": ""},
        {"job_id": "main", "analysis": "without_pca_outliers", "units": 3, "outliers": "SAMA4", "region": "gene_body",
         "genes_tested": 4, "coef": "heatH", "n_sig": 0, "up": 0, "down": 0, "large": 0, "lambda": 1.1, "note": ""}])
    (dd / "summary.yaml").write_text(yaml.safe_dump({"study_id": "ALPHA2020_GILL_WGBS", "generated_on": "2026-11-20",
                                                     "snp_mask": MASK, "gene_tables_from_run": "brave_turing",
                                                     "jobs": [{"job_id": "main", "formula": "~ heat"}]}))
    _gz(dd / "main.gene_body.tsv.gz", [
        ["gene_id", "symbol", "description", "mean_pct", "mean_heat=C", "mean_heat=H",
         "diff_heatH", "stat_heatH", "p_heatH", "fdr_heatH"],
        ["1003", "LOC1003", GENES["1003"][1], "60", "50", "70", "20", "5.1", "1e-6", "0.004"],
        ["1001", "dnmt1", GENES["1001"][1], "80", "81", "79", "-2", "-1.2", "0.2", "0.4"],
        ["1002", "LOC1002", GENES["1002"][1], "3", "3", "3", "0", "0.1", "0.9", "0.9"]])
    (a / "replication.json").write_text(json.dumps({
        "study_id": "ALPHA2020_GILL_WGBS", "generated_on": "2026-11-21",
        "sources": [{"id": "P1", "citation": "Alpha 2020", "doi": "10.1/x"}],
        "caveats": ["test caveat"], "based_on": {"nextflow_run": "brave_turing", "snp_mask": MASK},
        "counts": {"replicated": 1, "not_testable": 1},
        "claims": [{"id": "C1", "source": "P1", "kind": "global_level", "statement": "Global CpG methylation 12%.",
                    "verdict": "replicated", "re_derived": "12.5%", "rule": "within 2 points"},
                   {"id": "C2", "source": "P1", "kind": "named_gene", "statement": "HSP70 is hypermethylated.",
                    "verdict": "not_testable", "reason": "no such contrast"}]}))
    return Store(repo)


A, B = "ALPHA2020_GILL_WGBS", "BETA2021_MANTLE_EMSEQ"


def test_overview_reports_atlas_exclusions_and_outputs(store):
    o = tools.overview(store)
    assert o["atlas"]["studies"] == [A, B]
    assert o["atlas"]["units"] == 6 and o["atlas"]["genes"] == 4
    assert o["atlas"]["excluded"] == [{"study_id": "GAMMA2022_NOMASK_WGBS", "part": None,
                                       "reason": "gene tables not SNP-masked"}]
    assert o["outputs_present"] == {"differential": [A], "sample_curation": [A], "replication": [A]}
    assert {s["study_id"]: s["nextflow_run"] for s in o["sources"]} == {A: "brave_turing", B: "calm_hopper"}
    assert all(s["snp_mask"] == MASK for s in o["sources"])


def test_list_and_describe_studies(store):
    out = tools.list_studies(store, tissue="gill")
    assert [s["study_id"] for s in out["studies"]] == [A, "GAMMA2022_NOMASK_WGBS"]
    alpha = out["studies"][0]
    assert alpha["qc"] == {"libraries": 5, "passed": 4, "pass_rate": 0.8}
    assert alpha["quality_flags"] == ["sex_unknown"] and alpha["stage"] == "summarized"
    assert tools.list_studies(store, chemistry="emseq")["count"] == 1

    d = tools.describe_study(store, "GAMMA2022_NOMASK_WGBS")
    assert d["atlas"] == {"included": False, "units": 0, "excluded_reason": "gene tables not SNP-masked"}
    d = tools.describe_study(store, A)
    assert d["available"]["differential"] and d["available"]["replication"]["counts"]["replicated"] == 1
    assert d["runs"][0]["nextflow_run"] == "brave_turing"
    with pytest.raises(tools.QueryError, match="Unknown study_id"):
        tools.describe_study(store, "NOPE")


def test_samples_curated_and_paginated(store):
    s = tools.get_samples(store, A, limit=3)
    assert s["page"] == {"total": 4, "offset": 0, "returned": 3, "limit": 3, "cap": tools.MAX_LIMIT,
                         "truncated": True, "next_offset": 3}
    first = s["samples"][0]
    assert first["factors"] == {"heat": "C"} and first["curated"]["animal_id"] == "anA1"
    assert "sources" not in first["curated"] and first["qc"]["qc_pass"] == "true"
    assert first["atlas_unit"] == "SAMA1"
    assert "sources" in tools.get_samples(store, A, full=True)["samples"][0]["curated"]
    assert tools.get_samples(store, B)["samples"][0].get("curated") is None
    assert tools.get_samples(store, A, limit=10_000)["page"]["limit"] == tools.MAX_LIMIT


def test_gene_lookup_by_id_symbol_and_description(store):
    assert [g["gene_id"] for g in tools.find_genes(store, "1004")["genes"]] == ["1004"]
    assert [g["gene_id"] for g in tools.find_genes(store, "DNMT1")["genes"]] == ["1001"]
    assert [g["gene_id"] for g in tools.find_genes(store, "heat shock 70")["genes"]] == ["1002", "1003"]

    g = tools.get_gene(store, "LOC1003")
    gene = g["genes"][0]
    alpha, beta = gene["studies"]
    assert alpha["gene_body"]["mean_pct"] == 60.0 and alpha["gene_body"]["units_measured"] == 4
    assert alpha["promoter"]["mean_pct"] == 30.0
    assert beta["gene_body"]["mean_pct"] is None and beta["gene_body"]["units_measured"] == 0
    groups = {(x["group_by"], x["level"]): x["mean_pct"] for x in alpha["gene_body_groups"]}
    assert groups[("heat", "C")] == 50.0 and groups[("heat", "H")] == 70.0
    assert g["differential"][0]["study_id"] == A and g["differential"][0]["fdr"] == 0.004
    assert {s["study_id"] for s in g["sources"]} == {A, B}

    many = tools.get_gene(store, "heat shock", max_genes=1)
    assert many["matched"] == 2 and many["shown"] == 1 and many["truncated"]


def test_rank_genes_with_filters(store):
    r = tools.rank_genes(store, min_study_pct=50)
    # 1003 is not measured in BETA, so it is out when every study is required
    assert [g["gene_id"] for g in r["genes"]] == ["1001"]
    assert r["genes"][0]["per_study"] == {A: 80.0, B: 77.0}

    gill = tools.rank_genes(store, tissue="gill", limit=2)
    assert [s["study_id"] for s in gill["studies"]] == [A]
    assert [g["gene_id"] for g in gill["genes"]] == ["1004", "1001"]
    assert gill["page"]["truncated"] and gill["page"]["next_offset"] == 2

    low = tools.rank_genes(store, order="low", region="promoter", min_studies=1)
    assert low["genes"][0]["gene_id"] == "1002" and low["genes_passing"] == 4
    assert tools.rank_genes(store, chemistry="emseq")["studies"] == [{"study_id": B, "units": 2}]
    with pytest.raises(tools.QueryError, match="no atlas units"):
        tools.rank_genes(store, tissue="heart")
    with pytest.raises(tools.QueryError, match="region"):
        tools.rank_genes(store, region="exon")


def test_variable_genes(store):
    v = tools.variable_genes(store)
    assert v["genes"][0]["gene_id"] == "1004"
    assert v["genes"][0]["range"] == 79.0
    assert any("chemistries" in c for c in v["caveats"])
    with pytest.raises(tools.QueryError, match="two studies"):
        tools.variable_genes(store, study_ids=[A])


def test_differential_for_a_contrast_and_a_gene(store):
    c = tools.list_contrasts(store)
    assert c["without_differential"] == [B, "GAMMA2022_NOMASK_WGBS"]
    t = c["studies"][0]["tests"][0]
    assert t["lambda"] == 1.5 and t["inflated"] and t["n_sig_fdr_0.05"] == 1

    d = tools.get_differential(store, A)
    assert d["job_id"] == "main" and d["coef"] == "heatH" and d["inflated"]
    assert [g["gene_id"] for g in d["genes"]] == ["1003"]
    assert d["genes"][0]["mean_heat=H"] == 70 and d["genes"][0]["diff_pp"] == 20
    assert d["sensitivity_without_pca_outliers"]["outliers"] == "SAMA4"
    every = tools.get_differential(store, A, max_fdr=None, direction="down")
    assert [g["gene_id"] for g in every["genes"]] == ["1001"]
    assert tools.get_differential(store, B)["available"] is False

    g = tools.gene_differential(store, "dnmt1")
    assert [(r["study_id"], r["coef"], r["p"]) for r in g["results"]] == [(A, "heatH", 0.2)]
    assert g["studies_without_differential"] == [B, "GAMMA2022_NOMASK_WGBS"]


def test_replication(store):
    r = tools.get_replication(store)
    assert [s["study_id"] for s in r["studies"]] == [A]
    assert len(r["claims"]) == 2 and r["claims"][0]["study_id"] == A
    assert B in r["studies_without_replication"]
    assert [c["id"] for c in tools.get_replication(store, verdict="not_testable")["claims"]] == ["C2"]
    assert [c["id"] for c in tools.get_replication(store, text="hsp70")["claims"]] == ["C2"]
    with pytest.raises(tools.QueryError):
        tools.get_replication(store, verdict="maybe")


def test_tools_never_write(store):
    before = {p: p.stat().st_mtime_ns for p in store.repo.rglob("*")}
    tools.overview(store)
    tools.get_gene(store, "1001")
    tools.get_differential(store, A)
    assert {p: p.stat().st_mtime_ns for p in store.repo.rglob("*")} == before


def test_answers_are_json(store):
    for out in (tools.overview(store), tools.get_gene(store, "1001"), tools.variable_genes(store),
                tools.describe_study(store, A)):
        json.dumps(out, allow_nan=False)


def test_mcp_server_lists_tools_and_answers(store):
    pytest.importorskip("mcp")
    from whale.query.server import TOOL_NAMES, build_server

    server = build_server(store)

    async def exercise():
        names = {t.name for t in await server.list_tools()}
        result = await server.call_tool("get_gene", {"gene": "1001", "include_differential": False})
        bad = await server.call_tool("describe_study", {"study_id": "NOPE"})
        return names, result, bad

    names, result, bad = asyncio.run(exercise())
    assert names == TOOL_NAMES
    # Depending on the mcp 1.x release, call_tool returns the content blocks or (content, structured).
    content = result[0] if isinstance(result, tuple) else result
    assert json.loads(content[0].text)["genes"][0]["studies"][0]["gene_body"]["mean_pct"] == 80.0
    bad = bad[0] if isinstance(bad, tuple) else bad
    assert "Unknown study_id" in json.loads(bad[0].text)["error"]


def test_whale_mcp_starts_over_stdio(store):
    """The real command: start `whale mcp` as a subprocess and list its tools over stdio."""
    pytest.importorskip("mcp")
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    from whale.query.server import TOOL_NAMES

    params = StdioServerParameters(command=sys.executable,
                                   args=["-c", "from whale.cli import main; main()", "mcp",
                                         "--repo", str(store.repo)])

    async def exercise():
        async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
            await session.initialize()
            listed = await session.list_tools()
            res = await session.call_tool("list_studies", {"chemistry": "emseq"})
            return {t.name for t in listed.tools}, json.loads(res.content[0].text)

    names, studies = asyncio.run(asyncio.wait_for(exercise(), timeout=60))
    assert names == TOOL_NAMES
    assert [s["study_id"] for s in studies["studies"]] == [B]


def test_cli_without_mcp_says_how_to_install(monkeypatch):
    import builtins

    real = builtins.__import__

    def fake(name, *a, **k):
        if name.startswith("mcp"):
            raise ImportError(name)
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake)
    from whale.query import server

    with pytest.raises(SystemExit, match=r"\.\[mcp\]"):
        server._server_class()


def test_cli_registers_mcp_command():
    out = subprocess.run([sys.executable, "-c", "from whale.cli import main; main()", "mcp", "--help"],
                         capture_output=True, text=True, check=True).stdout
    assert "stdio" in out and "--repo" in out
