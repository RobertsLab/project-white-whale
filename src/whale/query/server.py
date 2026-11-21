"""`whale mcp`: White Whale as a read-only Model Context Protocol server (stdio).

The tools are thin wrappers over whale.query.tools, so every number an MCP client
(Claude Code, Claude Desktop, ...) reports comes from the committed results, with
the study IDs, run names and SNP mask it rests on. Nothing here writes.

Needs the optional extra: pip install -e ".[mcp]". Pattern from AREE's
src/query/mcp_server.py (see VENDORED_FROM_AREE.md). Setup: docs/mcp.md.
"""
from __future__ import annotations

import contextlib
import threading
from pathlib import Path

from whale.query import tools
from whale.query.store import Store

INSTRUCTIONS = """\
You are answering questions from White Whale: public Pacific oyster (Magallana gigas)
DNA-methylation studies re-processed from raw reads with one pipeline on one reference
genome (GCF_963853765.1), with C/T SNPs masked, and pooled into a gene-level atlas.

How to answer:
- Take numbers only from these tools; never from memory. Cite the study_id behind every
  number, and the nextflow_run and snp_mask (in each answer's `sources`) when precision matters.
- Start with whale_overview to see which studies and outputs exist. Outputs not yet produced
  for a study (differential, sample curation, replication) are listed as missing: say so
  rather than guessing.
- Gene symbols are mostly NCBI LOC ids. To find a gene by name, use find_genes with words
  from its product description (e.g. "DNA (cytosine-5)-methyltransferase", "heat shock 70").
  Many oyster gene families are large; say how many genes matched.
- A gene's study value is the mean over that study's units (libraries of one BioSample
  summed). A tissue, life stage or chemistry that comes from one study cannot be separated
  from that study: report design caveats, and the number of studies a statement rests on.
- WGBS and EM-seq are different chemistries; say when a comparison spans both.
- For differential results, report lambda: when `inflated` is true (lambda > 1.3), the
  significant-gene count is not a finding. Group sizes are small; never imply causation.
- Pages are capped; when `truncated` is true, say the list is partial or page with offset.
"""


def _server_class():
    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as exc:
        raise SystemExit(
            "`whale mcp` needs the 'mcp' package (Python >= 3.10). Install it with:\n"
            '  pip install -e ".[mcp]"'
        ) from exc
    return FastMCP


def _guard(fn, *args, **kwargs) -> dict:
    """A question the data cannot answer comes back as {'error': ...}, not a crash."""
    try:
        return fn(*args, **kwargs)
    except tools.QueryError as exc:
        return {"error": str(exc)}


def build_server(store: Store):
    """The MCP server with every White Whale tool registered."""
    server = _server_class()(name="whale", instructions=INSTRUCTIONS)

    @server.tool()
    def whale_overview() -> dict:
        """Start here. Studies registered and in the atlas (and any excluded, with the reason),
        units, libraries and genes pooled, the genome assembly, design caveats, and which
        outputs (differential, sample_curation, replication) exist for which studies."""
        return _guard(tools.overview, store)

    @server.tool()
    def list_studies(tissue: str | None = None, life_stage: str | None = None,
                     chemistry: str | None = None, stage: str | None = None) -> dict:
        """Studies with design summary: BioProject, chemistry (wgbs, emseq, ...), tissues, life
        stages, sample unit, factors, sample counts, QC pass rate, quality flags, pipeline stage,
        Nextflow run and SNP mask. Filters are registry terms, e.g. tissue='gill', life_stage='adult',
        chemistry='wgbs', stage='summarized'."""
        return _guard(tools.list_studies, store, tissue, life_stage, chemistry, stage)

    @server.tool()
    def describe_study(study_id: str) -> dict:
        """Everything about one study: citation, accessions, population, factors and levels,
        limitations, each recorded run's QC (mapping, conversion, CpGs at 5x), CpG methylation
        overall / by factor / by condition, gene-level summary, atlas inclusion, and whether
        differential, curation and replication outputs exist."""
        return _guard(tools.describe_study, store, study_id)

    @server.tool()
    def get_samples(study_id: str, limit: int = 100, offset: int = 0, full: bool = False) -> dict:
        """Per-library metadata for one study: accessions, condition, factor levels, deposit
        attributes, curated attributes (animal, tissue, life stage, sex, ploidy, pooling,
        population, treatment, confidence) when sample_curation.tsv exists, QC metrics and pass
        flag, and the atlas unit. full=true adds the curators' sources and notes."""
        return _guard(tools.get_samples, store, study_id, limit, offset, full)

    @server.tool()
    def find_genes(query: str, limit: int = 25, offset: int = 0) -> dict:
        """Find genes by NCBI GeneID (e.g. 105331391), symbol (e.g. LOC105331391) or words from
        the product description (all words must match, e.g. 'heat shock 70'). Returns GeneIDs,
        symbols, descriptions and locations."""
        return _guard(tools.find_genes, store, query, limit, offset)

    @server.tool()
    def get_gene(gene: str, include_groups: bool = True, include_differential: bool = True,
                 max_genes: int = 3) -> dict:
        """Methylation of a gene in every atlas study: gene-body and promoter mean %, units
        measured, range, gene-body means per factor level / condition / tissue within each
        study, and (when present) the shared-model differential results. `gene` is a GeneID,
        symbol or description words; at most max_genes (cap 10) matching genes are shown."""
        return _guard(tools.get_gene, store, gene, include_groups, include_differential, max_genes)

    @server.tool()
    def rank_genes(region: str = "gene_body", order: str = "high", by: str = "mean",
                   tissue: str | None = None, life_stage: str | None = None,
                   chemistry: str | None = None, study_ids: list[str] | None = None,
                   min_studies: int | None = None, min_study_pct: float | None = None,
                   max_study_pct: float | None = None, limit: int = 50, offset: int = 0) -> dict:
        """Genes ranked by methylation across studies, e.g. "hypermethylated in gill in every
        study": tissue='gill', min_study_pct=50. region: gene_body or promoter. order: high or
        low. by: mean, min or max of the study means. min_studies defaults to every selected
        study. min_study_pct / max_study_pct apply to each study's mean. Filters select units
        (tissue, life_stage, chemistry, study_ids)."""
        return _guard(tools.rank_genes, store, region, order, by, tissue, life_stage, chemistry,
                      study_ids, min_studies, min_study_pct, max_study_pct, limit, offset)

    @server.tool()
    def variable_genes(region: str = "gene_body", tissue: str | None = None,
                       life_stage: str | None = None, chemistry: str | None = None,
                       study_ids: list[str] | None = None, min_studies: int | None = None,
                       min_mean_pct: float | None = None, limit: int = 50, offset: int = 0) -> dict:
        """Genes whose methylation differs most between studies (SD of study means), with the
        range, each study's mean and the median within-study SD for comparison. Same filters
        as rank_genes; needs at least two studies."""
        return _guard(tools.variable_genes, store, region, tissue, life_stage, chemistry, study_ids,
                      min_studies, min_mean_pct, limit, offset)

    @server.tool()
    def list_contrasts(study_id: str | None = None) -> dict:
        """The shared DSS model's tests per study: job, formula, coefficient, region, genes
        tested, significant genes (FDR < 0.05, up/down, |diff| >= 10 pp), lambda (inflation) and
        the PCA-outlier sensitivity run. Lists studies without differential results."""
        return _guard(tools.list_contrasts, store, study_id)

    @server.tool()
    def get_differential(study_id: str, job_id: str | None = None, region: str = "gene_body",
                         coef: str | None = None, max_fdr: float | None = 0.05,
                         direction: str | None = None, limit: int = 50, offset: int = 0) -> dict:
        """Genes for one study's contrast, most significant first: group means, difference in
        percentage points, Wald statistic, p, FDR, plus lambda. job_id/coef from list_contrasts
        (default: job 'main', its first coefficient). max_fdr=null returns all tested genes.
        direction: 'up' or 'down'."""
        return _guard(tools.get_differential, store, study_id, job_id, region, coef, max_fdr,
                      direction, limit, offset)

    @server.tool()
    def gene_differential(gene: str) -> dict:
        """One gene (GeneID, symbol, description words, or comma-separated GeneIDs) in every
        study's differential tables: per job, region and coefficient, diff, p, FDR and lambda."""
        return _guard(tools.gene_differential, store, gene)

    @server.tool()
    def get_replication(study_id: str | None = None, verdict: str | None = None,
                        kind: str | None = None, text: str | None = None,
                        limit: int = 50, offset: int = 0) -> dict:
        """Did published findings replicate on the re-derived data? Per study: verdict counts,
        caveats, papers. Per claim: statement, where in the paper, verdict (replicated, partly,
        not_replicated, inconclusive, not_testable, pending), the re-derived numbers and the
        rule applied. Filter by study, verdict, claim kind or text in the statement."""
        return _guard(tools.get_replication, store, study_id, verdict, kind, text, limit, offset)

    return server


TOOL_NAMES = {"whale_overview", "list_studies", "describe_study", "get_samples", "find_genes", "get_gene",
              "rank_genes", "variable_genes", "list_contrasts", "get_differential", "gene_differential",
              "get_replication"}


def run(repo: Path, transport: str = "stdio") -> None:
    store = Store(repo)
    server = build_server(store)
    # Pooling the gene tables takes several seconds; start it now so the first question
    # does not wait for all of it. Errors surface again on the first call that needs it.
    threading.Thread(target=_warm, args=(store,), daemon=True).start()
    server.run(transport=transport)


def _warm(store: Store) -> None:
    # A failure here is reported by the tool call that needs the atlas.
    with contextlib.suppress(Exception):
        store.atlas()
