# Asking Claude about White Whale: `whale mcp`

`whale mcp` is a read-only [Model Context Protocol](https://modelcontextprotocol.io)
server over the summarized results (PLAN.md §3, §8). Connected to Claude Code or Claude
Desktop, it lets Claude answer questions such as "which genes are hypermethylated in
gill across every study?" from the committed tables, with the study IDs, Nextflow runs
and SNP mask each number rests on.

It reads files and nothing else. It does not write, call git or use the network.

## Install

The MCP SDK is an optional dependency, so the core install stays light. You need
Python 3.10 or newer.

```bash
pip install -e ".[mcp]"
```

Check that it starts. It waits for a client on stdin; press Ctrl-C to stop it.

```bash
whale mcp --repo /path/to/project-white-whale
```

Without `--repo`, it serves the checkout around the current directory, or `$WHALE_REPO`.
MCP clients start servers from their own working directory, so pass `--repo` in the
client configuration.

## Add it to Claude Code

Use the absolute path of the `whale` executable in the environment where you installed
the extra (`which whale`):

```bash
claude mcp add whale --scope user -- /path/to/.venv/bin/whale mcp --repo /path/to/project-white-whale
```

`--scope user` makes it available in every project; leave it out to add it to the
current project only. `claude mcp list` should then show `whale` as connected, and
`/mcp` in a session lists its tools.

## Add it to Claude Desktop

Open **Settings → Developer → Edit Config**, which opens
`claude_desktop_config.json`. On macOS the file is in
`~/Library/Application Support/Claude/`. Add the server, then restart Claude Desktop:

```json
{
  "mcpServers": {
    "whale": {
      "command": "/path/to/.venv/bin/whale",
      "args": ["mcp", "--repo", "/path/to/project-white-whale"]
    }
  }
}
```

## What it serves

| Tool | Answers |
|---|---|
| `whale_overview` | Studies registered and in the atlas (and any excluded, with the reason); units, libraries and genes; design caveats; which optional outputs exist |
| `list_studies` | Design summary per study: chemistry, tissues, life stages, sample unit, factors, counts, QC pass rate, quality flags, stage, run, mask. Filter by tissue, life stage, chemistry or stage |
| `describe_study` | One study in full: citation, accessions, factors, limitations, each run's QC and CpG methylation by factor, atlas inclusion, available outputs |
| `get_samples` | Per-library metadata: accessions, factor levels, deposit attributes, curated attributes (`sample_curation.tsv`), QC, atlas unit |
| `find_genes` | Genes by NCBI GeneID, symbol or product-description words |
| `get_gene` | A gene in every atlas study: gene-body and promoter mean, units measured, per-group means within each study, and differential results |
| `rank_genes` | Genes ranked by methylation across studies, with filters (tissue, life stage, chemistry, studies, region) and per-study thresholds |
| `variable_genes` | Genes whose methylation differs most between studies (SD of study means), beside the within-study SD |
| `list_contrasts` | The shared DSS model's tests per study: counts, λ, PCA-outlier sensitivity |
| `get_differential` | Genes for one study's contrast, most significant first |
| `gene_differential` | One gene in every study's differential tables |
| `get_replication` | Replication verdicts per study and per published claim, with the rule and numbers behind each |

Each answer has a `sources` list with the `study_id`, `nextflow_run`, `snp_mask`,
`genome_assembly` and summary file of every study it draws on. Lists are paginated:
`limit` (at most 200) and `offset`, with a `page` block reporting the total, the cap and
whether the list was truncated. `get_gene` shows at most 10 genes per call.

### Where the numbers come from

| Source | Used for |
|---|---|
| `registry/studies/*.yaml` | Design, chemistry, tissues, life stages, quality flags, limitations |
| `data/studies/<ID>/study_summary.yaml`, `status.yaml` | Recorded run, QC, SNP mask, stage |
| `data/studies/<ID>/samplesheet.tsv`, `qc_summary.tsv` | Per-library metadata and QC |
| `data/studies/<ID>/genes/*.tsv.gz` | Gene methylation, pooled in memory the same way as `whale build-atlas` (`whale.atlas.pool` and `units`) |
| `data/studies/<ID>/sample_curation.tsv` | Curated sample attributes, when present; they take precedence over the deposit's for filtering |
| `data/studies/<ID>/differential/` | `whale differential` results, when present |
| `data/studies/<ID>/replication.json` | `whale replicate` verdicts, when present |
| `data/reference/gene_info/<ASSEMBLY>.tsv.gz` | Gene symbols, descriptions, locations |
| `site/atlas/data/atlas.json` | The atlas headline, if you have run `whale build-atlas` locally |

A study without one of the optional outputs is listed as missing it, for example under
`without_differential` or `studies_without_replication`. The server never reports an
estimate in its place.

Gene values follow the atlas's rules (docs/atlas.md):

- Libraries of one BioSample are summed into one unit.
- A unit's value for a gene needs at least 5 covered CpGs.
- A study's value is the mean over its units. It is reported only when at least half of
  the units measure the gene.

The server notices when the committed files change, for example after a `git pull`, and
re-pools on the next question. Pooling takes several seconds. It starts when the server
starts, so the first question usually does not wait for all of it.

## Example

The question "which genes are hypermethylated in gill across every study?" becomes
`rank_genes(tissue="gill", min_study_pct=50, limit=5)`. This is the answer from `main`
at the time of writing, shortened to the third of the five genes and the first source:

```json
{
  "region": "gene_body", "ranked_by": "mean of study means",
  "filters": {"tissue": "gill", "min_study_pct": 50.0},
  "studies": [{"study_id": "ROBERTS2022_PH_PLOIDY_WGBS", "units": 24},
              {"study_id": "WANG2020_TIDAL_HEAT_WGBS", "units": 12},
              {"study_id": "WANG2023_TIDAL_TRANSGEN_WGBS", "units": 24}],
  "min_studies": 3, "genes_passing": 4221,
  "genes": [{"gene_id": "136269422", "symbol": "LOC136269422",
             "description": "transcriptional adapter 2-alpha-like", "mean_pct": 88.99,
             "per_study": {"ROBERTS2022_PH_PLOIDY_WGBS": 92.17, "WANG2020_TIDAL_HEAT_WGBS": 86.36,
                           "WANG2023_TIDAL_TRANSGEN_WGBS": 88.43}}],
  "page": {"total": 4221, "returned": 5, "cap": 200, "truncated": true, "next_offset": 5},
  "notes": ["VALDIVIESO2025_POMS_RS_EMSEQ: 40 unit(s) left out because the filtered attribute is not recorded"],
  "sources": [{"study_id": "ROBERTS2022_PH_PLOIDY_WGBS", "nextflow_run": "curious_mestorf",
               "snp_mask": "GCF_963853765.1.snpmask.v1", "genome_assembly": "GCF_963853765.1"}]
}
```

The note matters. VALDIVIESO2025 sampled gill too, but its tissue is not recorded per
library on `main`, so it is left out and the note says so. Once its
`sample_curation.tsv` is merged, its gill libraries are included.

## Development

The tools are plain functions in `src/whale/query/tools.py` that take a
`whale.query.Store`. `src/whale/query/server.py` registers them with FastMCP, following
AREE's `src/query/mcp_server.py`. To add a question, write the function, register it in
`build_server`, add its name to `TOOL_NAMES`, and test it in `tests/test_query.py`. Those
tests build a small synthetic repository; one starts `whale mcp` over stdio and lists its
tools.
