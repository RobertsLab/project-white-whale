# Sample metadata curation — brief for one curator

You are one of two curators working on the same study **independently**. A
third person compares the two curations afterwards, so the value of your work
is that it is your own reading of the primary sources. Precision beats
completeness: an honest "unknown" with a note is a good answer; a plausible
guess is a bad one.

## The project

White Whale re-processes every public DNA-methylation dataset of the Pacific
oyster (*Crassostrea / Magallana gigas*) and pools them. Its analyses split
methylation variance by tissue, life stage, sex, ploidy, population and
treatment across studies, so each sequencing library needs those attributes
in a shared vocabulary. Today they are often missing or wrong.

## What you may read

- **The study's sample sheet**:
  `data/studies/<STUDY_ID>/samplesheet.tsv`.
  One row per library (`sample_id`); its columns were copied from ENA/SRA and
  BioSample metadata. It is deposit data, not curation.
- **The vocabularies**:
  `registry/controlled_vocabularies/`
  (`tissue.yaml`, `life_stage.yaml`, `sample_unit.yaml`, `stressor.yaml`).
- **Primary sources on the web**:
  - ENA run records:
    `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<BIOPROJECT>&result=read_run&fields=run_accession,sample_accession,library_name,sample_title,sample_alias,library_strategy,library_selection,library_construction_protocol,instrument_model,first_public&format=tsv`
  - BioSample records (every attribute): `https://www.ebi.ac.uk/biosamples/samples/<SAMN… or SAMEA…>`
    (JSON with `Accept: application/json`), or NCBI E-utilities for SAMN accessions.
  - The BioProject description (ENA `https://www.ebi.ac.uk/ena/browser/api/xml/<BIOPROJECT>`).
  - The paper, its methods and its supplementary tables: PubMed / PMC, the
    publisher page, preprints (bioRxiv, HAL, SSRN, Research Square), theses.
    If no paper is given, find it: search by BioProject accession, by the
    authors and title in the BioProject record, and by the design.

## What you must not read

Anything else in the repository: not `registry/studies/*.yaml`, not `docs/`,
not `data/studies/*` other than the sample sheet, not git history, not other
curators' output folders. Those hold an earlier curation, and reading them
would make yours a copy.

Do not modify anything in the repository. Write only to your output folder.

## Output

Write two files to your output folder (given in your task).

### `samples.tsv`

Tab-separated, header row, **one row per row of the sample sheet, in the
same order, every library** (including ones you suspect are bad).

| column | content |
|---|---|
| `sample_id` | exactly as in the sample sheet |
| `animal_id` | an identifier for the animal (or pool) the library came from. Libraries from the same animal (two tissues, two time points, technical replicates) get the **same** value. Use the paper's or deposit's code when there is one; otherwise empty |
| `tissue` | a `tissue.yaml` term |
| `life_stage` | a `life_stage.yaml` term (`unknown` if not stated) |
| `age` | as stated, free text (`10 months`, `48 hpf`), else empty |
| `sex` | `female`, `male`, `hermaphrodite` or `unknown` |
| `ploidy` | `diploid`, `triploid`, `tetraploid` or `unknown` |
| `sample_unit` | a `sample_unit.yaml` term |
| `pool_size` | animals per library when pooled, else empty |
| `population` | short label of the source population, line or family |
| `treatment` | short label of this library's experimental group in the paper's terms (`pH 7.7`, `heat 35 C 6 h`, `control`) |
| `confidence` | `high` (stated per sample in the deposit or paper), `medium` (stated for the study as a whole, applied to every library), `low` (inferred) |
| `sources` | where each non-empty value came from, `field: source` pairs separated by `; ` (`tissue: BioSample SAMN17095178 attribute tissue; age: paper Methods 2.1`) |
| `notes` | anything a reviewer must know about this row |

Use only vocabulary terms. If no term fits, use the closest one (or `other`)
and list the case in `study.json` → `unmappable_terms`. Never invent a term.

### `study.json`

```json
{
  "study_id": "...",
  "publication": {"citation": "...", "doi": "...", "pmid": "...", "how_found": "...",
                  "full_text_read": true, "supplements_read": true},
  "library_kit": {"name": "...", "directional": "directional|non_directional|pbat|unknown",
                  "synthetic_tail": "yes|no|unknown", "source": "..."},
  "design_summary": "2-5 sentences: animals, groups, replicates, what one library is",
  "published_dml_dmr": {"available": true, "where": "Supplementary Table S5 (xlsx)",
                        "what": "DMRs per contrast, DSS, FDR < 0.05", "assembly": "...",
                        "url": "...", "source": "..."},
  "paired_rnaseq": {"bioproject": "...", "relationship": "same_animals|same_design_different_animals|none|unknown",
                    "source": "..."},
  "discrepancies": ["deposit says X, paper says Y (where)"],
  "unmappable_terms": [{"value": "gastrula", "used": "embryo", "why": "..."}],
  "open_questions": ["..."],
  "sources_consulted": ["every URL or document you read"]
}
```

## Rules

1. **Cite every value.** A value with no source is not allowed.
2. **Never guess.** If neither the deposit nor the paper states it, write
   `unknown` (or leave the free-text column empty) and say so in `notes`. An
   inference is allowed only when it is labeled as one in `notes` and
   `confidence` is `low`.
3. **Deposit and paper disagree?** Record both in `discrepancies`, and use
   the value the per-sample evidence supports, saying which in `notes`.
4. **Pooling, repeated sampling and shared animals matter most.** They
   decide how many independent animals a study has. Read the methods for
   them specifically.
5. **The library kit matters.** Look in the deposit's
   `library_construction_protocol` and the paper's methods for the kit
   (e.g. Zymo Pico Methyl-Seq is non-directional and adds a synthetic tail;
   NEB EM-seq and classic TruSeq-style WGBS are directional).
6. If the paper cannot be found or read, say so in `publication`, curate
   from the deposit alone, and lower `confidence` accordingly.
7. Finish with a short reply: the output folder, the publication you used,
   how many rows have each confidence level, and the three most important
   problems you found.
