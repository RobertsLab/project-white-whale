# Replication: WANG2023_TIDAL_TRANSGEN_WGBS

Written by `whale replicate` (src/whale/replication.py) from `claims.yaml` and the study's committed results. Do not edit by hand; edit the claims or re-run.

## Sources

- **P1**: Wang X, Cong R, Li A, Wang W, Zhang G, Li L (2023). Transgenerational effects of intertidal environment on physiological phenotypes and DNA methylation in Pacific oysters. Science of the Total Environment 871:162112. doi:10.1016/j.scitotenv.2023.162112 — read: abstract only (PubMed record, PMID 36764539). Full text is closed access: no PMC copy, Europe PMC lists only the subscription DOI link, Semantic Scholar has no open-access PDF. Methods, figures, tables and supplements were not seen, so DMG calling method, thresholds and region (gene body vs promoter) are unknown.

## Re-derived from

- Pipeline run `elated_pesquet` on GCF_963853765.1, SNP mask `GCF_963853765.1.snpmask.v1`, 24 QC-passing samples
- Shared differential model: run 2026-11-20

## Caveats

- These results come from run elated_pesquet, which did not clip the synthetic tail on read 2; global methylation levels are about 1.5 points too high and are not compared here. Re-run with the clipping in PR #61.
- Only the abstract could be read (closed access), so the claims are few and the paper's methods, thresholds and gene lists are unknown.

## Summary

| replicated | partly | not_replicated | inconclusive | not_testable | pending |
|---|---|---|---|---|---|
| 0 | 2 | 0 | 2 | 8 | 0 |

## Claims

| id | where | claim (paraphrased) | verdict | re-derived | rule |
|---|---|---|---|---|---|
| C1 | P1 Abstract | Intertidal and subtidal wild (F0) oysters differ in gill DNA methylation at roughly 3,864 differentially methylated genes; the number is not stated, but follows from 1655 DMGs being 42.83% of the F0 set. | **partly** | job main, tidal_levelIn: gene_body: 3645 genes at FDR < 0.05 (lambda 1.937); promoter: 451 genes at FDR < 0.05 (lambda 1.266); paper 3864 DMGs | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over generation (paper: F0 only), so at most partly; test statistics are inflated (lambda 1.937 > 1.3), so the count is not trusted: at most partly |
| C2 | P1 Abstract | The intertidal vs subtidal methylation divergence is heritable, so In vs Sub DMGs are still present in the F1 generation bred in the subtidal area; no F1 count is given in the abstract. | **inconclusive** | job main, tidal_levelIn: gene_body: 3645 genes at FDR < 0.05 (lambda 1.937); promoter: 451 genes at FDR < 0.05 (lambda 1.266) | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over generation (paper: F1 only), so at most partly; test statistics are inflated (lambda 1.937 > 1.3), so the count is not trusted: at most partly |
| C3 | P1 Abstract | Intertidal vs subtidal DMGs persist to F2, at least 1655 genes (those shared with the F0 set), so the F2 In vs Sub contrast should give at least about 1655 DMGs under the paper's method. | **partly** | job main, tidal_levelIn: gene_body: 3645 genes at FDR < 0.05 (lambda 1.937); promoter: 451 genes at FDR < 0.05 (lambda 1.266); paper 1655 DMGs | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over generation (paper: F2 only), so at most partly; test statistics are inflated (lambda 1.937 > 1.3), so the count is not trusted: at most partly |
| C4 | P1 Abstract | 1655 of the F0 intertidal vs subtidal DMGs (42.83%) are again differentially methylated in F2, i.e. the overlap of the F0 and F2 In vs Sub DMG lists is 1655 genes. | **not_testable** | Requires overlapping per-generation In vs Sub gene lists (tidal_level stratified by generation); no claim kind covers list overlap. Could be approximated if stratified contrasts exist. |  |
| C5 | P1 Abstract | Laboratory heat shock changes gill DNA methylation in intertidal F0 oysters (heat vs control DMGs exist); count not stated in the abstract. | **inconclusive** | job main, treatmentH: gene_body: 8 genes at FDR < 0.05 (lambda 0.723); promoter: 2 genes at FDR < 0.05 (lambda 0.841) | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over tidal_level (paper: In only), generation (paper: F0 only), so at most partly |
| C6 | P1 Abstract | About 2,171 thermal-response DMGs are specific to the intertidal F0 group (not found in subtidal F0); the number follows from 320 being 14.74% of that set. | **not_testable** | 'Specific' is a set difference between H vs C DMGs in In-F0 and in Sub-F0, not a single registered contrast; the count 2171 = 320 / 0.1474 is derived. |  |
| C7 | P1 Abstract | 320 (14.74%) of the intertidal-F0-specific thermal-response DMGs are also found among thermal-response DMGs in F1 and F2. | **not_testable** | Needs overlap of heat-response DMG lists across tidal levels and generations; no single contrast. |  |
| C8 | P1 Abstract | Transgenerational plasticity alters physiological phenotypes related to oxidative status and energy in unstressed subtidal-reared offspring for at least two generations. | **not_testable** | Physiological assays (oxidation and energy markers), not methylation data. |  |
| C9 | P1 Abstract | Intertidal parental experience alters the thermal plasticity of physiological phenotypes within and across generations. | **not_testable** | Physiological phenotype data, not methylation. |  |
| C10 | P1 Abstract | The 1655 retained F0-to-F2 DMGs are enriched for macromolecule metabolism, cellular response to stress, positive regulation of molecular function, and especially fatty acid metabolism. | **not_testable** | GO/functional enrichment of a cross-generation overlap set; no named genes given in the abstract. |  |
| C11 | P1 Abstract | The 320 transgenerationally retained thermal-response DMGs map to carbohydrate, lipid and energy metabolism, signal transduction, and organismal immune system pathways. | **not_testable** | KEGG pathway enrichment of an overlap set; no named genes in the abstract. |  |
| C12 | P1 Abstract | DNA methylation is an epigenetic mechanism mediating transgenerational plasticity for stress adaptation in oysters. | **not_testable** | Interpretive conclusion linking methylation to phenotype. |  |

## Reading this

Differential results here are per gene (gene body or promoter, summed calls, DSS beta-binomial model, FDR < 0.05), on GCF_963853765.1. Papers count CpGs or regions on other assemblies, so counts are compared by presence, and named genes are matched by symbol or product description. See docs/replication.md.
