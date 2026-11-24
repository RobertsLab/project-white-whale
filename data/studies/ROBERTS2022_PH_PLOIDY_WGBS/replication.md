# Replication: ROBERTS2022_PH_PLOIDY_WGBS

Written by `whale replicate` (src/whale/replication.py) from `claims.yaml` and the study's committed results. Do not edit by hand; edit the claims or re-run.

## Sources

- **P1**: Venkataraman YR (2021). Ocean acidification influences on physiology and epigenetics in the Pacific oyster (Crassostrea gigas). PhD dissertation, University of Washington. Chapter 4: Polyploidy and environmental stress have distinct impacts on Pacific oyster (Crassostrea gigas) ctenidia methylomes (pp. 86-112). — read: full text of chapter 4 (abstract, methods, results, discussion, Table 4.1, Figs 4.1-4.4) and the Appendix A.4 data-accessibility list, from the ResearchWorks PDF
- **R1**: RobertsLab/project-oyster-oa GitHub repository, analyses/Haws_04-DSS/DML, Haws_07-DML-characterization, Haws_09-functional-enrichment (authors' own output cited in Appendix A.4). — read: DSS DML tables (DML-pH-DSS.csv, DML-ploidy-DSS.csv, DML-ploidypH-DSS.csv) and enriched-GO-with-transcript tables; used only to make P1's claims precise (gene identities behind Table 4.1). NCBI Gene esummary used to resolve LOC gene IDs to product names.

## Re-derived from

- Pipeline run `curious_mestorf` on GCF_963853765.1, SNP mask `GCF_963853765.1.snpmask.v1`, 24 QC-passing samples
- Shared differential model: run 2026-11-20

## Caveats

- These results come from run curious_mestorf, which aligned these non-directional (Zymo Pico Methyl-Seq) libraries in directional mode and so used about half of the read pairs (30% mapping). Calls are unbiased but at half depth; re-run with --non_directional (PR #61) and re-derive before treating the differential verdicts as final.
- The paper is a dissertation chapter; the authors' own DML tables (RobertsLab/project-oyster-oa) label the pH and ploidy sets the other way round from the chapter's text. The claims follow the chapter.

## Summary

| replicated | partly | not_replicated | inconclusive | not_testable | pending |
|---|---|---|---|---|---|
| 7 | 2 | 2 | 0 | 13 | 0 |

## Claims

| id | where | claim (paraphrased) | verdict | re-derived | rule |
|---|---|---|---|---|---|
| C1 | P1 Results 4.4.1 | Across all 24 samples, 11.6% of cytosines in CpG context were called methylated (Bismark, pooled over samples). | **replicated** | 10.94% (mean over 24 samples, all samples; CpG calls pooled per library); paper 11.6% | \|difference\| <= 2.0 points: replicated; <= 5.0: partly |
| C2 | P1 Abstract; Results 4.4.3; Discussion 4.5.1 | There was no global methylation difference between diploid and triploid oysters; both were predominantly lowly methylated (no difference value or test statistic reported). | **replicated** | ploidy diploid 10.81% vs triploid 11.07% (difference +0.26 points, permutation p = 0.142) | the paper reports no difference: p >= 0.05 replicated; inconclusive when the group sizes cannot reach p < 0.05 |
| C3 | P1 Abstract | There were no global methylation differences between samples, including between low (7.7) and high (8.2) pH treatments (no values reported). | **replicated** | ph 7.7 10.90% vs 8.2 10.98% (difference +0.09 points, permutation p = 0.649) | the paper reports no difference: p >= 0.05 replicated; inconclusive when the group sizes cannot reach p < 0.05 |
| C4 | P1 Results 4.4.3; Fig. 4.1E | A global PCA of CpG percent methylation (PC1 11.0%, PC2 6.63% of variance) did not separate samples by ploidy. | **replicated** | ploidy explains 4% (adjusted R^2) of PC1 (12.41% of variance); PCA of job main | the paper reports no separation: adjusted R^2 < 0.2 replicated; < 0.5 partly |
| C5 | P1 Results 4.4.3; Fig. 4.1E | The same global PCA did not separate samples by pH treatment. | **replicated** | ph explains 0% (adjusted R^2) of PC1 (12.41% of variance); PCA of job main | the paper reports no separation: adjusted R^2 < 0.2 replicated; < 0.5 partly |
| C6 | P1 Results 4.4.3; Fig. 4.1C-D | Of 10,497,320 CpGs with >=5x coverage in every sample, 9.4% were highly methylated (mean >=50%). | **not_testable** | our tables are per gene, not per CpG |  |
| C7 | P1 Results 4.4.3; Fig. 4.2 | Methylation is concentrated in genes: 94.3% of highly methylated CpGs lay in genes (53.2% intronic, 35.7% CDS), a distribution significantly different from all 5x CpGs (chi-squared). | **replicated** | gene_body 14.73% vs genome 10.94% (CpG calls pooled over all genes, mean over samples) | direction only |
| C8 | P1 Results 4.4.3; Fig. 4.2 | Upstream 1 kb flanking regions held only 0.5% of highly methylated CpGs versus 94.3% in genes, implying low promoter-proximal methylation relative to gene bodies. | **replicated** | gene_body 14.73% vs promoter 3.66% (CpG calls pooled over all genes, mean over samples) | direction only |
| C9 | P1 Abstract; Results 4.4.4; Fig. 4.3A | DSS multifactor model identified 154 ploidy-DML (ploidyT coefficient) out of 12,794,074 CpGs with >=1x coverage. | **partly** | job main, ploidytriploid: gene_body: 7 genes at FDR < 0.05 (lambda 1.03); promoter: 1 genes at FDR < 0.05 (lambda 0.862); paper 154 DMCs | units differ (the paper counts DMCs, we count genes), so only presence is compared: both some partly, both none replicated, otherwise not replicated; approximate: our job pools over ph (paper: 8.2 only), so at most partly |
| C10 | P1 Abstract; Results 4.4.4; Fig. 4.3B | DSS multifactor model identified 178 pH-DML (pHL coefficient). | **partly** | job main, ph7.7: gene_body: 0 genes at FDR < 0.05 (lambda 1.006); promoter: 3 genes at FDR < 0.05 (lambda 0.883); paper 178 DMCs | units differ (the paper counts DMCs, we count genes), so only presence is compared: both some partly, both none replicated, otherwise not replicated; approximate: our job pools over ploidy (paper: diploid only), so at most partly |
| C11 | P1 Results 4.4.4; Fig. 4.3C | The ploidy x pH interaction term gave 53 interaction-DML. | **not_testable** | Interaction contrast is not a registered contrast (registry has ploidy and ph main factors only); also CpG-level. |  |
| C12 | P1 Results 4.4.4; Fig. 4.4 | Ploidy-DML were mostly genic: 145 of 154 (94.2%) fell in 109 genes, chiefly introns (74.0%), then CDS (16.9%) and UTR (5.2%). | **not_testable** | Counts genes containing a significant CpG, not genes tested as units; intron/CDS split also needs other regions. |  |
| C13 | P1 Results 4.4.4; Fig. 4.4 | pH-DML were mostly genic: 123 of 178 (69.1%) fell in 94 genes, chiefly introns (58.4%), with 48.3% in transposable elements and 15.2% intergenic. | **not_testable** | Per-CpG locations in introns/TEs/intergenic are not available; genes-with-DML is not a per-gene test. |  |
| C14 | P1 Results 4.4.4 | DML sets overlapped little: 21 loci shared by pH and ploidy, 11 by pH and interaction, 17 by ploidy and interaction, and 25 in all three. | **not_testable** | Locus-level overlaps. |  |
| C15 | P1 Results 4.4.4; Fig. 4.3D | Chromosome 9 carried the most DML per 10,000 CpGs despite having the fewest genes. | **not_testable** | Per-chromosome DML density on the Roslin assembly. |  |
| C16 | P1 Results 4.4.4; Discussion 4.5.1 | Few DML overlapped C/T SNPs (5 pH-DML, 3 ploidy-DML, 1 interaction-DML), so DML were not SNP artefacts. | **not_testable** | Locus-level SNP overlap. |  |
| C17 | P1 Results 4.4.5; Table 4.1 | Genes with ploidy-DML were enriched for protein dephosphorylation and cell-cell adhesion (3 BP terms, 2 genes each); the dephosphorylation gene carries 3 ploidy-DML. | **not_replicated** | job main, ploidytriploid: 0 replicated, 0 partly, 1 not, of 1 tested (0 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over ph (paper: 8.2 only), so at most partly |
| C18 | P1 Abstract; Results 4.4.5; Table 4.1; Discussion 4.5.2 | Genes with pH-DML were enriched for apoptosis (10 BP terms) and for myosin/cytoskeletal binding, ubiquitin ligase activity and zinc ion binding (8 MF terms); the genes behind these terms are listed. | **not_replicated** | job main, ph7.7: 0 replicated, 0 partly, 16 not, of 16 tested (1 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over ploidy (paper: diploid only), so at most partly |
| C19 | P1 Results 4.4.5; Discussion 4.5.2; Conclusion | Interaction-DML genes showed no significantly enriched GO terms, so ploidy and pH together showed no clear functional signature. | **not_testable** | Interaction contrast not registered; GO enrichment not part of re-derived results. |  |
| C20 | P1 Abstract; Discussion 4.5.1 | Functions of ploidy-DML genes do not resemble those reported for naturally polyploid organisms (e.g. transposon silencing, transcriptional control). | **not_testable** | Interpretive comparison with other taxa. |  |
| C21 | P1 Abstract; Results 4.4.4-4.4.5 | Both ploidy- and pH-DML were enriched in genes with multiple transcripts, suggesting regulation via alternative splicing. | **not_testable** | Transcript-structure/splicing inference; no RNA-seq. |  |
| C22 | P1 Results 4.4.3 | Pairwise Pearson correlations of CpG methylation between samples were uniformly high (0.87-0.93). | **not_testable** | Computed on per-CpG methylKit profiles (>=5x, >=8 samples per group). |  |
| C23 | P1 Results 4.4.1; Fig. 4.1A-B | Coverage was similarly distributed in diploids and triploids (of 12.79 M CpGs with >=1x, 5.40 M had higher mean coverage in triploids). | **not_testable** | Per-CpG coverage comparison; per-library CpGs covered is available but the claim is about the coverage distribution. |  |
| C24 | P1 Results 4.4.1 | 61.0% of trimmed read pairs aligned uniquely to the Roslin genome plus mitochondrion (7.4 x 10^8 of ~2.4 x 10^9 reads). | **not_testable** | Alignment metric on a different genome assembly and aligner settings (non-directional, score_min L,0,-0.6). |  |

## Named genes

**C17**

| gene (as published) | status | re-derived |
|---|---|---|
| receptor-type tyrosine-protein phosphatase epsilon | not_replicated | LOC105335144 (receptor-type tyrosine-protein phosphatase epsilon-like), gene_body: diff -0.65 points, p 0.015, FDR 0.732; best of 41 tests (41 matching genes), so p and FDR are multiplied by 41 for the verdict |

**C18**

| gene (as published) | status | re-derived |
|---|---|---|
| tripartite motif-containing protein 5 | not_replicated | LOC117682664 (tripartite motif-containing protein 5-like), gene_body: diff -2.73 points, p 0.144, FDR 0.97; best of 8 tests (8 matching genes), so p and FDR are multiplied by 8 for the verdict |
| uncharacterized LOC105330472 | not_replicated | LOC105330472 (uncharacterized LOC105330472), gene_body: diff -0.58 points, p 0.314, FDR 0.996 |
| uncharacterized LOC117683735 | not_found | no gene with this symbol or description in the annotation |
| mitogen-activated protein kinase kinase kinase 15 | not_replicated | LOC105338968 (mitogen-activated protein kinase kinase kinase 15), gene_body: diff -0.10 points, p 0.958, FDR 0.997 |
| mitogen-activated protein kinase kinase kinase 4 | not_replicated | LOC105334041 (mitogen-activated protein kinase kinase kinase 4-like), gene_body: diff -5.41 points, p 0.348, FDR 0.997; best of 10 tests (10 matching genes), so p and FDR are multiplied by 10 for the verdict |
| serine protease HTRA2, mitochondrial | not_replicated | LOC105347711 (serine protease HTRA2, mitochondrial), gene_body: diff -2.38 points, p 0.12, FDR 0.954 |
| apoptosis-inducing factor 1, mitochondrial | not_replicated | LOC105348828 (apoptosis-inducing factor 1, mitochondrial), gene_body: diff +2.40 points, p 0.429, FDR 0.997 |
| capping protein inhibiting regulator of actin dynamics | not_replicated | LOC117684293 (capping protein inhibiting regulator of actin dynamics), gene_body: diff -2.44 points, p 0.273, FDR 0.992 |
| FERM domain-containing protein 5 | not_replicated | LOC105322131 (FERM domain-containing protein 5), gene_body: diff -0.15 points, p 0.897, FDR 0.997 |
| spectrin beta chain | not_replicated | LOC105340614 (spectrin beta chain), gene_body: diff -0.57 points, p 0.581, FDR 0.997; best of 2 tests (2 matching genes), so p and FDR are multiplied by 2 for the verdict |
| E3 ubiquitin-protein ligase RNF213-like | not_replicated | LOC136271845 (E3 ubiquitin-protein ligase rnf213-alpha-like), gene_body: diff -7.24 points, p 0.021, FDR 0.799; best of 7 tests (7 matching genes), so p and FDR are multiplied by 7 for the verdict |
| E3 ubiquitin-protein ligase TRIM71-like | not_replicated | LOC105332187 (E3 ubiquitin-protein ligase TRIM71), gene_body: diff -7.74 points, p 0.00371, FDR 0.524; best of 260 tests (260 matching genes), so p and FDR are multiplied by 260 for the verdict |
| prostaglandin reductase-3 | not_replicated | LOC105328838 (prostaglandin reductase-3), gene_body: diff -0.39 points, p 0.885, FDR 0.997 |
| probable nuclear hormone receptor HR3 | not_replicated | LOC105335191 (probable nuclear hormone receptor HR3), gene_body: diff -1.55 points, p 0.253, FDR 0.992 |
| histone deacetylase 6 | not_replicated | LOC105348104 (histone deacetylase 6), gene_body: diff +2.90 points, p 0.585, FDR 0.997; best of 2 tests (2 matching genes), so p and FDR are multiplied by 2 for the verdict |
| S phase cyclin A-associated protein in the endoplasmic reticulum | not_replicated | LOC105348300 (S phase cyclin A-associated protein in the endoplasmic reticulum), gene_body: diff -0.10 points, p 0.954, FDR 0.997 |
| uncharacterized LOC117687382 | not_replicated | LOC117687382 (uncharacterized LOC117687382), gene_body: diff +0.58 points, p 0.477, FDR 0.997 |


## Reading this

Differential results here are per gene (gene body or promoter, summed calls, DSS beta-binomial model, FDR < 0.05), on GCF_963853765.1. Papers count CpGs or regions on other assemblies, so counts are compared by presence, and named genes are matched by symbol or product description. See docs/replication.md.
