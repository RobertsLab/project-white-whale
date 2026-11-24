# Replication: VALDIVIESO2025_POMS_RS_EMSEQ

Written by `whale replicate` (src/whale/replication.py) from `claims.yaml` and the study's committed results. Do not edit by hand; edit the claims or re-run.

## Sources

- **P1**: Valdivieso A, Morga B, Degremont L, Mege M, Courtay G, Dorant Y, Escoubas J-M, Gawra J, de Lorgeril J, Mitta G, Cosseau C, Vidal-Dupiol J (2025). DNA methylation landscapes before and after Pacific Oyster Mortality Syndrome are different within and between resistant and susceptible Magallana gigas. Science of the Total Environment 962:178385. doi:10.1016/j.scitotenv.2025.178385 — read: full text (HAL accepted manuscript, hal-04901984); supplementary workbook mmc10.xlsx Tables S1-S18 (S4, S5, S7, S9, S11, S13-S18 tallied); supplementary Figs S3 (global methylation) and S4 (PCA) viewed

## Re-derived from

- Pipeline run `lethal_goodall` on GCF_963853765.1, SNP mask `GCF_963853765.1.snpmask.v1`, 40 QC-passing samples
- Shared differential model: run 2026-11-20

## Caveats

- The paper compares T1 with T0 within each phenotype; the shared model's time jobs pool both phenotypes (with each oyster as a block), so those verdicts are approximate and at most partly.

## Summary

| replicated | partly | not_replicated | inconclusive | not_testable | pending |
|---|---|---|---|---|---|
| 5 | 9 | 6 | 0 | 8 | 0 |

## Claims

| id | where | claim (paraphrased) | verdict | re-derived | rule |
|---|---|---|---|---|---|
| C1 | P1 Results 3.2; Fig. S3 | Genome-wide CpG methylation is low, about 12-13 % per library (individual values roughly 11.5-13.6 % in Fig. S3 across tissues, phenotypes and times); the paper gives no single number in the text, value read from the figure. | **replicated** | 11.31% (mean over 40 samples, all samples; CpG calls pooled per library); paper 12.6% | \|difference\| <= 2.0 points: replicated; <= 5.0: partly |
| C2 | P1 Results 3.2; Discussion 4.1; Fig. S3A | In gills, global methylation is consistently higher after POMS infection (T1) than before (T0), by about 0.54 percentage points. | **replicated** | time T1 11.46% vs T0 11.05% (difference +0.41 points, n = 10 vs 10, permutation p = 0.0065, within tissue_type = G) | same direction and p < 0.05 (or, when group sizes cannot reach 0.05, complete separation): replicated; same direction: partly; opposite: not replicated |
| C3 | P1 Results 3.2; Discussion 4.1; Fig. S3B | In mantle, global methylation is consistently higher at T1 than at T0, by about 0.48 percentage points. | **replicated** | time T1 11.54% vs T0 11.20% (difference +0.35 points, n = 10 vs 10, permutation p = 0.001, within tissue_type = M) | same direction and p < 0.05 (or, when group sizes cannot reach 0.05, complete separation): replicated; same direction: partly; opposite: not replicated |
| C4 | P1 Results 3.2; Fig. S3A | Global methylation in gills is higher in susceptible than resistant oysters (one-way ANOVA across the four phenotype x time groups, F = 4.62, P = 0.005); in Fig. S3A the only group differing significantly is susceptible T1 (~13.2 % vs ~12.4-12.8 % for the others). | **partly** | phenotype_class S 11.41% vs R 11.10% (difference +0.32 points, n = 10 vs 10, permutation p = 0.0524, within tissue_type = G) | same direction and p < 0.05 (or, when group sizes cannot reach 0.05, complete separation): replicated; same direction: partly; opposite: not replicated |
| C5 | P1 Results 3.2; Fig. S3B | Global methylation in mantle is higher in susceptible than resistant oysters (one-way ANOVA, F = 7.49, P = 0.002); in Fig. S3B the difference is driven by susceptible T1 (~12.9-13 %), while at T0 susceptible medians are slightly below resistant ones. | **partly** | phenotype_class S 11.43% vs R 11.31% (difference +0.12 points, n = 10 vs 10, permutation p = 0.351, within tissue_type = M) | same direction and p < 0.05 (or, when group sizes cannot reach 0.05, complete separation): replicated; same direction: partly; opposite: not replicated |
| C6 | P1 Results 3.2; Fig. S4 | In PCAs of CpGs common to all samples (1,886,331 in gills, 1,806,224 in mantle), samples cluster mainly by source population (each oyster's T0 and T1 sit together; population 5, SC18, separates on PC1), PC1 and PC2 each explaining about 8-11 % of variance. | **replicated** | population explains 97% (adjusted R^2) of PC1 (14.43% of variance), the larger of PC1 and PC2; PCA of job time.G | adjusted R^2 >= 0.5 replicated; >= 0.2 partly |
| C7 | P1 Results 3.2; Fig. S4 | Infection (time) has only a subtle effect on the methylation landscape: T0 and T1 samples do not separate in the PCA, nor do phenotypes. | **replicated** | time explains 0% (adjusted R^2) of PC2 (11.87% of variance), the larger of PC1 and PC2; PCA of job time.G | the paper reports no separation: adjusted R^2 < 0.2 replicated; < 0.5 partly |
| C8 | P1 Results 3.3.1; Tables S5-S7; Fig. 2A | In susceptible gills, T1 vs T0 gave 3069 DMRs located in 2014 genes. | **partly** | job time.G, timeT1: gene_body: 1468 genes at FDR < 0.05 (lambda 1.442); promoter: 24 genes at FDR < 0.05 (lambda 0.737); paper 2014 DMGs | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over phenotype_class (paper: S only), so at most partly; test statistics are inflated (lambda 1.442 > 1.3), so the count is not trusted: at most partly |
| C9 | P1 Results 3.3.2; Tables S5, S8, S9; Fig. 2B | In resistant gills, T1 vs T0 gave 2594 DMRs located in 1702 genes (Fig. 2 legend says 1027 genes, apparently a typo; Table S9 lists 1702). | **partly** | job time.G, timeT1: gene_body: 1468 genes at FDR < 0.05 (lambda 1.442); promoter: 24 genes at FDR < 0.05 (lambda 0.737); paper 1702 DMGs | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over phenotype_class (paper: R only), so at most partly; test statistics are inflated (lambda 1.442 > 1.3), so the count is not trusted: at most partly |
| C10 | P1 Results 3.3.4; Tables S5, S10, S11; Fig. 3A | In susceptible mantle, T1 vs T0 gave 2836 DMRs in 1866 genes. | **partly** | job time.M, timeT1: gene_body: 217 genes at FDR < 0.05 (lambda 1.326); promoter: 27 genes at FDR < 0.05 (lambda 0.79); paper 1866 DMGs | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over phenotype_class (paper: S only), so at most partly; test statistics are inflated (lambda 1.326 > 1.3), so the count is not trusted: at most partly |
| C11 | P1 Results 3.3.5; Tables S5, S12, S13; Fig. 3B | In resistant mantle, T1 vs T0 gave 5469 DMRs, the largest response, in 2977 genes (1683 hyper + 1294 hypo; Table S13 lists 2977, Fig. 3 legend says 2997). | **partly** | job time.M, timeT1: gene_body: 217 genes at FDR < 0.05 (lambda 1.326); promoter: 27 genes at FDR < 0.05 (lambda 0.79); paper 2977 DMGs | genes vs genes: within 3x of the paper's count replicated; fewer or more but some: partly; none: not replicated; approximate: our job pools over phenotype_class (paper: R only), so at most partly; test statistics are inflated (lambda 1.326 > 1.3), so the count is not trusted: at most partly |
| C12 | P1 Results 3.3.1; Table S7 | In susceptible gills, DMR genes are split almost evenly after infection: 1014 hypermethylated vs 1000 hypomethylated at T1 (50.3 % hyper). | **partly** | job time.G, timeT1: 1469 higher and 23 lower in T1 (FDR < 0.05, both regions) | majority direction of >= 5 significant genes; approximate: our job pools over phenotype_class (paper: S only), so at most partly |
| C13 | P1 Results 3.3.2; Table S9 | In resistant gills, slightly more DMR genes gain methylation after infection: 923 hyper vs 779 hypo (54.2 % hyper). | **partly** | job time.G, timeT1: 1469 higher and 23 lower in T1 (FDR < 0.05, both regions) | majority direction of >= 5 significant genes; approximate: our job pools over phenotype_class (paper: R only), so at most partly |
| C14 | P1 Results 3.3.4; Table S11 | In susceptible mantle, most DMR genes lose methylation after infection: 1176 hypo vs 690 hyper (63.0 % hypo). | **not_replicated** | job time.M, timeT1: 232 higher and 12 lower in T1 (FDR < 0.05, both regions) | majority direction of >= 5 significant genes; approximate: our job pools over phenotype_class (paper: S only), so at most partly |
| C15 | P1 Results 3.3.5; Table S13 | In resistant mantle, most DMR genes gain methylation after infection: 1683 hyper vs 1294 hypo (56.5 % hyper). | **partly** | job time.M, timeT1: 232 higher and 12 lower in T1 (FDR < 0.05, both regions) | majority direction of >= 5 significant genes; approximate: our job pools over phenotype_class (paper: R only), so at most partly |
| C16 | P1 Results 3.3.3, 3.4; Table S14; Discussion 4.1 | Highlighted immune genes are differentially methylated after infection in susceptible gills (directions from Table S14 gene means, ST1-ST0); MAPKAP1 falls from 35.1 % to 31.7 %. | **not_replicated** | job time.G, timeT1: 0 replicated, 0 partly, 7 not, of 7 tested (3 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over phenotype_class (paper: S only), so at most partly |
| C17 | P1 Results 3.3.3, 3.4; Table S14 | The same highlighted immune genes are differentially methylated after infection in resistant gills (directions RT1-RT0 from Table S14); MAPKAP1 falls sharply from 36.9 % to 5.0 %. | **not_replicated** | job time.G, timeT1: 0 replicated, 0 partly, 7 not, of 7 tested (3 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over phenotype_class (paper: R only), so at most partly |
| C18 | P1 Results 3.3.6, 3.4; Table S14 | Highlighted immune genes are differentially methylated after infection in susceptible mantle (directions ST1-ST0 from Table S14); interferon-induced protein 44 falls from 43.6 % to 12.9 %. | **not_replicated** | job time.M, timeT1: 0 replicated, 0 partly, 7 not, of 7 tested (3 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over phenotype_class (paper: S only), so at most partly |
| C19 | P1 Results 3.3.6, 3.4; Table S14 | Highlighted immune genes are differentially methylated after infection in resistant mantle (directions RT1-RT0 from Table S14); interferon-induced protein 44 falls from 47.5 % to 27.1 %. | **not_replicated** | job time.M, timeT1: 0 replicated, 0 partly, 7 not, of 7 tested (3 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over phenotype_class (paper: R only), so at most partly |
| C20 | P1 Results 3.5.1; Fig. 4A; Table S15 | Before infection, a 151 bp, 10-CpG DMR in the Hedgehog-interacting protein gene body is about 20.0 points less methylated in resistant than susceptible gills (22.6 points lower in resistant T1 vs susceptible T0; only 2.6 points between resistant T0 and T1). | **not_testable** | job phenotype.G.T0, phenotype_classR: 0 replicated, 0 partly, 0 not, of 0 tested (1 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly |
| C21 | P1 Results 3.5.2; Fig. 4B; Table S17 | In resistant gills, a 194 bp, 5-CpG DMR in the THO complex subunit 3 gene body loses methylation after infection (resistant T0 is 2.4-fold and susceptible T0 2.8-fold the resistant T1 level). | **not_replicated** | job time.G, timeT1: 0 replicated, 0 partly, 1 not, of 1 tested (0 not found in the annotation, 0 not tested) | per gene: FDR < 0.05 in the claimed direction replicated, p < 0.05 partly (both corrected for the number of matching genes and regions); claim: most tested genes replicated, any replicated or partly: partly; approximate: our job pools over phenotype_class (paper: R only), so at most partly |
| C22 | P1 Abstract; Discussion 4.1; Conclusion | Resistant and susceptible oysters already carry distinct methylation patterns before exposure (T0), irrespective of population origin; the paper reports no plain R-vs-S DMR count at T0, only the filtered epibiomarker sets (C24). | **not_testable** | phenotype_class R vs S within time T0 is a registered contrast, but the paper gives no count or direction for it; only gene-level example is C20 |  |
| C23 | P1 Abstract; Results 3.4; Fig. S7; Table S14 | 31 immune-related genes carry DMRs after infection in all four tissue x phenotype comparisons (from 222, 185, 181 and 284 immune-related DMR genes in S gills, R gills, S mantle, R mantle), with contrasting directions between phenotypes. | **not_testable** | intersection of four stratified T1-vs-T0 gene lists filtered by the authors' immune-gene annotation; could be approximated as an overlap of our four per-gene time contrasts, but no kind covers overlap counts and the immune flag is author-curated |  |
| C24 | P1 Abstract; Results 3.5.1; Tables 1, S15, S16 | 1998 candidate phenotype epibiomarker DMRs (1204 gills, 1030 in gene bodies; 794 mantle, 718 in gene bodies) meet \|RT0-ST0\| >= 20 %, \|RT1-ST0\| >= 20 % and \|RT0-RT1\| <= 5 %. | **not_testable** | multi-criteria DMR filter on group means, not a single DSS contrast; per-region on paper genome coordinates |  |
| C25 | P1 Abstract; Results 3.5.2; Tables 2, S17, S18 | 164 gene-body candidate site (POMS exposure) epibiomarkers (66 of 68 DMRs in gills, 98 of 128 in mantle) meet \|ST0-RT0\| <= 5 %, \|ST0-RT1\| >= 20 % and \|RT0-RT1\| >= 20 %; note the abstract's 164 counts gene-body DMRs while its 1998 counts all DMRs. | **not_testable** | multi-criteria DMR filter on group means; per-region on paper genome coordinates |  |
| C26 | P1 Results 3.3.3, 3.3.6; Fig. S6C-H | Genes with DMRs after infection overlap between phenotypes: 527 in gills (277 with the same direction: 157 hyper, 120 hypo) and 788 in mantle (384 same direction: 213 hyper, 171 hypo). | **not_testable** | overlap between two stratified time contrasts; no kind for gene-list overlap, though it could be approximated from our per-gene results |  |
| C27 | P1 Abstract; Results 3.3.1-3.3.6; Figs. 2-3 | Genes with DMRs after infection are enriched for immunity, signalling, metabolism and transport GO terms and KEGG pathways (e.g. response to virus and DNA damage response in hypermethylated S-gill genes; mitophagy in hypermethylated R-gill genes). | **not_testable** | GO_MWU / DAVID enrichment against the paper's annotation; enrichment is not part of the re-derived outputs |  |
| C28 | P1 Results 3.1; Fig. 1C; Fig. S2; Tables S1-S3 | The cohabitation challenge worked: survival differed among populations (15 % in SC18/#5 to 85 % in Thau/#1, log-rank P <= 0.001), and at T1 OsHV-1 load was higher in susceptible than resistant tissues (Mann-Whitney W = 23, P = 0.032). | **not_testable** | survival and qPCR viral load, not methylation |  |

## Named genes

**C16**

| gene (as published) | status | re-derived |
|---|---|---|
| bak1 / Bcl-2-like protein antagonist/killer | not_found | no gene with this symbol or description in the annotation |
| birc3 / Baculoviral IAP repeat-containing protein 7 | not_replicated | LOC136271523 (baculoviral IAP repeat-containing protein 7-like), gene_body: diff +3.23 points, p 0.0396, FDR 0.262; best of 5 tests (5 matching genes), so p and FDR are multiplied by 5 for the verdict; opposite to the paper's hypo |
| trim3 / Tripartite motif-containing protein 3 | not_replicated | LOC109620450 (Tripartite motif-containing protein 3), gene_body: diff -1.70 points, p 0.00257, FDR 0.0527; best of 103 tests (103 matching genes), so p and FDR are multiplied by 103 for the verdict; opposite to the paper's hyper |
| trim36 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +2.42 points, p 0.00122, FDR 0.0341; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict |
| trim71 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +2.42 points, p 0.00122, FDR 0.0341; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict; opposite to the paper's hypo |
| uba1 / Ubiquitin-activating enzyme E1 | not_found | no gene with this symbol or description in the annotation |
| itga2 / Integrin alpha-2 | not_found | no gene with this symbol or description in the annotation |
| Interferon alpha-inducible protein 27 | not_replicated | LOC109618716 (interferon alpha-inducible protein 27-like protein 2A), gene_body: diff +1.86 points, p 0.00534, FDR 0.0795; best of 7 tests (7 matching genes), so p and FDR are multiplied by 7 for the verdict; opposite to the paper's hypo |
| Macrophage mannose receptor 1 | not_replicated | LOC105324874 (macrophage mannose receptor 1), gene_body: diff +1.63 points, p 0.087, FDR 0.411; best of 11 tests (11 matching genes), so p and FDR are multiplied by 11 for the verdict; opposite to the paper's hypo |
| MAPKAP1 / Target of rapamycin complex 2 subunit MAPKAP1 | not_replicated | LOC105340981 (target of rapamycin complex 2 subunit MAPKAP1), gene_body: diff +2.44 points, p 0.00207, FDR 0.0462; opposite to the paper's hypo |

**C17**

| gene (as published) | status | re-derived |
|---|---|---|
| bak1 / Bcl-2-like protein antagonist/killer | not_found | no gene with this symbol or description in the annotation |
| birc3 / Baculoviral IAP repeat-containing protein 7 | not_replicated | LOC136271523 (baculoviral IAP repeat-containing protein 7-like), gene_body: diff +3.23 points, p 0.0396, FDR 0.262; best of 5 tests (5 matching genes), so p and FDR are multiplied by 5 for the verdict; opposite to the paper's hypo |
| trim3 / Tripartite motif-containing protein 3 | not_replicated | LOC109620450 (Tripartite motif-containing protein 3), gene_body: diff -1.70 points, p 0.00257, FDR 0.0527; best of 103 tests (103 matching genes), so p and FDR are multiplied by 103 for the verdict |
| trim36 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +2.42 points, p 0.00122, FDR 0.0341; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict |
| trim71 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +2.42 points, p 0.00122, FDR 0.0341; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict |
| uba1 / Ubiquitin-activating enzyme E1 | not_found | no gene with this symbol or description in the annotation |
| itga2 / Integrin alpha-2 | not_found | no gene with this symbol or description in the annotation |
| Interferon alpha-inducible protein 27 | not_replicated | LOC109618716 (interferon alpha-inducible protein 27-like protein 2A), gene_body: diff +1.86 points, p 0.00534, FDR 0.0795; best of 7 tests (7 matching genes), so p and FDR are multiplied by 7 for the verdict; opposite to the paper's hypo |
| Macrophage mannose receptor 1 | not_replicated | LOC105324874 (macrophage mannose receptor 1), gene_body: diff +1.63 points, p 0.087, FDR 0.411; best of 11 tests (11 matching genes), so p and FDR are multiplied by 11 for the verdict |
| MAPKAP1 / Target of rapamycin complex 2 subunit MAPKAP1 | not_replicated | LOC105340981 (target of rapamycin complex 2 subunit MAPKAP1), gene_body: diff +2.44 points, p 0.00207, FDR 0.0462; opposite to the paper's hypo |

**C18**

| gene (as published) | status | re-derived |
|---|---|---|
| bak1 / Bcl-2-like protein antagonist/killer | not_found | no gene with this symbol or description in the annotation |
| birc3 / Baculoviral IAP repeat-containing protein 7 | not_replicated | LOC117691872 (baculoviral IAP repeat-containing protein 7), gene_body: diff +2.78 points, p 0.0558, FDR 0.468; best of 5 tests (5 matching genes), so p and FDR are multiplied by 5 for the verdict; opposite to the paper's hypo |
| trim3 / Tripartite motif-containing protein 3 | not_replicated | LOC117682418 (tripartite motif-containing protein 3-like), gene_body: diff +0.64 points, p 0.00388, FDR 0.155; best of 103 tests (103 matching genes), so p and FDR are multiplied by 103 for the verdict |
| trim36 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +3.31 points, p 0.00024, FDR 0.0427; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict; opposite to the paper's hypo |
| trim71 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +3.31 points, p 0.00024, FDR 0.0427; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict; opposite to the paper's hypo |
| uba1 / Ubiquitin-activating enzyme E1 | not_found | no gene with this symbol or description in the annotation |
| itga2 / Integrin alpha-2 | not_found | no gene with this symbol or description in the annotation |
| Interferon alpha-inducible protein 27 | not_replicated | LOC109618716 (interferon alpha-inducible protein 27-like protein 2A), gene_body: diff +1.12 points, p 0.166, FDR 0.674; best of 7 tests (7 matching genes), so p and FDR are multiplied by 7 for the verdict; opposite to the paper's hypo |
| Macrophage mannose receptor 1 | not_replicated | LOC105324874 (macrophage mannose receptor 1), gene_body: diff +1.67 points, p 0.0928, FDR 0.561; best of 11 tests (11 matching genes), so p and FDR are multiplied by 11 for the verdict; opposite to the paper's hypo |
| Interferon-induced protein 44 | not_replicated | LOC136274310 (interferon-induced protein 44-like), gene_body: diff -3.72 points, p 0.00865, FDR 0.214; best of 36 tests (36 matching genes), so p and FDR are multiplied by 36 for the verdict |

**C19**

| gene (as published) | status | re-derived |
|---|---|---|
| bak1 / Bcl-2-like protein antagonist/killer | not_found | no gene with this symbol or description in the annotation |
| birc3 / Baculoviral IAP repeat-containing protein 7 | not_replicated | LOC117691872 (baculoviral IAP repeat-containing protein 7), gene_body: diff +2.78 points, p 0.0558, FDR 0.468; best of 5 tests (5 matching genes), so p and FDR are multiplied by 5 for the verdict |
| trim3 / Tripartite motif-containing protein 3 | not_replicated | LOC117682418 (tripartite motif-containing protein 3-like), gene_body: diff +0.64 points, p 0.00388, FDR 0.155; best of 103 tests (103 matching genes), so p and FDR are multiplied by 103 for the verdict |
| trim36 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +3.31 points, p 0.00024, FDR 0.0427; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict; opposite to the paper's hypo |
| trim71 / Tripartite motif-containing protein 2 (table description) | not_replicated | LOC117682723 (tripartite motif-containing protein 2), gene_body: diff +3.31 points, p 0.00024, FDR 0.0427; best of 247 tests (247 matching genes), so p and FDR are multiplied by 247 for the verdict; opposite to the paper's hypo |
| uba1 / Ubiquitin-activating enzyme E1 | not_found | no gene with this symbol or description in the annotation |
| itga2 / Integrin alpha-2 | not_found | no gene with this symbol or description in the annotation |
| Interferon alpha-inducible protein 27 | not_replicated | LOC109618716 (interferon alpha-inducible protein 27-like protein 2A), gene_body: diff +1.12 points, p 0.166, FDR 0.674; best of 7 tests (7 matching genes), so p and FDR are multiplied by 7 for the verdict |
| Macrophage mannose receptor 1 | not_replicated | LOC105324874 (macrophage mannose receptor 1), gene_body: diff +1.67 points, p 0.0928, FDR 0.561; best of 11 tests (11 matching genes), so p and FDR are multiplied by 11 for the verdict; opposite to the paper's hypo |
| Interferon-induced protein 44 | not_replicated | LOC136274310 (interferon-induced protein 44-like), gene_body: diff -3.72 points, p 0.00865, FDR 0.214; best of 36 tests (36 matching genes), so p and FDR are multiplied by 36 for the verdict |

**C20**

| gene (as published) | status | re-derived |
|---|---|---|
| Hedgehog-interacting protein | not_found | no gene with this symbol or description in the annotation |

**C21**

| gene (as published) | status | re-derived |
|---|---|---|
| THO complex subunit 3 | not_replicated | LOC105345709 (THO complex subunit 3), gene_body: diff +1.50 points, p 0.0505, FDR 0.304; best of 2 tests (2 matching genes), so p and FDR are multiplied by 2 for the verdict; opposite to the paper's hypo |


## Reading this

Differential results here are per gene (gene body or promoter, summed calls, DSS beta-binomial model, FDR < 0.05), on GCF_963853765.1. Papers count CpGs or regions on other assemblies, so counts are compared by presence, and named genes are matched by symbol or product description. See docs/replication.md.
