# Replication: WANG2014_FIRST_METHYLOMES_WGBS

Written by `whale replicate` (src/whale/replication.py) from `claims.yaml` and the study's committed results. Do not edit by hand; edit the claims or re-run.

## Sources

- **P1**: Wang X, Li Q, Lian J, Li L, Jin L, Cai H, Xu F, Qi H, Zhang L, Wu F, Meng J, Que H, Fang X, Guo X, Zhang G (2014). Genome-wide and single-base resolution DNA methylomes of the Pacific oyster Crassostrea gigas provide insight into the evolution of invertebrate CpG methylation. BMC Genomics 15:1119. doi:10.1186/1471-2164-15-1119 — read: full text (PMC4378273, via PubMed Central) and Additional file 1 (Figs S1-S5, Tables S1-S2, text extracted from the publisher PDF); Additional file 2 (GO tables S3-S5) not read

## Re-derived from

- Pipeline run `scruffy_colden` on GCF_963853765.1, SNP mask `GCF_963853765.1.snpmask.v1`, 2 QC-passing samples
- Shared differential model: not run yet

## Summary

| replicated | partly | not_replicated | inconclusive | not_testable | pending |
|---|---|---|---|---|---|
| 2 | 1 | 0 | 2 | 18 | 0 |

## Claims

| id | where | claim (paraphrased) | verdict | re-derived | rule |
|---|---|---|---|---|---|
| C1 | P1 Results, 'Genome-wide methylation level...'; Table 1 | In the inbred mantle methylome, 15.96% of covered CpGs (>=2x) were called methylated (binomial test, FDR < 0.01). | **not_testable** | our tables are per gene, not per CpG |  |
| C2 | P1 Results, 'Genome-wide methylation level...'; Table 1 | In the wild mantle methylome, 16.25% of covered CpGs (>=2x) were called methylated (binomial test, FDR < 0.01); the paper summarises both animals as up to ~16%. | **not_testable** | our tables are per gene, not per CpG |  |
| C3 | P1 Results, 'Genome-wide methylation level...'; Table 1 | About 2% of all covered genomic cytosines were methylated (3,055,909 mCs = 1.95% in the inbred animal; 2,391,673 mCs = 1.96% in the wild animal), an order of magnitude above insects (0.11-0.15%) and close to sea anemone, sea squirt and human. | **not_testable** | Requires per-cytosine significance calls in all contexts; the cross-species comparison uses other species' data. |  |
| C4 | P1 Abstract; Results; Table 1 | More than 99% of methylated cytosines were in CpG context (99.56% inbred, 99.69% wild), with CHG at 0.08/0.06% and CHH at 0.36/0.25% of mCs. | **not_testable** | Stated as a share of significant mCs by context. Our per-library non-CpG methylation % gives an indirect check: it should be near the ~0.4% bisulfite error rate while CpG methylation is far higher. |  |
| C5 | P1 Table 1 (values only; not an author conclusion) | The wild animal had a marginally higher share of methylated CpGs than the inbred animal (16.25% vs 15.96%, 0.29 points); the authors present the two as essentially equal. | **not_testable** | 1 vs 1 samples | needs >= 2 samples per group |
| C6 | P1 Results; Additional file 1 Fig. S2 | Most methylated CpGs were highly methylated (per-site methylation level distribution of mCGs skewed toward 1), i.e. CpG methylation is largely bimodal. | **not_testable** | Per-CpG level distribution; only per-gene summed calls are available. |  |
| C7 | P1 Results, 'Accumulation of mCGs...'; Fig. 2A | The oyster genome shows a mosaic methylation pattern, with methylated stretches interspersed among long unmethylated regions (mCG density in 1-kb windows along scaffold172). | **not_testable** | Needs windowed mCG density along the authors' scaffolds. |  |
| C8 | P1 Results, 'Accumulation of mCGs...'; Fig. 2B; Conclusions | Gene bodies are the most heavily methylated genomic element, well above the unannotated intergenic background ('other' category). | **replicated** | gene_body 14.75% vs genome 10.88% (CpG calls pooled over all genes, mean over samples) | direction only |
| C9 | P1 Results, 'Gene body methylation, gene length and gene expression'; Fig. 3A; Conclusions | Methylation is low upstream of the TSS and rises sharply after it, so promoters are generally unmethylated relative to gene bodies ('general lack of promoter methylation'). | **replicated** | gene_body 14.75% vs promoter 3.53% (CpG calls pooled over all genes, mean over samples) | direction only |
| C10 | P1 Table 1 (values only; not an author conclusion) | Mean gene-body methylation level was 0.20 in the inbred animal and 0.18 in the wild animal (genes >=70% covered, CpGs >=5x), so inbred slightly higher by 2 points. | **not_testable** | 1 vs 1 samples | needs >= 2 samples per group |
| C11 | P1 Results; Table 1 | Overall gene-body methylation level in oyster is ~0.19 (0.20 inbred, 0.18 wild; mean of per-CpG levels), an order of magnitude higher than insects (0.01-0.02) and close to sea anemone (0.16). | **partly** | 16.66% (mean gene_body methylation of 31889 genes with >= 5 covered CpGs in each of 2 samples, all samples); paper 19.0% | \|difference\| <= 2.0 points: replicated; <= 5.0: partly |
| C12 | P1 Table 1 | 54% of genes in the inbred animal were methylated (>=2 significant mCGs; only genes >=70% covered). | **inconclusive** | 35.6% of 32143 genes have mean gene_body methylation >= 10%, 15.4% >= 50% (1 samples, origin = inbred); paper 54.0% | the paper's definition has no methylation cut-off to compare with |
| C13 | P1 Table 1 | 58% of genes in the wild animal were methylated (>=2 significant mCGs; only genes >=70% covered). | **inconclusive** | 35.8% of 32263 genes have mean gene_body methylation >= 10%, 14.3% >= 50% (1 samples, origin = wild); paper 58.0% | the paper's definition has no methylation cut-off to compare with |
| C14 | P1 Results; Fig. 2B; Conclusions | Within genes, coding exons (CDS) are the most frequently methylated element, above UTRs and introns; tRNA loci are also often methylated, snRNA less, rRNA rarely. | **not_testable** | Needs CDS/UTR/intron and ncRNA-locus matrices. |  |
| C15 | P1 Results; Fig. 2C; Abstract | Repeats as a whole are methylated about twice the intergenic background; SINEs are the most methylated repeat class (near genic levels), while LINEs and LTRs are near background. | **not_testable** | Needs repeat-class annotation and per-repeat methylation. |  |
| C16 | P1 Results; Fig. 2D; Additional file 1 Fig. S3; Abstract | Repeat methylation correlates negatively with repeat divergence (Spearman r = -0.15 inbred, -0.16 wild), and methylated repeats are less diverged than unmethylated ones (Wilcoxon p < 1e-15): young repeats are preferential targets. | **not_testable** | Needs per-repeat methylation and divergence from consensus. |  |
| C17 | P1 Results; Fig. 3A | Along genes, methylation plateaus across the body and falls slowly back to background after the TTS, with no 5'-end bias of the kind seen in insects. | **not_testable** | Needs a binned metagene profile (100-bp bins around TSS/TTS). |  |
| C18 | P1 Results; Fig. 3B; Additional file 1 Fig. S4; Abstract | Gene-body methylation correlates positively with gene length (Spearman r = 0.36 inbred, 0.39 wild), and methylated genes are longer than unmethylated genes (Wilcoxon p < 1e-15). | **not_testable** | No check kind for length correlations; could be tested by joining our per-gene gene-body methylation to gene lengths from the GCF_963853765.1 annotation. |  |
| C19 | P1 Results; Fig. 3A, 3C; Abstract | Gene-body methylation correlates positively with expression (Spearman r = 0.46 over all genes, 0.50 without the top 20% expressed), with moderately expressed genes the most methylated. | **not_testable** | Uses mantle RNA-seq from the wild animal. |  |
| C20 | P1 Results; Fig. 3A | No genome-wide negative relationship between promoter methylation and expression was seen. | **not_testable** | Needs expression data. |  |
| C21 | P1 Results; Fig. 4; Additional file 1 Fig. S5; Abstract | Gene-body methylation peaks in genes originating in the eukaryote ancestor (phylostratum 2), not in the oldest genes (phylostratum 1), and declines with younger gene age, with lineage-specific genes (phylostrata 8-10) lowest. | **not_testable** | Needs a phylostratigraphic assignment of genes; not part of our outputs. |  |
| C22 | P1 Results; Fig. 1; Additional file 1 Fig. S1; Abstract | The genome encodes one ortholog each of DNMT1, DNMT2 and DNMT3 plus MBD2/3 and UHRF1, expressed across development (DNMT1 and UHRF1 co-expressed, Pearson r = 0.86). | **not_testable** | Gene-content and expression claim, not methylation. |  |
| C23 | P1 Results; Additional file 1 Tables S1-S2 | The wild animal's reads mapped less well than the inbred animal's (unique-mapped data ~11.8x vs ~15.7x per strand; 73.9% vs 95.0% of genomic Cs covered >=2x), attributed to sequence divergence from the inbred-derived reference; lambda-based error rate was ~0.4% per library. | **not_testable** | Mapping/QC claim tied to the authors' inbred-derived reference; our reference is a different animal, so CpGs covered per library is only an indirect check. |  |

## Reading this

Differential results here are per gene (gene body or promoter, summed calls, DSS beta-binomial model, FDR < 0.05), on GCF_963853765.1. Papers count CpGs or regions on other assemblies, so counts are compared by presence, and named genes are matched by symbol or product description. See docs/replication.md.
