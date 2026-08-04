# File Size Estimation Methodology and Calculations

## Estimation Approach

### Base Assumptions for File Size Calculations

#### RNA-seq Data
- **Average reads per sample**: 50-100 million paired-end reads
- **Read length**: 100-150 bp (modern studies)
- **Compression ratio**: ~0.3-0.5 (gzipped FASTQ vs raw)
- **File size per million reads**: ~100-150 MB (compressed)
- **Typical sample size**: 2-5 GB per sample

#### DNA Methylation Data

##### Whole Genome Bisulfite Sequencing (WGBS)
- **Coverage target**: 10-30X genome coverage
- **Oyster genome size**: ~650 Mb
- **Reads needed**: 200-600 million reads per sample
- **File size per sample**: 10-30 GB

##### Reduced Representation Bisulfite Sequencing (RRBS)
- **Target regions**: ~1-5% of genome
- **Reads needed**: 20-100 million reads per sample
- **File size per sample**: 2-6 GB

##### MeDIP-seq
- **Enrichment method**: Methylated DNA immunoprecipitation
- **Reads needed**: 30-100 million reads per sample
- **File size per sample**: 3-8 GB

## Detailed Size Calculations

### RNA-seq Dataset Size Estimates

| Study Type | Samples | Reads/Sample (M) | Size/Sample (GB) | Total Size (GB) |
|------------|---------|------------------|------------------|-----------------|
| Early developmental (2010-2015) | 50 | 30-50 | 1.5-2.5 | 75-125 |
| Stress response (2015-2020) | 100 | 50-100 | 2.5-5.0 | 250-500 |
| Population studies (2015-2020) | 80 | 40-80 | 2.0-4.0 | 160-320 |
| Recent studies (2020-2024) | 120 | 80-150 | 4.0-7.5 | 480-900 |
| Magallana studies (2018-2024) | 50 | 60-120 | 3.0-6.0 | 150-300 |
| **SUBTOTAL** | **400** | **Variable** | **Variable** | **1,115-2,145** |

### DNA Methylation Dataset Size Estimates (verified actuals, NCBI SRA 2026-08)

These are the real run counts and summed SRA download sizes for the verified
BioProjects (see
[`../ncbi-datasets/dna-methylation-datasets.md`](../ncbi-datasets/dna-methylation-datasets.md)),
not modelled estimates. The per-sample assumptions above still describe the
methods, but the totals below supersede them for these specific projects.

| Dataset (BioProject) | Method | Runs | Total Size (GB) |
|----------------------|--------|------|-----------------|
| POMS adaptation (PRJEB60400) | WGBS | 246 | 450-500 |
| GESTINOV POMS (PRJEB81880) | EM-seq | 40 | 400-450 |
| Aging / DECICOMP (PRJEB105019) | WGBS | 60 | 450-500 |
| PESTO (PRJEB58545) | WGBS | 48 | 350-400 |
| Transgen. infection (PRJNA609264) | WGBS | 47 | 600-650 |
| pH / ploidy (PRJNA682817) | WGBS | 24 | 100-110 |
| Developmental (PRJNA324546) | MeDIP-seq | 21 | 5-10 |
| Epigenomics series (PRJNA807732/562805/213124) | WGBS | 50 | 300-330 |
| **SUBTOTAL** | **WGBS/EM-seq/MeDIP** | **~536** | **~2,655-2,950** |

## Size Distribution by Research Theme

### Environmental Stress Studies
| Data Type | Samples | Size Range (GB) | Best Estimate (GB) |
|-----------|---------|-----------------|-------------------|
| Temperature stress RNA-seq | 80 | 160-400 | 280 |
| Temperature stress WGBS | 25 | 200-500 | 350 |
| pH stress RNA-seq | 60 | 120-300 | 210 |
| pH stress methylation | 20 | 160-400 | 280 |
| Salinity stress RNA-seq | 40 | 80-200 | 140 |
| Multiple stress RRBS | 30 | 60-180 | 120 |
| **Stress Studies Total** | **255** | **780-1,980** | **1,380** |

### Developmental Studies
| Data Type | Samples | Size Range (GB) | Best Estimate (GB) |
|-----------|---------|-----------------|-------------------|
| Embryonic RNA-seq | 30 | 45-150 | 100 |
| Larval RNA-seq | 50 | 75-250 | 160 |
| Adult development RNA-seq | 40 | 80-200 | 140 |
| Developmental WGBS | 20 | 200-500 | 350 |
| Reproductive RRBS | 25 | 50-150 | 100 |
| **Development Studies Total** | **165** | **450-1,250** | **850** |

## Storage and Transfer Considerations

### Data Storage Requirements
- **Verified methylation raw data**: ~2.6-3.0 TB (the 10 BioProjects above)
- **Processed data**: Additional 20-30% of raw data size
- **Metadata and documentation**: ~1-2% of total size
- **RNA-seq**: not included (accessions pending verification)

### Data Transfer Estimates
- **Typical internet speed**: 100 Mbps = 12.5 MB/s = 45 GB/hour
- **Time to download 1 TB**: ~22-25 hours
- **Full dataset download**: 75-200 hours (3-8 days continuous)

### Recommended Download Strategy
1. **Prioritize by research question**: Focus on specific subsets first
2. **Start with processed data**: Often 10-20% the size of raw data
3. **Use institutional networks**: Higher bandwidth at universities
4. **Parallel downloads**: Multiple simultaneous downloads where possible

## Quality and Coverage Metrics

### RNA-seq Quality Indicators
- **Read depth**: 20M+ reads for basic analysis, 50M+ for comprehensive
- **Mapping rate**: >80% to reference genome
- **rRNA contamination**: <10% for good libraries
- **Insert size**: 200-800 bp for paired-end libraries

### Methylation Quality Indicators
- **Bisulfite conversion**: >95% for reliable results
- **Coverage depth**: 5X minimum, 10X+ preferred for WGBS
- **CpG coverage**: 60%+ of genome CpGs for WGBS
- **Duplication rate**: <30% for good libraries

## Cost-Benefit Analysis for Data Access

### High-Value Datasets (Priority 1)
- POMS adaptation WGBS (PRJEB60400): largest verified WGBS set (246 runs)
- Transgenerational infection WGBS (PRJNA609264): disease-response epigenetics
- Developmental methylome MeDIP-seq (PRJNA324546): small, quick to download

### Medium-Value Datasets (Priority 2)
- Population-specific studies: Geographic relevance
- Developmental time series: Temporal dynamics
- Recent Magallana studies: Updated nomenclature

### Lower-Priority Datasets
- Single-tissue studies: Limited scope
- Low-coverage studies: Technical limitations
- Unpublished/proprietary: Access restrictions

## File Format Breakdown

### Expected File Types and Sizes
| File Type | Description | Typical Size | Purpose |
|-----------|-------------|--------------|---------|
| .fastq.gz | Raw sequencing reads | 2-20 GB | Primary data |
| .bam | Aligned reads | 1-10 GB | Processed alignments |
| .bw/.bigwig | Coverage tracks | 50-500 MB | Visualization |
| .bed | Genomic intervals | 1-100 MB | Annotations |
| .txt/.tsv | Count matrices | 10-100 MB | Analysis ready |

## Last Updated
August 2026 — methylation dataset sizes replaced with verified NCBI SRA actuals.