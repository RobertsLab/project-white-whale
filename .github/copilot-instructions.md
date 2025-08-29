# White Whale - Crassostrea gigas Public Datasets Repository

This repository contains research documentation for identifying and cataloging public datasets with RNA-seq and DNA methylation data for *Crassostrea gigas* (Pacific oyster) and *Magallana gigas* (updated genus name). This is a documentation-only repository focused on research methodology and dataset cataloging.

Always reference these instructions first and fallback to search or exploration only when you encounter information that does not match what is documented here.

## Working Effectively

### Repository Setup and Validation
- **No build system required**: This is a documentation repository with no code to compile or applications to run
- **No dependencies to install**: Repository contains only Markdown documentation files
- **No tests to execute**: Validation focuses on content accuracy and documentation quality
- **Version control**: Standard git workflow for documentation changes

### Content Validation Process
- **Markdown syntax**: Always validate Markdown formatting using standard linters
- **Link verification**: Check all internal links between documentation files work correctly
- **Data accuracy**: Cross-reference dataset information with source databases when possible
- **Citation consistency**: Ensure all scientific references follow consistent formatting
- **File organization**: Maintain the established directory structure and naming conventions

### Documentation Standards
- **File naming**: Use lowercase with hyphens (kebab-case) for all .md files
- **Headers**: Use consistent header hierarchy (# for main title, ## for sections, ### for subsections)
- **Tables**: Format data tables consistently with proper alignment
- **Scientific names**: Always italicize species names (*Crassostrea gigas*, *Magallana gigas*)
- **Numbers**: Use consistent formatting for file sizes (GB, TB) and sample counts
- **Dates**: Use format "Month YYYY" (e.g., "December 2024") for last updated fields

## Repository Structure and Navigation

### Directory Organization
```
/
├── README.md                    # Main repository overview
├── action-plan.md              # Research methodology and search strategies
├── ncbi-datasets/              # NCBI database search results
│   ├── README.md
│   ├── search-results.md       # Database search strategies and outcomes
│   ├── rna-seq-datasets.md     # RNA-seq dataset documentation
│   └── dna-methylation-datasets.md # DNA methylation dataset documentation
├── literature-review/          # Published literature analysis
│   ├── README.md
│   └── key-publications.md     # Major publications with associated datasets
├── dataset-summaries/          # Consolidated analysis
│   ├── README.md
│   └── consolidated-summary.md # Comprehensive dataset overview
└── file-size-estimates/        # Storage requirement calculations
    ├── README.md
    └── calculation-methodology.md # Size estimation methodology
```

### Key Files to Reference
- **action-plan.md**: Complete research methodology including search terms and database strategies
- **consolidated-summary.md**: Overall statistics and findings summary
- **calculation-methodology.md**: File size estimation formulas and assumptions

## Research Methodology

### Dataset Identification Process
1. **NCBI Database Mining**: Systematic searches of SRA, GEO, and BioProject databases
2. **Literature Review**: PubMed and journal-specific searches for relevant publications
3. **International Repository Search**: ENA, DDBJ, and institutional repositories
4. **Data Cataloging**: Comprehensive information collection for each dataset
5. **Size Estimation**: Calculate storage requirements using validated methodologies
6. **Access Strategy**: Prioritization and download planning

### Search Terms and Strategies
**Primary species searches**:
- "Crassostrea gigas"[Organism] AND "RNA-Seq"[Strategy]
- "Crassostrea gigas"[Organism] AND "Bisulfite-Seq"[Strategy]
- "Magallana gigas"[Organism] AND "RNA-Seq"[Strategy]
- "Magallana gigas"[Organism] AND "Bisulfite-Seq"[Strategy]

**Additional keyword searches**:
- ("Crassostrea gigas" OR "Magallana gigas") AND ("RNA-seq" OR "transcriptome")
- ("Crassostrea gigas" OR "Magallana gigas") AND ("methylation" OR "bisulfite")

### Data Collection Requirements
For each identified dataset, document:
- **Accession numbers**: SRA, GEO, BioProject IDs
- **Sample count**: Number of biological replicates
- **Tissue types**: Specific organs/tissues analyzed
- **Environmental conditions**: Experimental treatments
- **Technical details**: Sequencing platform, read length, depth
- **File sizes**: Actual or estimated data volume
- **Publication status**: Associated papers and citations

## Content Maintenance

### Adding New Datasets
1. **Verify dataset exists**: Confirm accession numbers and accessibility
2. **Collect required information**: Follow the data collection requirements above
3. **Update relevant files**: Add to appropriate category files and update summaries
4. **Recalculate totals**: Update size estimates and sample counts in consolidated-summary.md
5. **Update timestamps**: Modify "Last Updated" fields in affected files
6. **Cross-reference**: Ensure consistency across all documentation files

### File Size Estimation Guidelines
**RNA-seq data assumptions**:
- Typical sample size: 2-5 GB per sample (paired-end, 50-100M reads)
- High-depth samples: 5-10 GB per sample
- Compressed FASTQ: ~50% reduction from raw

**DNA methylation data assumptions**:
- WGBS samples: 10-30 GB per sample (high coverage)
- RRBS samples: 2-6 GB per sample (reduced representation)
- MeDIP-seq: 3-8 GB per sample

### Quality Assurance
- **Cross-reference sources**: Verify information against multiple databases
- **Check for duplicates**: Ensure no datasets are counted multiple times
- **Validate calculations**: Double-check all size estimates and totals
- **Update regularly**: Research is ongoing; documentation should reflect current status

## Common Tasks

### Updating Dataset Information
When adding new datasets or modifying existing entries:
1. **Start with specific category files** (rna-seq-datasets.md, dna-methylation-datasets.md)
2. **Update consolidated-summary.md** with new totals
3. **Modify calculation-methodology.md** if estimation methods change
4. **Update README files** in affected directories
5. **Check action-plan.md** for any methodology updates needed

### Validating Documentation
- **Markdown linting**: Use markdownlint or similar tools to check syntax
- **Link checking**: Verify all internal references work correctly
- **Table formatting**: Ensure all tables have proper alignment and headers
- **Consistent terminology**: Use established terms throughout documentation
- **Scientific accuracy**: Cross-check data against original sources when possible

### Research Scope
This research focuses on publicly available datasets containing:
- **RNA-seq data**: Transcriptome studies across various conditions
- **DNA methylation data**: Including WGBS, RRBS, MeDIP-seq, and targeted approaches

**Research themes covered**:
- Environmental stress response (temperature, pH, salinity, hypoxia)
- Developmental biology (embryogenesis, larval development, reproduction)
- Epigenetic mechanisms (methylation patterns, transgenerational effects)
- Population and comparative studies (geographic variation, species comparisons)

## Validation Scenarios

### Content Accuracy Validation
After making changes to documentation:
1. **Verify all internal links work**: Check references between files
2. **Confirm data consistency**: Ensure numbers match across related files
3. **Check scientific accuracy**: Validate species names, technical terms, and citations
4. **Review formatting**: Ensure consistent Markdown style and table formatting
5. **Update summaries**: Reflect changes in overview documents

### Documentation Completeness
- **All datasets have required fields**: Accession numbers, sample counts, file sizes
- **Cross-references are current**: Links between related documents work
- **Calculations are accurate**: Size estimates use correct methodologies
- **Timestamps are current**: "Last Updated" fields reflect recent changes

## Repository Statistics
- **Total documentation files**: ~1,040 lines across 13 Markdown files
- **Estimated datasets covered**: 300-830 samples across both data types
- **Estimated total data size**: 1.3-3.6 TB
- **Time span covered**: 2010-2024
- **Geographic scope**: Global studies

Always maintain these statistics when updating documentation and recalculate totals when adding new datasets.