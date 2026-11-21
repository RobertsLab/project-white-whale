// Per-sample enrichment in 1 kb tiles, gene bodies and promoters, and the
// study's QC table (workflows/enrichment/bin/enrichment_counts.py; method in
// docs/enrichment.md). `whale record-enrichment` uploads the per-sample tables
// and commits gene x sample matrices built from them.

process ENRICHMENT_COUNT {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/tiles", mode: params.publish_mode, pattern: "*.tiles_1kb.tsv.gz"
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/features", mode: params.publish_mode, pattern: "*.enrichment_features.tsv.gz"
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/qc", mode: params.publish_mode, pattern: "*.enrichment_qc.tsv"

    input:
    tuple val(sample_id), path(fragments)
    path genome_fasta
    path features

    output:
    path "${sample_id}.tiles_1kb.tsv.gz", emit: tiles
    path "${sample_id}.enrichment_features.tsv.gz", emit: features
    path "${sample_id}.enrichment_qc.tsv", emit: qc

    script:
    """
    enrichment_counts.py count ${sample_id} ${fragments} ${genome_fasta} ${features} . \\
        --window ${params.window_size} --fragment-length ${params.fragment_length} \\
        --min-cpgs ${params.enrichment_min_cpgs} --seed ${params.saturation_seed}
    """

    stub:
    """
    enrichment_counts.py count ${sample_id} ${fragments} ${genome_fasta} ${features} . \\
        --window ${params.window_size} --fragment-length ${params.fragment_length} \\
        --min-cpgs ${params.enrichment_min_cpgs} --seed ${params.saturation_seed}
    """
}

process ENRICHMENT_QC_GATE {
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}", mode: 'copy'

    input:
    path reports, stageAs: "reports/*"
    path samplesheet

    output:
    path "qc_summary.tsv", emit: summary

    script:
    """
    enrichment_counts.py qc-table --samplesheet ${samplesheet} --out qc_summary.tsv \\
        --min-mapping ${params.enrichment_min_mapping} --min-fragments ${params.enrichment_min_fragments} \\
        --min-relh ${params.enrichment_min_relh} reports/*
    """

    stub:
    """
    enrichment_counts.py qc-table --samplesheet ${samplesheet} --out qc_summary.tsv reports/*
    """
}
