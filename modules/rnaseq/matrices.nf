// Study-level matrices and the QC table (workflows/rnaseq/bin/rnaseq_quant.py).
// Published to <outdir>/<STUDY_ID>/rnaseq/; `whale record-rnaseq` commits the
// gene-level matrices and qc_summary.tsv and uploads the transcript-level ones.

process GENE_MATRICES {
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.outdir}/${params.study_id}/rnaseq", mode: 'copy'

    input:
    path quant_dirs, stageAs: "salmon/*"
    path samplesheet
    path tx2gene

    output:
    path "matrices", emit: dir

    script:
    """
    rnaseq_quant.py matrices --samplesheet ${samplesheet} --tx2gene ${tx2gene} --out-dir matrices salmon/*
    """

    stub:
    """
    rnaseq_quant.py matrices --samplesheet ${samplesheet} --tx2gene ${tx2gene} --out-dir matrices salmon/*
    """
}

process RNASEQ_QC {
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.outdir}/${params.study_id}/rnaseq", mode: 'copy'

    input:
    path reports, stageAs: "reports/*"
    path matrices
    path samplesheet

    output:
    path "qc_summary.tsv", emit: summary

    script:
    """
    rnaseq_quant.py qc-table --samplesheet ${samplesheet} --matrices ${matrices} --out qc_summary.tsv \\
        --min-mapping ${params.rnaseq_min_mapping} --min-mapped ${params.rnaseq_min_mapped} \\
        --min-genes ${params.rnaseq_min_genes} reports/*
    """

    stub:
    """
    rnaseq_quant.py qc-table --samplesheet ${samplesheet} --matrices ${matrices} --out qc_summary.tsv reports/*
    """
}
