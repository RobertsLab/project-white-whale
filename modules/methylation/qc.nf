// Per-sample CpG summaries and the QC gate (workflows/methylation/bin/
// methylome_qc.py). qc_summary.tsv is published next to run_manifest.json;
// `whale record-run` commits both to the study branch when the run ends.

process CPG_SUMMARY {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(sample_id), path(cov)

    output:
    path "${sample_id}.cpg_summary.tsv", emit: summary

    script:
    """
    methylome_qc.py cpg-summary ${sample_id} ${cov} ${sample_id}.cpg_summary.tsv
    """

    stub:
    """
    printf 'sample_id\\tcpgs_1x\\tcpgs_5x\\tcpgs_10x\\tmean_depth\\tcpg_meth_pooled_pct\\tcpg_meth_mean_5x_pct\\n${sample_id}\\t0\\t0\\t0\\t0\\t\\t\\n' > ${sample_id}.cpg_summary.tsv
    """
}

process QC_GATE {
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
    methylome_qc.py qc-table --samplesheet ${samplesheet} --out qc_summary.tsv \\
        --min-mapping ${params.qc_min_mapping} --min-conversion ${params.qc_min_conversion} \\
        --min-cpgs-5x ${params.qc_min_cpgs_5x} reports/*
    """

    stub:
    """
    methylome_qc.py qc-table --samplesheet ${samplesheet} --out qc_summary.tsv reports/*
    """
}
