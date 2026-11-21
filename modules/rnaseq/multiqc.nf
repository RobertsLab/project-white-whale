// One MultiQC report over FastQC, fastp and Salmon, for people. The numbers
// the recorder commits come from qc_summary.tsv (RNASEQ_QC), not from here.
process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'
    publishDir "${params.outdir}/${params.study_id}/rnaseq", mode: params.publish_mode, pattern: "*.html"

    input:
    path inputs, stageAs: "inputs/*"

    output:
    path "multiqc_report.html", emit: report
    path "versions.yml", emit: versions

    script:
    """
    multiqc --force --no-data-dir --filename multiqc_report.html --title "${params.study_id} RNA-seq" inputs
    printf 'MULTIQC:\\n  multiqc: "%s"\\n' "\$(multiqc --version | sed 's/.*version //')" > versions.yml
    """

    stub:
    """
    touch multiqc_report.html
    printf 'MULTIQC:\\n  multiqc: "stub"\\n' > versions.yml
    """
}
