process FASTQC {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/qc", mode: params.publish_mode, pattern: "*_fastqc.*"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), emit: reports
    path "versions.yml", emit: versions

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}
    printf 'FASTQC:\\n  fastqc: "%s"\\n' "\$(fastqc --version | sed 's/FastQC v//')" > versions.yml
    """

    stub:
    """
    for r in ${reads}; do b=\$(basename \$r .fastq.gz); touch \${b}_fastqc.zip \${b}_fastqc.html; done
    printf 'FASTQC:\\n  fastqc: "stub"\\n' > versions.yml
    """
}
