// Adapter and quality trimming with fastp. Adapted from AREE's
// modules/rnaseq/trim_fastp.nf (VENDORED_FROM_AREE.md): same container pin,
// explicit single-end adapter, --detect_adapter_for_pe for pairs, poly-X
// trimming for 3' tag libraries. Trimmed reads are not published: they are
// tens of GB per study and regenerable from the raw reads.
process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
    publishDir "${params.outdir}/${params.study_id}/rnaseq/qc", mode: params.publish_mode, pattern: "*.fastp.{json,html}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R*.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    path "${sample_id}.fastp.html", emit: html
    path "versions.yml", emit: versions

    script:
    // A single file stages as a Path, not a list. Single-end adapters are
    // given explicitly (params.se_adapter), not auto-detected: AREE found
    // fastp's single-end detection picking an abundant oyster transcript as
    // the "adapter" (DUPERRET2025_STARVE_RNASEQ), a different guess per
    // library. Paired-end detection uses read overlap and is reliable.
    def fq = reads instanceof Collection ? reads as List : [reads]
    def se_adapter = params.se_adapter == 'auto' ? '' : "--adapter_sequence ${params.se_adapter}"
    def io = fq.size() == 1 \
        ? "-i ${fq[0]} -o ${sample_id}_trimmed_R1.fastq.gz ${se_adapter}" \
        : "-i ${fq[0]} -I ${fq[1]} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz --detect_adapter_for_pe"
    def tag_seq = params.tag_seq ? '--trim_poly_x' : ''
    """
    fastp ${io} ${tag_seq} ${params.fastp_args} \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --report_title "${sample_id} fastp report"
    printf 'FASTP:\\n  fastp: "%s"\\n' "\$(fastp --version 2>&1 | sed -e 's/fastp //g')" > versions.yml
    """

    stub:
    def n = reads instanceof Collection ? (reads as List).size() : 1
    """
    touch ${sample_id}_trimmed_R1.fastq.gz ${n == 2 ? "${sample_id}_trimmed_R2.fastq.gz" : ''} ${sample_id}.fastp.html
    printf '{"summary": {"before_filtering": {"total_reads": 0}, "after_filtering": {"total_reads": 0}}${n == 2 ? ', "read2_before_filtering": {}' : ''}}\\n' > ${sample_id}.fastp.json
    printf 'FASTP:\\n  fastp: "stub"\\n' > versions.yml
    """
}
