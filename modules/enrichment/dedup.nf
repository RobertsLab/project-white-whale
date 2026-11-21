// Duplicate removal and the fragment table the counting step reads.
//
// samtools markdup -r removes PCR/optical duplicates (both mates of a pair,
// or single reads at the same 5' position and strand), as the assay_chemistry
// vocabulary asks for MeDIP and MBD (deduplicate: true). Its stats file gives
// the duplication rate in qc_summary.tsv.
//
// The fragment table keeps only what counting needs: FLAG, RNAME, POS, CIGAR,
// TLEN of primary, mapped, non-duplicate records with MAPQ >= params.min_mapq
// (and, for paired libraries, proper pairs only).
process MARKDUP {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/qc", mode: params.publish_mode, pattern: "*.markdup.txt"
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/bam", mode: params.publish_mode, pattern: "*.dedup.bam*", enabled: params.publish_bams

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fragments.tsv.gz"), emit: fragments
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    path "${sample_id}.markdup.txt", emit: stats

    script:
    """
    samtools markdup -r -@ ${task.cpus} -f ${sample_id}.markdup.txt ${bam} ${sample_id}.dedup.bam
    samtools index ${sample_id}.dedup.bam
    # Paired library if the first record has FLAG bit 1; then proper pairs only.
    FLAG=\$(samtools view ${sample_id}.dedup.bam | head -1 | cut -f2)
    FILTER="-F 0x904"
    if [ -n "\$FLAG" ] && [ \$(( FLAG & 1 )) -eq 1 ]; then FILTER="-f 0x2 -F 0x904"; fi
    samtools view \$FILTER -q ${params.min_mapq} ${sample_id}.dedup.bam \\
        | cut -f2,3,4,6,9 | gzip -c > ${sample_id}.fragments.tsv.gz
    """

    stub:
    """
    echo -n | gzip -c > ${sample_id}.fragments.tsv.gz
    touch ${sample_id}.dedup.bam ${sample_id}.dedup.bam.bai
    printf 'EXAMINED: 0\\nDUPLICATE TOTAL: 0\\n' > ${sample_id}.markdup.txt
    """
}
