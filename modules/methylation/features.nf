// Per-sample methylation of every gene's body, promoter, exons and introns
// (workflows/methylation/bin/methylome_features.py), against the reference
// bundle's features/. Published per study under features/; `whale
// record-run` uploads each table with the sample's methylome and commits
// gene x sample matrices built from them.
process FEATURE_SUMMARY {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/features", mode: params.publish_mode

    input:
    tuple val(sample_id), path(cov)
    path features
    path mask       // a SNP mask's cpg_mask.bed.gz, or [] for none

    output:
    path "${sample_id}.features.tsv.gz", emit: table

    script:
    """
    methylome_features.py ${sample_id} ${cov} ${features} ${sample_id}.features.tsv.gz ${mask ?: ''}
    """

    stub:
    """
    methylome_features.py ${sample_id} ${cov} ${features} ${sample_id}.features.tsv.gz ${mask ?: ''}
    """
}
