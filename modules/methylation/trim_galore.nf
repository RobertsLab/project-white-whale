// Adapter and quality trimming. `extra_args` carries the chemistry preset
// (e.g. --rrbs) plus params.trim_args; M-bias-driven clipping belongs in
// params.extractor_args, not here, so trimmed reads stay comparable.
process TRIM_GALORE {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
    publishDir "${params.outdir}/${params.study_id}/qc", mode: params.publish_mode, pattern: "*report.txt"

    input:
    tuple val(sample_id), path(reads)
    val extra_args

    output:
    tuple val(sample_id), path("*_val_*.fq.gz"), emit: reads
    path "*report.txt", emit: reports
    path "versions.yml", emit: versions

    script:
    def paired = reads instanceof List && reads.size() == 2
    def cores = Math.max(1, (task.cpus as int) / 4 as int)
    if (paired) {
        """
        trim_galore --paired --cores ${cores} --gzip --basename ${sample_id} ${extra_args} ${reads[0]} ${reads[1]}
        printf 'TRIM_GALORE:\\n  trim_galore: "%s"\\n' "\$(trim_galore --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)" > versions.yml
        """
    } else {
        """
        trim_galore --cores ${cores} --gzip --basename ${sample_id} ${extra_args} ${reads}
        mv ${sample_id}_trimmed.fq.gz ${sample_id}_val_1.fq.gz
        printf 'TRIM_GALORE:\\n  trim_galore: "%s"\\n' "\$(trim_galore --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)" > versions.yml
        """
    }

    stub:
    def paired = reads instanceof List && reads.size() == 2
    """
    touch ${sample_id}_val_1.fq.gz ${paired ? "${sample_id}_val_2.fq.gz" : ''} ${sample_id}_trimming_report.txt
    printf 'TRIM_GALORE:\\n  trim_galore: "stub"\\n' > versions.yml
    """
}
