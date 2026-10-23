// Bismark genome preparation, alignment, deduplication and methylation
// extraction. Container pins as used by AREE's methylation workflow.

process BISMARK_GENOME_PREPARATION {
    tag "${fasta.simpleName}"
    label 'process_high'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'

    input:
    path fasta

    output:
    path "bismark_genome", emit: genome_dir

    script:
    // Bismark only recognises .fa/.fasta(.gz) files, and NCBI's genome is
    // genomic.fna, so it goes in under a name Bismark accepts. It is copied,
    // not linked: later tasks see this folder inside their own container,
    // where the staged input's real path (the reference drive) is not mounted.
    """
    mkdir bismark_genome
    cp -L ${fasta} bismark_genome/genome.fa
    bismark_genome_preparation --parallel ${Math.max(1, (task.cpus as int) / 2 as int)} bismark_genome
    """

    stub:
    """
    mkdir -p bismark_genome/Bisulfite_Genome
    cp -L ${fasta} bismark_genome/genome.fa
    """
}

process BISMARK_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'
    publishDir "${params.outdir}/${params.study_id}/qc", mode: params.publish_mode, pattern: "*report.txt"

    input:
    tuple val(sample_id), path(reads)
    path genome_dir

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "*report.txt", emit: report
    path "versions.yml", emit: versions

    script:
    def paired = reads instanceof List && reads.size() == 2
    def input = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "${reads}"
    // Each Bismark instance runs ~3 bowtie2 threads; --multicore = cpus / 3.
    def multicore = Math.max(1, (task.cpus as int) / 3 as int)
    // No --basename: Bismark 0.24 refuses it together with --multicore, so
    // outputs are named after the first input file and renamed here.
    """
    bismark --genome ${genome_dir} --multicore ${multicore} --temp_dir . ${params.bismark_align_args} \\
        ${input}
    mv *_bismark_bt2*.bam ${sample_id}.bam
    for r in *_bismark_bt2*_report.txt; do mv "\$r" "${sample_id}_\${r#*_bismark_bt2_}"; done
    printf 'BISMARK_ALIGN:\\n  bismark: "%s"\\n' "\$(bismark --version | grep -Eo 'v[0-9.]+' | head -1)" > versions.yml
    """

    stub:
    """
    touch ${sample_id}.bam ${sample_id}_PE_report.txt
    printf 'BISMARK_ALIGN:\\n  bismark: "stub"\\n' > versions.yml
    """
}

process BISMARK_DEDUPLICATE {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'
    publishDir "${params.outdir}/${params.study_id}/qc", mode: params.publish_mode, pattern: "*report.txt"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: bam
    path "*report.txt", emit: report

    script:
    """
    deduplicate_bismark --bam --output_dir . --outfile ${sample_id}.dedup ${bam}
    mv ${sample_id}.dedup*.bam ${sample_id}.dedup.bam
    """

    stub:
    """
    touch ${sample_id}.dedup.bam ${sample_id}.deduplication_report.txt
    """
}

process BISMARK_METHYLATION_EXTRACTOR {
    tag "${sample_id}"
    label 'process_high'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'
    publishDir "${params.outdir}/${params.study_id}/qc", mode: params.publish_mode, pattern: "*.txt"
    publishDir "${params.outdir}/${params.study_id}/cpg", mode: params.publish_mode, pattern: "*.cov.gz"

    input:
    tuple val(sample_id), path(bam)
    path genome_dir

    output:
    tuple val(sample_id), path("*.CpG_report.merged_CpG_evidence.cov.gz"), emit: cpg
    tuple val(sample_id), path("*.M-bias.txt"), path("*splitting_report.txt"), emit: qc
    path "versions.yml", emit: versions

    script:
    // Strand-merged CpG calls are the per-sample unit of the shared matrices
    // (PLAN.md §6.2). Non-CpG context files are kept for the CHH-based
    // conversion estimate when a study has no spike-in.
    """
    # Paired-end if the first read's SAM flag has bit 1 set. Plain shell
    # arithmetic: awk's and() is GNU-only and the container's awk may not be.
    PE=""
    FLAG=\$(samtools view ${bam} 2>/dev/null | head -1 | cut -f2)
    if [ -n "\$FLAG" ] && [ \$(( FLAG & 1 )) -eq 1 ]; then PE="--paired-end --no_overlap"; fi
    bismark_methylation_extractor \$PE --gzip --bedGraph --comprehensive --merge_non_CpG \\
        --genome_folder ${genome_dir} --multicore ${Math.max(1, (task.cpus as int) / 3 as int)} \\
        ${params.extractor_args} ${bam}
    coverage2cytosine --merge_CpG --zero_based --gzip --genome_folder ${genome_dir} \\
        -o ${sample_id} *.bismark.cov.gz
    printf 'BISMARK_METHYLATION_EXTRACTOR:\\n  bismark: "%s"\\n' "\$(bismark --version | grep -Eo 'v[0-9.]+' | head -1)" > versions.yml
    """

    stub:
    """
    echo "" | gzip > ${sample_id}.CpG_report.merged_CpG_evidence.cov.gz
    touch ${sample_id}.M-bias.txt ${sample_id}_splitting_report.txt
    printf 'BISMARK_METHYLATION_EXTRACTOR:\\n  bismark: "stub"\\n' > versions.yml
    """
}
