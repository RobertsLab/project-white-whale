// Bowtie2 index and alignment for enrichment-tier libraries (MeDIP, MBD-seq):
// plain, unconverted reads against the reference bundle's genome. One pinned
// image carries bowtie2 2.5.2 and samtools 1.18 (the nf-core bowtie2 image),
// so alignment, sorting and mate fixing stream without an intermediate SAM.

// Built once per assembly and kept (storeDir) as <index dir>/bowtie2_genome,
// by default inside the reference bundle, beside bismark_genome. Its
// genome.sha256 records the FASTA it was built from, which main.nf checks.
process BOWTIE2_BUILD {
    tag "${fasta.simpleName}"
    label 'process_high'
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
    storeDir "${params.bowtie2_index_dir ?: "${params.reference_dir}/${params.genome_assembly}"}"

    input:
    path fasta

    output:
    path "bowtie2_genome", emit: index

    script:
    """
    mkdir bowtie2_genome
    sha256sum ${fasta} | cut -d' ' -f1 > bowtie2_genome/genome.sha256
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_genome/genome
    """

    stub:
    """
    mkdir bowtie2_genome
    sha256sum ${fasta} | cut -d' ' -f1 > bowtie2_genome/genome.sha256
    touch bowtie2_genome/genome.1.bt2 bowtie2_genome/STUB   # main.nf refuses this index in a real run
    """
}

// End-to-end alignment (bowtie2's default --sensitive) plus params.bowtie2_args;
// name-collated, mate-fixed (for markdup) and coordinate-sorted in one stream.
process BOWTIE2_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
    publishDir "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}/qc", mode: params.publish_mode, pattern: "*.bowtie2.log"

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    path "${sample_id}.bowtie2.log", emit: log
    path "versions.yml", emit: versions

    script:
    def paired = reads instanceof List && reads.size() == 2
    def input = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"
    def threads = Math.max(1, (task.cpus as int) - 2)
    """
    bowtie2 -p ${threads} -x ${index}/genome ${input} ${params.bowtie2_args} 2> ${sample_id}.bowtie2.log \\
        | samtools collate -O -u - tmp_collate \\
        | samtools fixmate -m -u - - \\
        | samtools sort -@ 2 -m 1G -o ${sample_id}.sorted.bam -
    cat ${sample_id}.bowtie2.log
    printf 'BOWTIE2_ALIGN:\\n  bowtie2: "%s"\\n  samtools: "%s"\\n' \\
        "\$(bowtie2 --version 2>&1 | head -1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)" \\
        "\$(samtools --version | head -1 | cut -d' ' -f2)" > versions.yml
    """

    stub:
    """
    touch ${sample_id}.sorted.bam
    printf '0 reads; of these:\\n0.00%% overall alignment rate\\n' > ${sample_id}.bowtie2.log
    printf 'BOWTIE2_ALIGN:\\n  bowtie2: "stub"\\n  samtools: "stub"\\n' > versions.yml
    """
}
