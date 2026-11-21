// Salmon: a decoy-aware transcriptome index, built once per assembly, and
// per-library quantification. Adapted from AREE's
// modules/rnaseq/align_quant_salmon.nf (VENDORED_FROM_AREE.md), which built
// a plain index and left decoys out; here the whole genome is the decoy.

// The index is built from the reference bundle's rna.fna with the bundle's
// genomic.fna appended as decoy sequence (Srivastava et al. 2020, Genome
// Biology 21:239): a read from an intron, an unannotated locus or genomic
// DNA that matches the genome better than any transcript is discarded,
// instead of being forced onto a similar transcript. Kept (storeDir) as
// <index dir>/salmon_index, by default inside the reference bundle like the
// Bismark index, so later runs on the machine reuse it. Its
// whale_index.json records the inputs' SHA-256, k and Salmon version, which
// main.nf checks against the bundle and params before a run; transcripts.txt
// lists the transcript names (the stub quantification reads them).
process SALMON_INDEX {
    tag "${transcripts.simpleName}"
    label 'process_high'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'
    storeDir "${params.salmon_index_dir ?: "${params.reference_dir}/${params.genome_assembly}"}"

    input:
    path transcripts
    path genome          // genomic.fna for the decoys, or [] for a transcript-only index
    val index_name

    output:
    path "${index_name}", emit: index

    script:
    def decoys = genome ? true : false
    def sha_genome = decoys ? "\"\$(sha256sum ${genome} | cut -d' ' -f1)\"" : 'null'
    def opts = "-k ${params.salmon_kmer} -p ${task.cpus} -i ${index_name} ${params.salmon_index_args}"
    def build = decoys
        ? "grep '^>' ${genome} | cut -d' ' -f1 | sed 's/^>//' > decoys.txt && " +
          "cat ${transcripts} ${genome} > gentrome.fa && salmon index -t gentrome.fa -d decoys.txt ${opts} && rm gentrome.fa"
        : "salmon index -t ${transcripts} ${opts}"
    """
    ${build}
    grep '^>' ${transcripts} | cut -d' ' -f1 | sed 's/^>//' > ${index_name}/transcripts.txt
    cat > ${index_name}/whale_index.json <<EOF_JSON
    {
      "salmon": "\$(salmon --version | sed 's/salmon //')",
      "kmer": ${params.salmon_kmer},
      "decoys": ${decoys},
      "index_args": "${params.salmon_index_args}",
      "transcripts": "${transcripts.name}",
      "transcripts_sha256": "\$(sha256sum ${transcripts} | cut -d' ' -f1)",
      "genome_sha256": ${sha_genome},
      "decoy_sequences": ${decoys ? "\$(wc -l < decoys.txt | tr -d ' ')" : 0}
    }
    EOF_JSON
    """

    stub:
    def decoys = genome ? true : false
    """
    mkdir ${index_name}
    grep '^>' ${transcripts} | cut -d' ' -f1 | sed 's/^>//' > ${index_name}/transcripts.txt
    cat > ${index_name}/whale_index.json <<EOF_JSON
    {"salmon": "stub", "kmer": ${params.salmon_kmer}, "decoys": ${decoys}, "index_args": "${params.salmon_index_args}",
     "transcripts": "${transcripts.name}", "transcripts_sha256": "\$(sha256sum ${transcripts} | cut -d' ' -f1)",
     "genome_sha256": ${decoys ? "\"\$(sha256sum ${genome} | cut -d' ' -f1)\"" : 'null'}, "decoy_sequences": 0}
    EOF_JSON
    """
}

// Library type is detected per library (-l A) and reported in qc_summary.tsv,
// so a study that mixes stranded and unstranded libraries shows it. Selective
// alignment (Salmon's default since 1.0) with fragment GC-bias correction for
// paired-end libraries only: single-end reads do not observe fragment length,
// which the model needs. --tag_seq (3' tag libraries such as QuantSeq, about
// one read per transcript whatever its length) switches length correction and
// GC bias off.
process SALMON_QUANT {
    tag "${sample_id}"
    label 'process_high'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'
    publishDir "${params.outdir}/${params.study_id}/rnaseq/salmon", mode: params.publish_mode, saveAs: { it == "versions.yml" ? null : it }

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant
    path "versions.yml", emit: versions

    script:
    def fq = reads instanceof Collection ? reads as List : [reads]
    def input = fq.size() == 1 ? "-r ${fq[0]}" : "-1 ${fq[0]} -2 ${fq[1]}" + (params.tag_seq ? '' : ' --gcBias')
    def length_opt = params.tag_seq ? '--noLengthCorrection' : ''
    """
    salmon quant -i ${index} -l A ${input} -p ${task.cpus} --validateMappings ${length_opt} \\
        ${params.salmon_quant_args} -o ${sample_id}
    printf 'SALMON_QUANT:\\n  salmon: "%s"\\n' "\$(salmon --version | sed 's/salmon //')" > versions.yml
    """

    stub:
    def paired = (reads instanceof Collection ? (reads as List).size() : 1) == 2
    """
    mkdir -p ${sample_id}/aux_info
    { printf 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n'
      while read t; do printf '%s\\t1000\\t800\\t0\\t0\\n' "\$t"; done < ${index}/transcripts.txt; } > ${sample_id}/quant.sf
    printf '{"num_processed": 0, "num_mapped": 0, "num_decoy_fragments": 0, "percent_mapped": 0.0, "library_types": ["${paired ? 'ISR' : 'SR'}"]}\\n' > ${sample_id}/aux_info/meta_info.json
    printf '{"expected_format": "${paired ? 'ISR' : 'SR'}", "compatible_fragment_ratio": 0.0}\\n' > ${sample_id}/lib_format_counts.json
    printf 'SALMON_QUANT:\\n  salmon: "stub"\\n' > versions.yml
    """
}
