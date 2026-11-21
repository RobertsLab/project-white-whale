#!/usr/bin/env nextflow
/*
 * White Whale :: RNA-seq workflow — FASTQ in, gene and transcript abundance out,
 * for the studies whose methylation has paired RNA-seq (PLAN.md §7.2, Goal 2).
 *
 *   FASTQC
 *   FASTP -> SALMON_QUANT (decoy-aware index: rna.fna + genomic.fna as decoys,
 *            built once per assembly by SALMON_INDEX and cached like the Bismark index)
 *         -> GENE_MATRICES (tximport-style: transcript -> gene with tx2gene.tsv)
 *         -> RNASEQ_QC (qc_summary.tsv, the gate) -> MULTIQC -> MANIFEST
 *
 * Inputs come from the reference bundle (`whale fetch-reference --rnaseq`
 * adds rna.fna and tx2gene.tsv; docs/reference_bundle.md) and the study's
 * RNA-seq sample sheet, data/studies/<STUDY_ID>/rnaseq/samplesheet.tsv
 * (`whale fetch-rnaseq-samplesheet`), with columns sample_id, run_accession,
 * fastq_1, fastq_2 (empty for single-end). FASTQ paths are relative to
 * ${params.rnaseq_raw_dir} (default ${params.raw_dir}/<STUDY_ID>/rnaseq/) unless absolute.
 *
 * Everything is published under ${params.outdir}/<STUDY_ID>/rnaseq/:
 *   salmon/<sample>/          quant.sf, aux_info/, lib_format_counts.json, logs/
 *   matrices/                 gene_{tpm,counts,length}.tsv.gz, transcript_{tpm,counts}.tsv.gz
 *   qc/                       FastQC and fastp reports
 *   qc_summary.tsv, multiqc_report.html, run_manifest.json
 *
 * When a run succeeds it records itself with `whale record-rnaseq` (gene
 * matrices, QC table and manifest committed under data/studies/<STUDY_ID>/rnaseq/,
 * per-sample files uploaded to the data host); --record_status false turns that off.
 *
 *   nextflow run workflows/rnaseq -c config/<STUDY_ID>.config -profile raven
 *   nextflow run workflows/rnaseq -profile test -stub         # smoke test
 *   nextflow run workflows/rnaseq -profile test,docker        # every tool, simulated reads
 */

nextflow.enable.dsl = 2

include { FASTQC }        from '../../modules/rnaseq/fastqc.nf'
include { FASTP }         from '../../modules/rnaseq/fastp.nf'
include { SALMON_INDEX }  from '../../modules/rnaseq/salmon.nf'
include { SALMON_QUANT }  from '../../modules/rnaseq/salmon.nf'
include { GENE_MATRICES } from '../../modules/rnaseq/matrices.nf'
include { RNASEQ_QC }     from '../../modules/rnaseq/matrices.nf'
include { MULTIQC }       from '../../modules/rnaseq/multiqc.nf'
include { MANIFEST }      from '../../modules/rnaseq/manifest.nf'

def bundleDir() {
    return "${params.reference_dir}/${params.genome_assembly}"
}

// From the command line (--salmon_decoys false) the value arrives as a string.
def decoysOn() {
    return !(params.salmon_decoys.toString().toLowerCase() in ['false', '0', 'no'])
}

def indexName() {
    return decoysOn() ? 'salmon_index' : 'salmon_index_nodecoy'
}

def rawDir() {
    return params.rnaseq_raw_dir ?: "${params.raw_dir}/${params.study_id}/rnaseq"
}

def resolveFastq(String p) {
    if (!p) return null
    def f = p.startsWith('/') ? file(p) : file("${rawDir()}/${p}")
    if (!f.exists()) error "FASTQ not found: ${f}"
    return f
}

def bundleFile(String name) {
    def f = file("${bundleDir()}/${name}")
    if (!f.exists()) {
        error "${f} is missing. RNA-seq needs the reference bundle's transcript files: " +
              "whale fetch-reference --rnaseq (docs/reference_bundle.md)."
    }
    return f
}

// A cached Salmon index must have been built from the bundle's current
// rna.fna (and genomic.fna, for the decoys) with this run's settings.
def checkCachedIndex() {
    def index = file("${params.salmon_index_dir ?: bundleDir()}/${indexName()}")
    def info  = index.resolve('whale_index.json')
    if (!info.exists()) return
    def built = new groovy.json.JsonSlurper().parse(info.toFile())
    def stale = []
    if (built.kmer != params.salmon_kmer) stale << "k ${built.kmer}, not ${params.salmon_kmer}"
    if ((built.index_args ?: '') != (params.salmon_index_args ?: '')) stale << "index args '${built.index_args}'"
    def prov = file("${bundleDir()}/reference_provenance.json")
    if (prov.exists()) {
        def files = new groovy.json.JsonSlurper().parse(prov.toFile()).files
        def rna_sha = files.find { it.file == 'rna.fna' }?.sha256
        def genome_sha = files.find { it.file == 'genomic.fna' }?.sha256
        if (rna_sha && built.transcripts_sha256 != rna_sha) stale << 'a different rna.fna'
        if (built.decoys && genome_sha && built.genome_sha256 != genome_sha) stale << 'a different genomic.fna'
    }
    if (stale) {
        error "The Salmon index ${index} was built with ${stale.join(', ')} than this run uses. " +
              "Delete ${index} to rebuild it."
    }
}

def gitCommit() {
    try {
        def proc = ['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute()
        proc.waitFor()
        return proc.exitValue() == 0 ? proc.text.trim() : null
    } catch (Exception ignored) {
        return null
    }
}

// Run `whale record-rnaseq` in the checkout the run was launched from. A
// failure is reported but does not fail the run; it can be repeated by hand.
def recordRun() {
    def results = "${params.outdir}/${params.study_id}/rnaseq"
    def cmd = [params.whale_cmd, 'record-rnaseq', '--results', results.toString(), '--run-name', workflow.runName]
    log.info "Recording the run: ${cmd.join(' ')}"
    try {
        def pb = new ProcessBuilder(cmd).directory(workflow.launchDir.toFile()).redirectErrorStream(true)
        def proc = pb.start()
        proc.inputStream.eachLine { log.info "  ${it}" }
        if (proc.waitFor() != 0) {
            log.warn "whale record-rnaseq failed (exit ${proc.exitValue()}); fix the cause and run: ${cmd.join(' ')}"
        }
    } catch (Exception e) {
        log.warn "Could not run ${params.whale_cmd} (${e.message}); is it on PATH in the head job? Then run: ${cmd.join(' ')}"
    }
}

workflow {
    main:
    if (!params.study_id) error "params.study_id is not set (study config)"
    def samplesheet = file(params.rnaseq_samplesheet ?: "data/studies/${params.study_id}/rnaseq/samplesheet.tsv")
    if (!samplesheet.exists()) {
        error "No RNA-seq sample sheet at ${samplesheet}; build it with whale fetch-rnaseq-samplesheet (docs/rnaseq.md)."
    }
    def transcripts = bundleFile('rna.fna')
    def tx2gene     = bundleFile('tx2gene.tsv')
    def genome      = decoysOn() ? bundleFile('genomic.fna') : []
    checkCachedIndex()

    def rows = samplesheet.splitCsv(header: true, sep: '\t')
    def ids = rows*.sample_id
    def dupes = ids.findAll { id -> ids.count(id) > 1 }.unique()
    if (dupes) error "sample_id repeated in ${samplesheet}: ${dupes}"

    reads_ch = Channel
        .fromList(rows)
        .map { row ->
            def reads = [resolveFastq(row.fastq_1), resolveFastq(row.fastq_2)].findAll { it != null }
            if (!reads) error "no FASTQ for sample ${row.sample_id}"
            tuple(row.sample_id, reads)
        }

    FASTQC(reads_ch)
    FASTP(reads_ch)
    SALMON_INDEX(transcripts, genome, indexName())
    def index = SALMON_INDEX.out.index   // a value channel: every input is a value
    SALMON_QUANT(FASTP.out.reads, index)

    quant_dirs = SALMON_QUANT.out.quant.map { sid, dir -> dir }.collect()
    GENE_MATRICES(quant_dirs, samplesheet, tx2gene)
    RNASEQ_QC(FASTP.out.json.map { sid, f -> f }.mix(SALMON_QUANT.out.quant.map { sid, dir -> dir }).collect(),
              GENE_MATRICES.out.dir, samplesheet)
    MULTIQC(FASTQC.out.zip.map { sid, f -> f }.flatten()
            .mix(FASTP.out.json.map { sid, f -> f }, SALMON_QUANT.out.quant.map { sid, dir -> dir })
            .collect())

    versions_ch = FASTQC.out.versions
        .mix(FASTP.out.versions, SALMON_QUANT.out.versions, MULTIQC.out.versions)
        .collect()
    def run_params = groovy.json.JsonOutput.toJson([
        study_id: params.study_id, workflow: 'rnaseq', genome_assembly: params.genome_assembly,
        salmon_decoys: decoysOn(),
        salmon_kmer: params.salmon_kmer, salmon_index_args: params.salmon_index_args,
        salmon_quant_args: params.salmon_quant_args, fastp_args: params.fastp_args,
        se_adapter: params.se_adapter, tag_seq: params.tag_seq,
        rnaseq_min_mapping: params.rnaseq_min_mapping, rnaseq_min_mapped: params.rnaseq_min_mapped,
        rnaseq_min_genes: params.rnaseq_min_genes, profile: workflow.profile,
        nextflow_version: nextflow.version.toString(), run_name: workflow.runName,
        commit: workflow.commitId ?: gitCommit(), session: workflow.sessionId.toString(),
    ])
    MANIFEST(versions_ch, samplesheet, transcripts, tx2gene, index, GENE_MATRICES.out.dir, run_params)

    onComplete:
    if (workflow.success && params.record_status && !workflow.stubRun) {
        recordRun()
    }
}
