#!/usr/bin/env nextflow
/*
 * White Whale :: enrichment workflow — MeDIP and MBD-seq libraries in, per-sample
 * relative enrichment out (PLAN.md §4.1, §9 Phase 4; docs/enrichment.md).
 *
 * Enrichment assays capture methylated fragments; they give read density over
 * methylated regions, not a per-CpG methylation proportion. So this workflow
 * aligns plain reads (no bisulfite conversion), counts fragments, and reports
 * enrichment normalised for library size and CpG density. Its outputs never
 * enter the base-resolution atlas tables (whale build-atlas leaves them out).
 *
 *   FASTQC -> TRIM_GALORE -> BOWTIE2_ALIGN -> MARKDUP (-r, MAPQ filter)
 *     -> ENRICHMENT_COUNT (1 kb tiles, gene bodies, promoters; QC: relH, GoGe,
 *        CpG coverage, saturation) -> ENRICHMENT_QC_GATE -> MANIFEST
 *
 * When a run on a study branch succeeds it records itself with `whale
 * record-enrichment` (qc_summary.tsv, run_manifest.json, gene x sample
 * enrichment matrices and a study summary marked tier: enrichment; per-sample
 * tables to the data host). Turn it off with --record_status false.
 *
 * Sample sheet and FASTQ locations as for workflows/methylation.
 *
 *   nextflow run workflows/enrichment -c config/<STUDY_ID>.config -profile raven
 *   nextflow run workflows/enrichment -profile test -stub        # DAG smoke test
 *   nextflow run workflows/enrichment -profile test,docker       # every tool, simulated MeDIP reads
 */

nextflow.enable.dsl = 2

include { FASTQC }             from '../../modules/methylation/fastqc.nf'
include { TRIM_GALORE }        from '../../modules/methylation/trim_galore.nf'
include { MANIFEST }           from '../../modules/methylation/manifest.nf'
include { BOWTIE2_BUILD }      from '../../modules/enrichment/bowtie2.nf'
include { BOWTIE2_ALIGN }      from '../../modules/enrichment/bowtie2.nf'
include { MARKDUP }            from '../../modules/enrichment/dedup.nf'
include { ENRICHMENT_COUNT }   from '../../modules/enrichment/counts.nf'
include { ENRICHMENT_QC_GATE } from '../../modules/enrichment/counts.nf'

// Mirrors registry/controlled_vocabularies/assay_chemistry.yaml (tier: enrichment).
def chemistries() {
    return [medip: [deduplicate: true], mbd: [deduplicate: true]]
}

// What the counts mean; copied into run_manifest.json and the study summary.
def enrichmentMethod() {
    return [
        name: 'cpg_density_normalised_window_counts',
        version: '1',
        counting: 'deduplicated fragments with MAPQ >= min_mapq (proper pairs; single-end reads extended to fragment_length), each counted once at its midpoint',
        windows: "${params.window_size} bp genome tiles from position 0 of each sequence",
        features: 'gene_body (genes.bed) and promoter (promoters.bed) of the reference bundle, intervals merged per gene',
        cpm: 'count * 1e6 / fragments counted in the library',
        obs_exp: '(count / fragments) / (region CpGs / genome CpGs); empty below min_cpgs',
        not: 'not a methylation proportion; no MEDIPS/QSEA calibration (mammalian CpG-density curve)',
    ]
}

def resolveFastq(String p) {
    if (!p) return null
    def f = p.startsWith('/') ? file(p) : file("${params.raw_dir}/${params.study_id}/${p}")
    if (!f.exists()) error "FASTQ not found: ${f}"
    return f
}

// A cached bowtie2 index must have been built from the bundle's current genomic.fna.
def checkCachedIndex() {
    def index = file("${params.bowtie2_index_dir ?: "${params.reference_dir}/${params.genome_assembly}"}/bowtie2_genome")
    if (index.resolve('STUB').exists() && !workflow.stubRun) {
        error "${index} is an empty index left by a -stub run; delete it to build a real one."
    }
    def built = index.resolve('genome.sha256')
    def prov  = file("${params.reference_dir}/${params.genome_assembly}/reference_provenance.json")
    if (!built.exists() || !prov.exists()) return
    def want = new groovy.json.JsonSlurper().parse(prov.toFile()).files.find { it.file == 'genomic.fna' }?.sha256
    if (want && built.text.trim() != want) {
        error "The bowtie2 index ${index} was built from a different genomic.fna than " +
              "${params.genome_assembly}'s reference bundle now holds. Delete ${index} to rebuild it."
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

// `whale record-enrichment` in the launch checkout. A failure is reported but
// does not fail the run; it can be repeated by hand.
def recordRun() {
    def results = "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}"
    def cmd = [params.whale_cmd, 'record-enrichment', '--results', results.toString(), '--run-name', workflow.runName]
    log.info "Recording the run: ${cmd.join(' ')}"
    try {
        def pb = new ProcessBuilder(cmd).directory(workflow.launchDir.toFile()).redirectErrorStream(true)
        def proc = pb.start()
        proc.inputStream.eachLine { log.info "  ${it}" }
        if (proc.waitFor() != 0) {
            log.warn "whale record-enrichment failed (exit ${proc.exitValue()}); fix the cause and run: ${cmd.join(' ')}"
        }
    } catch (Exception e) {
        log.warn "Could not run ${params.whale_cmd} (${e.message}); is it on PATH in the head job? Then run: ${cmd.join(' ')}"
    }
}

workflow {
    main:
    if (!params.study_id)    error "params.study_id is not set (study config)"
    if (!params.samplesheet) error "params.samplesheet is not set (study config)"
    def known = chemistries()
    if (!known.containsKey(params.chemistry)) {
        error "params.chemistry must be one of ${known.keySet()} (base-resolution assays use workflows/methylation); got ${params.chemistry}"
    }

    def genome_fasta = file("${params.reference_dir}/${params.genome_assembly}/genomic.fna", checkIfExists: true)
    checkCachedIndex()
    def samplesheet  = file(params.samplesheet, checkIfExists: true)
    def features     = file("${params.reference_dir}/${params.genome_assembly}/features", checkIfExists: true)

    reads_ch = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def reads = [resolveFastq(row.fastq_1), resolveFastq(row.fastq_2)].findAll { it != null }
            if (!reads) error "no FASTQ for sample ${row.sample_id}"
            tuple(row.sample_id, reads)
        }

    FASTQC(reads_ch)
    TRIM_GALORE(reads_ch, params.trim_args ?: '')
    BOWTIE2_BUILD(genome_fasta)
    BOWTIE2_ALIGN(TRIM_GALORE.out.reads, BOWTIE2_BUILD.out.index.first())
    MARKDUP(BOWTIE2_ALIGN.out.bam)
    ENRICHMENT_COUNT(MARKDUP.out.fragments, genome_fasta, features)

    qc_inputs = BOWTIE2_ALIGN.out.log
        .mix(MARKDUP.out.stats, ENRICHMENT_COUNT.out.qc)
        .collect()
    ENRICHMENT_QC_GATE(qc_inputs, samplesheet)

    versions_ch = FASTQC.out.versions
        .mix(TRIM_GALORE.out.versions, BOWTIE2_ALIGN.out.versions)
        .collect()

    def run_params = groovy.json.JsonOutput.toJson([
        study_id: params.study_id, part: params.part, chemistry: params.chemistry, tier: 'enrichment',
        workflow: 'enrichment',
        genome_assembly: params.genome_assembly, bowtie2_args: params.bowtie2_args, trim_args: params.trim_args,
        min_mapq: params.min_mapq, window_size: params.window_size, fragment_length: params.fragment_length,
        enrichment_min_cpgs: params.enrichment_min_cpgs, saturation_seed: params.saturation_seed,
        enrichment_method: enrichmentMethod(),
        qc_min_mapping: params.enrichment_min_mapping, qc_min_fragments: params.enrichment_min_fragments,
        qc_min_relh: params.enrichment_min_relh,
        deduplicate: known[params.chemistry].deduplicate, profile: workflow.profile,
        nextflow_version: nextflow.version.toString(), run_name: workflow.runName,
        commit: workflow.commitId ?: gitCommit(), session: workflow.sessionId.toString(),
    ])
    MANIFEST(versions_ch, samplesheet, genome_fasta, features, run_params)

    onComplete:
    if (workflow.success && params.record_status && !workflow.stubRun) {
        recordRun()
    }
}
