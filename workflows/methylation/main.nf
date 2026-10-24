#!/usr/bin/env nextflow
/*
 * White Whale :: methylation workflow — FASTQ in, per-sample methylome out.
 *
 * Still to come (PLAN.md §5): the C/T SNP mask, chunked alignment for klone,
 * per-sample feature summaries and the study report.
 *
 *   FASTQC -> TRIM_GALORE -> BISMARK_ALIGN -> BISMARK_DEDUPLICATE (not RRBS)
 *     -> BISMARK_METHYLATION_EXTRACTOR -> CPG_SUMMARY -> QC_GATE -> MANIFEST
 *
 * When a run on a study branch succeeds, it records itself: `whale
 * record-run` commits qc_summary.tsv, run_manifest.json, the methylome
 * manifest and the study summary, uploads the methylomes to the data host,
 * moves status.yaml on (aligned, called, qc_passed, summarized) and pushes.
 * The dashboard then shows the stage and the study's findings. Turn it off
 * with --record_status false.
 *
 * The sample sheet (data/studies/<STUDY_ID>/samplesheet.tsv) has columns
 * sample_id, run_accession, fastq_1, fastq_2 (empty for single-end). FASTQ
 * paths are relative to ${params.raw_dir}/<STUDY_ID>/ unless absolute.
 *
 *   nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile raven,docker
 *   nextflow run workflows/methylation -profile test -stub        # smoke test
 *   nextflow run workflows/methylation -profile test,docker       # every tool, simulated reads (make real-run)
 */

nextflow.enable.dsl = 2

include { FASTQC }                        from '../../modules/methylation/fastqc.nf'
include { TRIM_GALORE }                   from '../../modules/methylation/trim_galore.nf'
include { BISMARK_GENOME_PREPARATION }    from '../../modules/methylation/bismark.nf'
include { BISMARK_ALIGN }                 from '../../modules/methylation/bismark.nf'
include { BISMARK_DEDUPLICATE }           from '../../modules/methylation/bismark.nf'
include { BISMARK_METHYLATION_EXTRACTOR } from '../../modules/methylation/bismark.nf'
include { CPG_SUMMARY }                   from '../../modules/methylation/qc.nf'
include { QC_GATE }                       from '../../modules/methylation/qc.nf'
include { MANIFEST }                      from '../../modules/methylation/manifest.nf'

// Mirrors registry/controlled_vocabularies/assay_chemistry.yaml (base-resolution terms).
def chemistries() {
    return [
        wgbs : [deduplicate: true,  trim: '',       align: ''],
        emseq: [deduplicate: true,  trim: '',       align: ''],
        pbat : [deduplicate: true,  trim: '',       align: '--pbat'],
        rrbs : [deduplicate: false, trim: '--rrbs', align: ''],
    ]
}

def resolveFastq(String p) {
    if (!p) return null
    def f = p.startsWith('/') ? file(p) : file("${params.raw_dir}/${params.study_id}/${p}")
    if (!f.exists()) error "FASTQ not found: ${f}"
    return f
}

// A cached Bismark index must have been built from the bundle's current
// genomic.fna (fetch-reference may have replaced it since).
def checkCachedIndex() {
    def index = file("${params.bismark_index_dir ?: "${params.reference_dir}/${params.genome_assembly}"}/bismark_genome")
    def built = index.resolve('genome.sha256')
    def prov  = file("${params.reference_dir}/${params.genome_assembly}/reference_provenance.json")
    if (!built.exists() || !prov.exists()) return
    def want = new groovy.json.JsonSlurper().parse(prov.toFile()).files.find { it.file == 'genomic.fna' }?.sha256
    def have = built.text.trim()
    if (want && have != want) {
        error "The Bismark index ${index} was built from a different genomic.fna than " +
              "${params.genome_assembly}'s reference bundle now holds. Delete ${index} to rebuild it."
    }
}

// The commit the pipeline code came from. workflow.commitId is only set for
// runs of a remote repository (nextflow run <github repo>), never for a run
// from a checkout, which is how every study runs.
def gitCommit() {
    try {
        def proc = ['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute()
        proc.waitFor()
        return proc.exitValue() == 0 ? proc.text.trim() : null
    } catch (Exception ignored) {
        return null
    }
}

// Run `whale record-run` in the checkout the run was launched from. A failure
// here is reported but does not fail the run: its outputs are complete, and
// record-run can be repeated by hand.
def recordRun() {
    def results = "${params.outdir}/${params.study_id}${params.part ? '/parts/' + params.part : ''}"
    def cmd = [params.whale_cmd, 'record-run', '--results', results.toString(), '--run-name', workflow.runName]
    log.info "Recording the run: ${cmd.join(' ')}"
    try {
        def pb = new ProcessBuilder(cmd).directory(workflow.launchDir.toFile()).redirectErrorStream(true)
        def proc = pb.start()
        proc.inputStream.eachLine { log.info "  ${it}" }
        if (proc.waitFor() != 0) {
            log.warn "whale record-run failed (exit ${proc.exitValue()}); fix the cause and run: ${cmd.join(' ')}"
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
        error "params.chemistry must be one of ${known.keySet()} (enrichment assays use their own workflow); got ${params.chemistry}"
    }
    def chem = known[params.chemistry]

    def genome_fasta = file("${params.reference_dir}/${params.genome_assembly}/genomic.fna", checkIfExists: true)
    checkCachedIndex()
    def samplesheet  = file(params.samplesheet, checkIfExists: true)

    reads_ch = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def reads = [resolveFastq(row.fastq_1), resolveFastq(row.fastq_2)].findAll { it != null }
            if (!reads) error "no FASTQ for sample ${row.sample_id}"
            tuple(row.sample_id, reads)
        }

    FASTQC(reads_ch)
    TRIM_GALORE(reads_ch, [chem.trim, params.trim_args].findAll { it }.join(' '))
    BISMARK_GENOME_PREPARATION(genome_fasta)
    BISMARK_ALIGN(TRIM_GALORE.out.reads, BISMARK_GENOME_PREPARATION.out.genome_dir.first(),
                  [chem.align, params.bismark_align_args].findAll { it }.join(' '))

    def bam_ch
    def dedup_reports
    if (chem.deduplicate) {
        BISMARK_DEDUPLICATE(BISMARK_ALIGN.out.bam)
        bam_ch = BISMARK_DEDUPLICATE.out.bam
        dedup_reports = BISMARK_DEDUPLICATE.out.report
    } else {
        bam_ch = BISMARK_ALIGN.out.bam
        dedup_reports = Channel.empty()
    }

    BISMARK_METHYLATION_EXTRACTOR(bam_ch, BISMARK_GENOME_PREPARATION.out.genome_dir.first())

    CPG_SUMMARY(BISMARK_METHYLATION_EXTRACTOR.out.cpg)
    qc_inputs = BISMARK_ALIGN.out.report
        .mix(dedup_reports,
             BISMARK_METHYLATION_EXTRACTOR.out.qc.map { sid, mbias, splitting -> splitting },
             CPG_SUMMARY.out.summary)
        .collect()
    QC_GATE(qc_inputs, samplesheet)

    versions_ch = FASTQC.out.versions
        .mix(TRIM_GALORE.out.versions, BISMARK_ALIGN.out.versions, BISMARK_METHYLATION_EXTRACTOR.out.versions)
        .collect()

    def run_params = groovy.json.JsonOutput.toJson([
        study_id: params.study_id, part: params.part, chemistry: params.chemistry,
        genome_assembly: params.genome_assembly, bismark_align_args: params.bismark_align_args,
        chemistry_align_args: chem.align,
        extractor_args: params.extractor_args, trim_args: params.trim_args,
        qc_min_mapping: params.qc_min_mapping, qc_min_conversion: params.qc_min_conversion,
        qc_min_cpgs_5x: params.qc_min_cpgs_5x,
        deduplicate: chem.deduplicate, profile: workflow.profile,
        nextflow_version: nextflow.version.toString(), run_name: workflow.runName,
        commit: workflow.commitId ?: gitCommit(), session: workflow.sessionId.toString(),
    ])
    MANIFEST(versions_ch, samplesheet, genome_fasta, run_params)

    onComplete:
    if (workflow.success && params.record_status && !workflow.stubRun) {
        recordRun()
    }
}
