#!/usr/bin/env nextflow
/*
 * White Whale :: methylation workflow — FASTQ in, per-sample methylome out.
 *
 * Still to come (PLAN.md §5): the C/T SNP mask, chunked alignment for klone,
 * 1 kb tiles and the study report.
 *
 *   FASTQC -> TRIM_GALORE -> BISMARK_ALIGN -> BISMARK_DEDUPLICATE (not RRBS)
 *     -> BISMARK_METHYLATION_EXTRACTOR -> CPG_SUMMARY -> QC_GATE
 *                                      -> FEATURE_SUMMARY (genes) -> MANIFEST
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
 * Libraries in params.swap_mates were deposited with read 1 and read 2
 * exchanged; their fastq_2 is used as read 1 (the samplesheet stays as ENA
 * has it).
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
include { FEATURE_SUMMARY }               from '../../modules/methylation/features.nf'
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

// A library in params.swap_mates, read 1 and read 2 exchanged. Nextflow's task
// hash does not see the order of a task's input files, so the pair is linked
// under its own names (work/swap_mates/<id>/<id>_swapped_{1,2}.fastq.gz): its
// tasks then miss the -resume cache of a run that read it unswapped.
def swapMates(String sample_id, List reads) {
    def dir = file("${workflow.workDir}/swap_mates/${sample_id}")
    dir.mkdirs()
    return [reads[1], reads[0]].withIndex().collect { f, i ->
        f.mklink(dir.resolve("${sample_id}_swapped_${i + 1}.fastq.gz"), overwrite: true)
    }
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
    def features     = file("${params.reference_dir}/${params.genome_assembly}/features", checkIfExists: true)
    // C/T SNP mask (docs/snp_mask.md): a mask folder with cpg_mask.bed.gz and mask_manifest.json.
    def snp_mask = null
    def mask_bed = []
    if (params.snp_mask) {
        mask_bed = file("${params.snp_mask}/cpg_mask.bed.gz", checkIfExists: true)
        def mm = new groovy.json.JsonSlurper().parse(file("${params.snp_mask}/mask_manifest.json", checkIfExists: true).toFile())
        if (mm.params.genome_assembly != params.genome_assembly) {
            error "SNP mask ${mm.params.mask_id} is for ${mm.params.genome_assembly}, not ${params.genome_assembly}"
        }
        snp_mask = [mask_id: mm.params.mask_id, sha256: mm.outputs['cpg_mask.bed.gz'], cpgs_masked: mm.counts.cpgs_masked]
    }

    def sample_ids = samplesheet.splitCsv(header: true, sep: '\t')*.sample_id
    def unknown = params.swap_mates.findAll { !(it in sample_ids) }
    if (unknown) error "swap_mates lists samples not in the sample sheet: ${unknown}"

    reads_ch = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def reads = [resolveFastq(row.fastq_1), resolveFastq(row.fastq_2)].findAll { it != null }
            if (!reads) error "no FASTQ for sample ${row.sample_id}"
            if (row.sample_id in params.swap_mates) {
                if (reads.size() != 2) error "swap_mates lists ${row.sample_id}, which is not paired-end"
                reads = swapMates(row.sample_id, reads)
            }
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
    FEATURE_SUMMARY(BISMARK_METHYLATION_EXTRACTOR.out.cpg, features, mask_bed)

    versions_ch = FASTQC.out.versions
        .mix(TRIM_GALORE.out.versions, BISMARK_ALIGN.out.versions, BISMARK_METHYLATION_EXTRACTOR.out.versions)
        .collect()

    def run_params = groovy.json.JsonOutput.toJson([
        study_id: params.study_id, part: params.part, chemistry: params.chemistry,
        genome_assembly: params.genome_assembly, bismark_align_args: params.bismark_align_args,
        chemistry_align_args: chem.align,
        extractor_args: params.extractor_args, trim_args: params.trim_args,
        swap_mates: params.swap_mates,
        qc_min_mapping: params.qc_min_mapping, qc_min_conversion: params.qc_min_conversion,
        qc_min_cpgs_5x: params.qc_min_cpgs_5x, snp_mask: snp_mask,
        deduplicate: chem.deduplicate, profile: workflow.profile,
        nextflow_version: nextflow.version.toString(), run_name: workflow.runName,
        commit: workflow.commitId ?: gitCommit(), session: workflow.sessionId.toString(),
    ])
    MANIFEST(versions_ch, samplesheet, genome_fasta, features, run_params)

    onComplete:
    if (workflow.success && params.record_status && !workflow.stubRun) {
        recordRun()
    }
}
