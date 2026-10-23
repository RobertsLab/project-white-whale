#!/usr/bin/env nextflow
/*
 * White Whale :: methylation workflow — FASTQ in, per-sample methylome out.
 *
 * STATUS (Phase 0 scaffold): the DAG below runs end to end with -stub on the
 * test fixtures (CI). The real commands are written but have not yet been run
 * on a study. Still to come (PLAN.md §5): conversion efficiency, the C/T SNP
 * mask, chunked alignment for klone, per-sample feature summaries, the QC
 * gate, and the study report.
 *
 *   FASTQC -> TRIM_GALORE -> BISMARK_ALIGN -> BISMARK_DEDUPLICATE (not RRBS)
 *     -> BISMARK_METHYLATION_EXTRACTOR -> MANIFEST
 *
 * The sample sheet (data/studies/<STUDY_ID>/samplesheet.tsv) has columns
 * sample_id, run_accession, fastq_1, fastq_2 (empty for single-end). FASTQ
 * paths are relative to ${params.raw_dir}/<STUDY_ID>/ unless absolute.
 *
 *   nextflow run workflows/methylation -c config/<STUDY_ID>.config -profile raven,docker
 *   nextflow run workflows/methylation -profile test -stub        # smoke test
 */

nextflow.enable.dsl = 2

include { FASTQC }                        from '../../modules/methylation/fastqc.nf'
include { TRIM_GALORE }                   from '../../modules/methylation/trim_galore.nf'
include { BISMARK_GENOME_PREPARATION }    from '../../modules/methylation/bismark.nf'
include { BISMARK_ALIGN }                 from '../../modules/methylation/bismark.nf'
include { BISMARK_DEDUPLICATE }           from '../../modules/methylation/bismark.nf'
include { BISMARK_METHYLATION_EXTRACTOR } from '../../modules/methylation/bismark.nf'
include { MANIFEST }                      from '../../modules/methylation/manifest.nf'

// Mirrors registry/controlled_vocabularies/assay_chemistry.yaml (base-resolution terms).
def chemistries() {
    return [
        wgbs : [deduplicate: true,  trim: ''],
        emseq: [deduplicate: true,  trim: ''],
        pbat : [deduplicate: true,  trim: ''],
        rrbs : [deduplicate: false, trim: '--rrbs'],
    ]
}

def resolveFastq(String p) {
    if (!p) return null
    def f = p.startsWith('/') ? file(p) : file("${params.raw_dir}/${params.study_id}/${p}")
    if (!f.exists()) error "FASTQ not found: ${f}"
    return f
}

workflow {
    if (!params.study_id)    error "params.study_id is not set (study config)"
    if (!params.samplesheet) error "params.samplesheet is not set (study config)"
    def known = chemistries()
    if (!known.containsKey(params.chemistry)) {
        error "params.chemistry must be one of ${known.keySet()} (enrichment assays use their own workflow); got ${params.chemistry}"
    }
    def chem = known[params.chemistry]

    def genome_fasta = file("${params.reference_dir}/${params.genome_assembly}/genomic.fna", checkIfExists: true)
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
    BISMARK_ALIGN(TRIM_GALORE.out.reads, BISMARK_GENOME_PREPARATION.out.genome_dir.first())

    bam_ch = chem.deduplicate
        ? BISMARK_DEDUPLICATE(BISMARK_ALIGN.out.bam).bam
        : BISMARK_ALIGN.out.bam

    BISMARK_METHYLATION_EXTRACTOR(bam_ch, BISMARK_GENOME_PREPARATION.out.genome_dir.first())

    versions_ch = FASTQC.out.versions
        .mix(TRIM_GALORE.out.versions, BISMARK_ALIGN.out.versions, BISMARK_METHYLATION_EXTRACTOR.out.versions)
        .collect()

    def run_params = groovy.json.JsonOutput.toJson([
        study_id: params.study_id, part: params.part, chemistry: params.chemistry,
        genome_assembly: params.genome_assembly, bismark_align_args: params.bismark_align_args,
        extractor_args: params.extractor_args, trim_args: params.trim_args,
        deduplicate: chem.deduplicate, profile: workflow.profile,
        nextflow_version: nextflow.version.toString(), run_name: workflow.runName,
        commit: workflow.commitId ?: null, session: workflow.sessionId.toString(),
    ])
    MANIFEST(versions_ch, samplesheet, genome_fasta, run_params)
}
