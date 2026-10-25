#!/usr/bin/env nextflow
/*
 * White Whale :: C/T SNP mask (PLAN.md §6.1, docs/snp_mask.md).
 *
 * A C/T SNP reads like an unmethylated C in bisulfite data, so CpGs on SNPs
 * are masked before methylation is compared. The SNPs are called from the
 * WGBS libraries themselves, pooled:
 *
 *   SORT_BAM (each library) -> PILEUP (BISCUIT, all libraries, per genome chunk)
 *     -> BUILD_MASK (snps.tsv.gz, cpg_mask.bed.gz, mask_manifest.json)
 *
 * Input: deduplicated Bismark BAMs, as a sheet (params.bam_sheet: sample_id,
 * study_id, bam) or a glob (params.bams, all from params.bam_study). The
 * methylation workflow keeps them with --publish_bams true.
 *
 *   nextflow run workflows/snp_mask -profile klone --mask_id <ID> --bam_sheet bams.tsv
 *   nextflow run workflows/snp_mask -profile test -stub
 *
 * Then `whale record-mask --dir <mask dir>` uploads it to the data host and
 * registers it in data/reference/snp_masks.yaml.
 */

nextflow.enable.dsl = 2

include { FAIDX; SORT_BAM; PILEUP; BUILD_MASK } from '../../modules/snp_mask/snp_mask.nf'

def gitCommit() {
    try {
        def proc = ['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute()
        proc.waitFor()
        return proc.exitValue() == 0 ? proc.text.trim() : null
    } catch (Exception ignored) {
        return null
    }
}

// Genome chunks for BISCUIT's -g: 1-based, inclusive.
def chunks(String name, long length, long size) {
    def n = (length + size - 1).intdiv(size) as int
    return (0..<n).collect { i ->
        def start = (i as long) * size
        "${name}:${start + 1}-${Math.min(start + size, length)}".toString()
    }
}

workflow {
    main:
    if (!params.mask_id) error "params.mask_id is not set, e.g. --mask_id GCF_963853765.1.snpmask.v1"
    if (!params.bam_sheet && !params.bams) error "give the BAMs: --bam_sheet <tsv: sample_id, study_id, bam> or --bams '<glob>'"

    def genome_fasta = file("${params.reference_dir}/${params.genome_assembly}/genomic.fna", checkIfExists: true)

    def samples
    if (params.bam_sheet) {
        samples = file(params.bam_sheet, checkIfExists: true).splitCsv(header: true, sep: '\t')
            .collect { r -> [sample_id: r.sample_id, study_id: r.study_id, bam: file(r.bam, checkIfExists: true)] }
    } else {
        samples = files(params.bams, checkIfExists: true)
            .collect { f -> [sample_id: f.name.tokenize('.')[0], study_id: params.bam_study, bam: f] }
    }
    def ids = samples*.sample_id
    if (ids.size() != ids.unique(false).size()) error "sample_id repeats in the BAM list: ${ids.countBy { it }.findAll { it.value > 1 }.keySet()}"

    FAIDX(genome_fasta)
    SORT_BAM(Channel.fromList(samples.collect { s -> tuple(s.sample_id, s.bam) }))
    regions = FAIDX.out.genome.flatMap { fa, fai ->
        fai.readLines().collectMany { line ->
            def f = line.tokenize('\t')
            chunks(f[0], f[1] as long, params.snp_region_size as long)
        }
    }
    bams = SORT_BAM.out.bam.flatMap { bam, bai -> [bam, bai] }.collect()
    PILEUP(regions, FAIDX.out.genome.first(), bams)

    def mask_params = groovy.json.JsonOutput.toJson([
        mask_id: params.mask_id, genome_assembly: params.genome_assembly,
        samples: samples.collect { s -> [sample_id: s.sample_id, study_id: s.study_id, bam: s.bam.name] },
        studies: samples*.study_id.unique(false),
        min_mapq: params.snp_min_mapq, min_baseq: params.snp_min_baseq,
        prefilter_alt_reads: params.snp_prefilter_alt, region_size: params.snp_region_size,
        nextflow_version: nextflow.version.toString(), run_name: workflow.runName,
        commit: workflow.commitId ?: gitCommit(), profile: workflow.profile,
    ])
    BUILD_MASK(PILEUP.out.sites.collect(), FAIDX.out.genome.first(), PILEUP.out.versions.first(), mask_params)
}
