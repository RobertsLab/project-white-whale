// The C/T SNP mask (PLAN.md §6.1): sort and index each library's deduplicated
// Bismark BAM, pile all of them up together with BISCUIT in genome chunks,
// keep sites with alternative-allele evidence (pileup_sites.awk), and turn the
// SNPs into a mask of reference CpGs (snp_mask.py). See docs/snp_mask.md.

process FAIDX {
    tag "${fasta.simpleName}"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    input:
    path fasta, stageAs: 'genome.fa'

    output:
    tuple path('genome.fa', includeInputs: true), path('genome.fa.fai'), emit: genome

    script:
    """
    samtools faidx genome.fa
    """

    stub:
    """
    samtools faidx genome.fa 2>/dev/null || awk '/^>/{if(n)print n"\\t"l"\\t0\\t0\\t0"; n=substr(\$1,2); l=0; next}{l+=length(\$0)} END{print n"\\t"l"\\t0\\t0\\t0"}' genome.fa > genome.fa.fai
    """
}

process SORT_BAM {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -m 1500M -T tmp_${sample_id} -o ${sample_id}.sorted.bam ${bam}
    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
    """
}

// One genome chunk, every library at once. Only sites with >= snp_prefilter_alt
// alternative reads leave the task, so the multi-sample VCF never touches disk.
process PILEUP {
    tag "${region}"
    label 'process_high'
    container 'quay.io/biocontainers/biscuit:1.10.2.20260818--hdf5d79d_0'

    input:
    val region
    tuple path(genome), path(fai)
    path bams, stageAs: 'bams/*'

    output:
    path "*.sites.tsv.gz", emit: sites
    path "versions.yml", emit: versions

    script:
    def name = region.replaceAll(/[:-]/, '_')
    """
    set -o pipefail
    biscuit pileup -g ${region} -@ ${task.cpus} -m ${params.snp_min_mapq} -b ${params.snp_min_baseq} \\
        -w stats genome.fa \$(ls bams/*.bam | sort) \\
      | awk -v min_alt=${params.snp_prefilter_alt} -f \$(command -v pileup_sites.awk) \\
      | gzip > ${name}.sites.tsv.gz
    printf 'PILEUP:\\n  biscuit: "%s"\\n' "\$(biscuit version 2>&1 | grep -Eo '[0-9][0-9.]+' | head -1)" > versions.yml
    """

    stub:
    def name = region.replaceAll(/[:-]/, '_')
    """
    printf '' | gzip > ${name}.sites.tsv.gz
    printf 'PILEUP:\\n  biscuit: "stub"\\n' > versions.yml
    """
}

process BUILD_MASK {
    label 'process_medium'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.mask_outdir ?: "${params.reference_dir}/${params.genome_assembly}/snp_mask/${params.mask_id}"}", mode: 'copy'

    input:
    path sites, stageAs: 'sites/*'
    tuple path(genome), path(fai)
    path versions
    val mask_params

    output:
    path "snps.tsv.gz", emit: snps
    path "cpg_mask.bed.gz", emit: mask
    path "mask_manifest.json", emit: manifest

    script:
    """
    cat <<'EOF_PARAMS' > params.json
    ${mask_params}
    EOF_PARAMS
    python3 - <<'EOF_PY'
    import json, pathlib
    p = json.loads(pathlib.Path('params.json').read_text())
    tools = {}
    for line in pathlib.Path('${versions}').read_text().splitlines():
        if ':' in line and line.startswith(' '):
            k, v = line.strip().split(':', 1)
            tools[k] = v.strip().strip('"')
    p['tool_versions'] = tools
    pathlib.Path('params.json').write_text(json.dumps(p))
    EOF_PY
    snp_mask.py --genome genome.fa --params params.json \\
        --min-alt ${params.snp_min_alt} --min-af ${params.snp_min_af} sites/*
    """

    stub:
    """
    cat <<'EOF_PARAMS' > params.json
    ${mask_params}
    EOF_PARAMS
    snp_mask.py --genome genome.fa --params params.json sites/*
    """
}
