// The run manifest: parameters, tool versions and input checksums, committed
// to data/studies/<STUDY_ID>/ on the study branch (PLAN.md §2.5).
process MANIFEST {
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'
    publishDir "${params.outdir}/${params.study_id}", mode: 'copy'

    input:
    path versions, stageAs: "versions*.yml"
    path samplesheet
    path genome_fasta
    val run_params

    output:
    path "run_manifest.json", emit: manifest

    script:
    """
    cat <<'EOF_PARAMS' > params.json
    ${run_params}
    EOF_PARAMS
    python3 - <<'EOF_PY'
    import hashlib, json, pathlib
    def sha256(p):
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        return h.hexdigest()
    versions = {}
    for f in sorted(pathlib.Path('.').glob('versions*.yml')):
        proc = None
        for line in f.read_text().splitlines():
            if line and not line.startswith(' '):
                proc = line.rstrip(':')
            elif proc and ':' in line:
                k, v = line.strip().split(':', 1)
                versions.setdefault(proc, {})[k] = v.strip().strip('"')
    manifest = {
        'params': json.loads(pathlib.Path('params.json').read_text()),
        'tool_versions': versions,
        'inputs': {
            'samplesheet': {'name': '${samplesheet.name}', 'sha256': sha256('${samplesheet}')},
            'genome_fasta': {'name': '${genome_fasta.name}', 'sha256': sha256('${genome_fasta}')},
        },
    }
    pathlib.Path('run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\\n')
    EOF_PY
    """
}
