process BED_TO_REGION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(bed), val(scatter_count)

    output:
    tuple val(meta), path("*.region.txt"), val(scatter_count), emit: region
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - <<'PY'
from pathlib import Path
chrom = None
start = None
end = None
for line in Path("${bed}").read_text().splitlines():
    if not line.strip() or line.startswith(('#', 'track', 'browser')):
        continue
    parts = line.split()
    if len(parts) < 3:
        continue
    c, s, e = parts[0], int(parts[1]), int(parts[2])
    if chrom is None:
        chrom = c
        start = s
        end = e
    else:
        start = min(start, s)
        end = max(end, e)
if chrom is None:
    chrom, start, end = "1", 1, 2
Path("${prefix}.region.txt").write_text(f"{chrom}:{start}-{end}\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '1:1-2' > ${prefix}.region.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
