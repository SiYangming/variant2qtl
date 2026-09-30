process HET_OUTLIERS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(het)

    output:
    tuple val(meta), path("*.het_outliers.txt"), emit: outliers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sd_thresh = task.ext.sd_thresh ?: params.genotype_qc_het_sd ?: 3.0
    """
    python3 - <<'PY'
from pathlib import Path
het_path = Path("${het}")
sd_thresh = float("${sd_thresh}")
rows = []
with het_path.open() as fh:
    header = fh.readline()
    for line in fh:
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            f = float(parts[5])
        except ValueError:
            continue
        rows.append((parts[0], parts[1], f))
out = Path("${prefix}.het_outliers.txt")
if not rows:
    out.write_text("")
else:
    vals = [r[2] for r in rows]
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / max(len(vals) - 1, 1)
    sd = var ** 0.5
    keep = []
    if sd > 0:
        for fid, iid, f in rows:
            if abs(f - mean) > sd_thresh * sd:
                keep.append(f"{fid}\\t{iid}\\n")
    out.write_text("".join(keep))
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.het_outliers.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
