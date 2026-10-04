process SOMALIER_OUTLIERS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path("*.related_outliers.txt"), emit: outliers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def relatedness = task.ext.relatedness ?: params.genotype_qc_somalier_relatedness ?: 0.2
    """
    python3 - <<'PY'
import csv
from pathlib import Path
pairs = Path("${pairs}")
thr = float("${relatedness}")
remove = set()
with pairs.open() as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    fields = {(k or "").lower(): k for k in (reader.fieldnames or [])}
    # somalier pairs.tsv: #sample_a sample_b relatedness ...
    a_key = fields.get("#sample_a") or fields.get("sample_a") or fields.get("sample1")
    b_key = fields.get("sample_b") or fields.get("sample2")
    r_key = fields.get("relatedness") or fields.get("relatedness_ibs0")
    if not a_key:
        # fallback: first three columns
        fh.seek(0)
        next(fh, None)
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) < 3:
                continue
            try:
                rel = float(parts[2])
            except ValueError:
                continue
            if rel >= thr:
                remove.add(parts[1])
    else:
        for row in reader:
            try:
                rel = float(row[r_key])
            except (KeyError, TypeError, ValueError):
                continue
            if rel >= thr:
                # drop second sample of each related pair
                remove.add(row[b_key])
out = Path("${prefix}.related_outliers.txt")
# FID IID (assume FID==IID for somalier sample ids)
out.write_text("".join(f"{s}\\t{s}\\n" for s in sorted(remove)))
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.related_outliers.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
