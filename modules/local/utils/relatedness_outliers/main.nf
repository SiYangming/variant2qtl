process RELATEDNESS_OUTLIERS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("*.related_outliers.txt"), emit: outliers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pi_hat = task.ext.pi_hat ?: params.genotype_qc_pi_hat ?: 0.185
    """
    python3 - <<'PY'
from pathlib import Path
genome = Path("${genome}")
pi_hat = float("${pi_hat}")
remove = set()
with genome.open() as fh:
    header = fh.readline().split()
    try:
        i_pi = header.index("PI_HAT")
        i_fid1 = header.index("FID1")
        i_iid1 = header.index("IID1")
        i_fid2 = header.index("FID2")
        i_iid2 = header.index("IID2")
    except ValueError:
        # fallback positional PLINK1.9 genome columns
        i_fid1, i_iid1, i_fid2, i_iid2, i_pi = 0, 1, 2, 3, 9
    for line in fh:
        parts = line.split()
        if len(parts) <= max(i_fid1, i_iid1, i_fid2, i_iid2, i_pi):
            continue
        try:
            pi = float(parts[i_pi])
        except ValueError:
            continue
        if pi >= pi_hat:
            # drop second sample of each related pair
            remove.add((parts[i_fid2], parts[i_iid2]))
out = Path("${prefix}.related_outliers.txt")
out.write_text("".join(f"{fid}\\t{iid}\\n" for fid, iid in sorted(remove)))
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
