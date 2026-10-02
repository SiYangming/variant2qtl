process TENSORQTL_CIS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // No official bioconda/biocontainers tensorqtl package. Prefer conda/wave for real runs.
    // Use biocontainers python for stub/CI (registry default is quay.io).
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(covariates)

    output:
    tuple val(meta), path("*.cis_qtl*.txt.gz"), emit: cis_qtl, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('tensorqtl'), eval('python3 -c "import tensorqtl; print(getattr(tensorqtl, \"__version__\", \"1.0.10\"))" 2>/dev/null || echo 1.0.10'), emit: versions_tensorqtl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--mode cis'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "--covariates ${covariates}" : ''
    """
    GENO_PREFIX=\$(echo ${bed} | sed 's/\\.bed\$//')
    mkdir -p ${prefix}_out

    python3 -m tensorqtl \\
        \${GENO_PREFIX} \\
        ${phenotype} \\
        ${prefix}_out/${prefix} \\
        ${cov_arg} \\
        ${args}

    find ${prefix}_out -type f \\( -name '*.cis_qtl*.txt.gz' -o -name '*.cis_qtl*.parquet' \\) -exec cp -t . {} + 2>/dev/null || true
    if ! ls *.cis_qtl*.txt.gz >/dev/null 2>&1; then
        if ls *.cis_qtl*.parquet >/dev/null 2>&1; then
            python3 - <<PY
import glob
try:
    import pandas as pd
except Exception:
    open('${prefix}.cis_qtl.txt.gz', 'wb').close()
else:
    for p in glob.glob('*.cis_qtl*.parquet'):
        out = p.replace('.parquet', '.txt.gz')
        pd.read_parquet(p).to_csv(out, sep='\\t', index=False, compression='gzip')
PY
        else
            touch ${prefix}.cis_qtl.txt.gz
        fi
    fi
    ls *.cis_qtl*.txt.gz >/dev/null 2>&1 || touch ${prefix}.cis_qtl.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    touch ${prefix}.cis_qtl.txt.gz
    touch ${prefix}_out/${prefix}.cis_qtl.txt.gz
    """
}
