process TORUS_ENRICH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(annot)

    output:
    tuple val(meta), path("*.prior.tsv"), emit: prior
    tuple val(meta), path("*.enrich.tsv"), emit: enrich
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_py = "${moduleDir}/templates/torus_enrich.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' torus_enrich.py

    python3 torus_enrich.py \\
        --annot ${annot} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.prior.tsv .
    cp ${prefix}_out/${prefix}.enrich.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "snp\\tz\\tprior" > ${prefix}.prior.tsv
    echo -e "rs1\\t4.5\\t0.4" >> ${prefix}.prior.tsv
    echo -e "term\\tlog_odds" > ${prefix}.enrich.tsv
    echo -e "coding\\t1.2" >> ${prefix}.enrich.tsv
    cp ${prefix}.prior.tsv ${prefix}_out/
    cp ${prefix}.enrich.tsv ${prefix}_out/
    """
}
