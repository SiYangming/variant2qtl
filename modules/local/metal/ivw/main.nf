process METAL_IVW {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(cohorts)

    output:
    tuple val(meta), path("*.metal.tsv"), emit: metal
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_py = "${moduleDir}/templates/metal_ivw.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' metal_ivw.py

    python3 metal_ivw.py \\
        --cohorts ${cohorts} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.metal.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "snp\\tbeta\\tse\\tz\\tp\\tn_studies\\tq_het\\tp_het" > ${prefix}.metal.tsv
    echo -e "rs1\\t0.44\\t0.077\\t5.7\\t1e-8\\t2\\t0.4\\t0.5" >> ${prefix}.metal.tsv
    cp ${prefix}.metal.tsv ${prefix}_out/
    """
}
