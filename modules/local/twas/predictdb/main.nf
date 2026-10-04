process TWAS_PREDICTDB {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.weights.tsv"), emit: weights
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_py = "${moduleDir}/templates/predictdb_extract.py"
    """
    cp '${template_py}' predictdb_extract.py
    python3 predictdb_extract.py \\
        --db ${db} \\
        --out ${prefix}.weights.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "gene\\tsnp\\tweight\\tgene_id" > ${prefix}.weights.tsv
    echo -e "GENE_A\\trs1\\t0.4\\tENSG00000001" >> ${prefix}.weights.tsv
    echo -e "GENE_A\\trs2\\t-0.2\\tENSG00000001" >> ${prefix}.weights.tsv
    """
}
