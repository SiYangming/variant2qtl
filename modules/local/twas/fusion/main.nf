process TWAS_FUSION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(weights)
    tuple val(meta2), path(gwas)
    tuple val(meta3), path(ld)

    output:
    tuple val(meta), path("*.twas.tsv"), emit: twas
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ld_arg = ld ? "--ld ${ld}" : ''
    def template_py = "${moduleDir}/templates/twas_fusion.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' twas_fusion.py

    python3 twas_fusion.py \\
        --weights ${weights} \\
        --gwas ${gwas} \\
        ${ld_arg} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.twas.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "gene\\tn_snps\\tz\\tp\\tld" > ${prefix}.twas.tsv
    echo -e "GENE_A\\t2\\t3.1\\t0.002\\tno" >> ${prefix}.twas.tsv
    cp ${prefix}.twas.tsv ${prefix}_out/
    """
}
