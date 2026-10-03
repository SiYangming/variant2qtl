process FINEMAP_EXTRA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(sumstats)

    output:
    tuple val(meta), path("*.pip.tsv"), emit: pip
    tuple val(meta), path("*.cs.tsv"), emit: credible_sets
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--method finemap'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_py = "${moduleDir}/templates/finemap_extra.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' finemap_extra.py

    python3 finemap_extra.py \\
        --sumstats ${sumstats} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.pip.tsv .
    cp ${prefix}_out/${prefix}.cs.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "variant_id\\tz\\tbf\\tpip\\tmethod" > ${prefix}.pip.tsv
    echo -e "snp1\\t5.0\\t10\\t0.8\\tfinemap" >> ${prefix}.pip.tsv
    echo -e "credible_set\\tvariant_id\\tpip" > ${prefix}.cs.tsv
    echo -e "CS1\\tsnp1\\t0.8" >> ${prefix}.cs.tsv
    cp ${prefix}.pip.tsv ${prefix}_out/
    cp ${prefix}.cs.tsv ${prefix}_out/
    """
}
