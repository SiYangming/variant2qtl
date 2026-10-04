process LDSC_H2 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta2), path(annot)

    output:
    tuple val(meta), path("*.h2.tsv"), emit: h2
    tuple val(meta), path("*.part.tsv"), emit: partitioned, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def annot_arg = annot ? "--annot ${annot}" : ''
    def template_py = "${moduleDir}/templates/ldsc_h2.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' ldsc_h2.py

    python3 ldsc_h2.py \\
        --sumstats ${sumstats} \\
        ${annot_arg} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.h2.tsv .
    cp ${prefix}_out/${prefix}.part.tsv . 2>/dev/null || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "trait\\tnsnp\\tn_bar\\th2\\th2_se\\tslope" > ${prefix}.h2.tsv
    echo -e "trait1\\t10\\t10000\\t0.2\\t0.05\\t0.001" >> ${prefix}.h2.tsv
    echo -e "annotation\\tnsnp\\th2\\tenrichment" > ${prefix}.part.tsv
    echo -e "base\\t10\\t0.2\\t1" >> ${prefix}.part.tsv
    cp ${prefix}.h2.tsv ${prefix}_out/
    cp ${prefix}.part.tsv ${prefix}_out/
    """
}
