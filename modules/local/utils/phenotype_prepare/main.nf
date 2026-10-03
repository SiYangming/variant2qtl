process PHENOTYPE_PREPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(matrix)
    tuple val(meta2), path(gene_map)
    tuple val(meta3), path(samples)

    output:
    tuple val(meta), path("*.phenotype.bed.gz"), emit: bed
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--transform invnorm --max-missing 0.2'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gene_arg = gene_map ? "--gene-map ${gene_map}" : ''
    def samples_arg = samples ? "--samples ${samples}" : ''
    def template_py = "${moduleDir}/templates/phenotype_prepare.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' phenotype_prepare.py

    python3 phenotype_prepare.py \\
        --matrix ${matrix} \\
        ${gene_arg} \\
        ${samples_arg} \\
        --out ${prefix}_out/${prefix}.phenotype.bed \\
        ${args}

    gzip -c ${prefix}_out/${prefix}.phenotype.bed > ${prefix}.phenotype.bed.gz
    cp ${prefix}.phenotype.bed.gz ${prefix}_out/
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "#chr\\tstart\\tend\\tpheno_id\\ts1\\ts2" > ${prefix}.phenotype.bed
    echo -e "1\\t100\\t200\\tGENE_A\\t0.1\\t-0.1" >> ${prefix}.phenotype.bed
    gzip -c ${prefix}.phenotype.bed > ${prefix}.phenotype.bed.gz
    cp ${prefix}.phenotype.bed.gz ${prefix}_out/
    """
}
