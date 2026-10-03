process LEAFCUTTER_PREPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(counts)
    tuple val(meta2), path(genes)

    output:
    tuple val(meta), path("*.phenotype.bed.gz"), emit: bed
    tuple val(meta), path("*.phenotype_group.tsv"), emit: phenotype_group
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--transform invnorm --min-ratio 0.01'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gene_arg = genes ? "--genes ${genes}" : ''
    def template_py = "${moduleDir}/templates/leafcutter_prepare.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' leafcutter_prepare.py

    python3 leafcutter_prepare.py \\
        --counts ${counts} \\
        ${gene_arg} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    gzip -c ${prefix}_out/${prefix}.phenotype.bed > ${prefix}.phenotype.bed.gz
    cp ${prefix}_out/${prefix}.phenotype_group.tsv ${prefix}.phenotype_group.tsv
    cp ${prefix}.phenotype.bed.gz ${prefix}_out/
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "#chr\\tstart\\tend\\tpheno_id\\ts1\\ts2" > ${prefix}.phenotype.bed
    echo -e "chr1\\t100\\t200\\tchr1:100:200:clu_1\\t0.1\\t-0.1" >> ${prefix}.phenotype.bed
    gzip -c ${prefix}.phenotype.bed > ${prefix}.phenotype.bed.gz
    echo -e "phenotype_id\\tgroup_id" > ${prefix}.phenotype_group.tsv
    echo -e "chr1:100:200:clu_1\\tGENE_A" >> ${prefix}.phenotype_group.tsv
    cp ${prefix}.phenotype.bed.gz ${prefix}_out/
    cp ${prefix}.phenotype_group.tsv ${prefix}_out/
    """
}
