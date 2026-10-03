process HYPRCOLOC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(traits)

    output:
    tuple val(meta), path("*.clusters.tsv"), emit: clusters
    tuple val(meta), path("*.pairs.tsv"), emit: pairs
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--pp4-threshold 0.5'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_py = "${moduleDir}/templates/hyprcoloc_cluster.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' hyprcoloc_cluster.py

    python3 hyprcoloc_cluster.py \\
        --traits ${traits} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.clusters.tsv .
    cp ${prefix}_out/${prefix}.pairs.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "cluster\\ttraits\\tn_traits" > ${prefix}.clusters.tsv
    echo -e "C1\\tQTL,GWAS\\t2" >> ${prefix}.clusters.tsv
    echo -e "trait1\\ttrait2\\tPP.H4\\tcandidate_snp" > ${prefix}.pairs.tsv
    echo -e "QTL\\tGWAS\\t0.8\\trs5" >> ${prefix}.pairs.tsv
    cp ${prefix}.clusters.tsv ${prefix}_out/
    cp ${prefix}.pairs.tsv ${prefix}_out/
    """
}
