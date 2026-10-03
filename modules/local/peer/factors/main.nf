process PEER_FACTORS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/peer:1.3--h503566f_1' :
        'quay.io/biocontainers/peer:1.3--h503566f_1' }"

    input:
    tuple val(meta), path(phenotype)
    tuple val(meta2), path(covariates)

    output:
    tuple val(meta), path("*.peer.tensor.tsv"), emit: cov_tensorqtl
    tuple val(meta), path("*.peer.omiga.tsv"), emit: cov_omiga
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('peer'), eval('Rscript -e \'cat(as.character(packageVersion("peer")))\' 2>/dev/null || echo 1.3'), emit: versions_peer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--nk 10'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "--covariates ${covariates}" : ''
    def template_r = "${moduleDir}/templates/peer_factors.R"
    """
    mkdir -p ${prefix}_out
    cp '${template_r}' peer_factors.R

    Rscript peer_factors.R \\
        --phenotype ${phenotype} \\
        ${cov_arg} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.peer.tensor.tsv .
    cp ${prefix}_out/${prefix}.peer.omiga.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "id\\tPEER1\\tPEER2" > ${prefix}.peer.tensor.tsv
    echo -e "s1\\t0.1\\t-0.2" >> ${prefix}.peer.tensor.tsv
    echo -e "s2\\t-0.1\\t0.2" >> ${prefix}.peer.tensor.tsv
    echo -e "id\\ts1\\ts2" > ${prefix}.peer.omiga.tsv
    echo -e "PEER1\\t0.1\\t-0.1" >> ${prefix}.peer.omiga.tsv
    echo -e "PEER2\\t-0.2\\t0.2" >> ${prefix}.peer.omiga.tsv
    cp ${prefix}.peer.tensor.tsv ${prefix}_out/
    cp ${prefix}.peer.omiga.tsv ${prefix}_out/
    """
}
