process EMMAX_ASSOC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/emmax:0.0.20100307' :
        'quay.io/bioinfortools/emmax:0.0.20100307' }"

    input:
    tuple val(meta), path(tped), path(tfam)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(kinship)
    tuple val(meta4), path(covariates)

    output:
    tuple val(meta), path("*.ps"), emit: assoc
    tuple val(meta), path("*.reml"), emit: reml, optional: true
    tuple val("${task.process}"), val('emmax'), eval('echo 0.0.20100307'), emit: versions_emmax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-v -d 10'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "-c ${covariates}" : ''
    """
    TPED_PREFIX=\$(echo ${tped} | sed 's/\\.tped\$//')
    if [ ! -f \${TPED_PREFIX}.tfam ]; then
        ln -sf ${tfam} \${TPED_PREFIX}.tfam
    fi

    emmax \\
        ${args} \\
        -t \${TPED_PREFIX} \\
        -p ${phenotype} \\
        -k ${kinship} \\
        ${cov_arg} \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ps
    touch ${prefix}.reml
    """
}
