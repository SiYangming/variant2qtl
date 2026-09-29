process GEMMA_LMM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gemma:0.98.5--h38cc83e_1' :
        'quay.io/biocontainers/gemma:0.98.5--h38cc83e_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(relatedness)
    tuple val(meta4), path(covariates)

    output:
    tuple val(meta), path("*.assoc.txt"), emit: assoc
    tuple val(meta), path("*.log.txt"), emit: log, optional: true
    tuple val("${task.process}"), val('gemma'), eval('gemma -h 2>&1 | sed -n "s/.*version \\([0-9.]*\\).*/\\1/p" | head -1'), emit: versions_gemma, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-lmm 1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "-c ${covariates}" : ''
    """
    BFILE=\$(echo ${bed} | sed 's/\\.bed\$//')

    gemma \\
        -bfile \${BFILE} \\
        -p ${phenotype} \\
        -k ${relatedness} \\
        ${cov_arg} \\
        ${args} \\
        -o ${prefix}

    if [ -d output ]; then
        mv output/${prefix}.* . 2>/dev/null || true
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.assoc.txt
    touch ${prefix}.log.txt
    """
}
