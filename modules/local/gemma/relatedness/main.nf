process GEMMA_RELATEDNESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gemma:0.98.5--h38cc83e_1' :
        'quay.io/biocontainers/gemma:0.98.5--h38cc83e_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.cXX.txt"), emit: relatedness
    tuple val(meta), path("*.log.txt"), emit: log, optional: true
    tuple val("${task.process}"), val('gemma'), eval('gemma -h 2>&1 | sed -n "s/.*version \\([0-9.]*\\).*/\\1/p" | head -1'), emit: versions_gemma, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-gk 1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bfile = bed.baseName.replaceFirst(/\\.bed$/, '')
    // bed path may be sample.bed -> baseName sample; use file without extension via getBaseName in Nextflow
    """
    # PLINK bfile prefix: strip .bed suffix from filename
    BFILE=\$(echo ${bed} | sed 's/\\.bed\$//')

    gemma \\
        -bfile \${BFILE} \\
        ${args} \\
        -o ${prefix}

    # gemma writes to ./output/ by default
    if [ -d output ]; then
        mv output/${prefix}.* . 2>/dev/null || true
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cXX.txt
    touch ${prefix}.log.txt
    """
}
