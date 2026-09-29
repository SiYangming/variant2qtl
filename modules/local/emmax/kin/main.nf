process EMMAX_KIN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/emmax:0.0.20100307' :
        'quay.io/bioinfortools/emmax:0.0.20100307' }"

    input:
    tuple val(meta), path(tped), path(tfam)

    output:
    tuple val(meta), path("*.kinf"), emit: kinship
    tuple val("${task.process}"), val('emmax'), eval('echo 0.0.20100307'), emit: versions_emmax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Default IBS kinship; override with task.ext.args (e.g. BN: -v -h -d 10)
    def args = task.ext.args ?: '-v -h -s -d 10'
    def prefix = task.ext.prefix ?: tped.baseName.replaceFirst(/\\.tped$/, '')
    """
    # Ensure tped/tfam share the expected prefix for emmax-kin
    TPED_PREFIX=\$(echo ${tped} | sed 's/\\.tped\$//')
    if [ ! -f \${TPED_PREFIX}.tfam ]; then
        ln -sf ${tfam} \${TPED_PREFIX}.tfam
    fi

    emmax-kin \\
        ${args} \\
        \${TPED_PREFIX}

    # Normalize output name if prefix override requested
    if [ "${prefix}" != "\${TPED_PREFIX}" ]; then
        for f in \${TPED_PREFIX}*.kinf; do
            [ -f "\$f" ] || continue
            mv "\$f" "${prefix}.\${f#\${TPED_PREFIX}.}"
        done
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hIBS.kinf
    """
}
