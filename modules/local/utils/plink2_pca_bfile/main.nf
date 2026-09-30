process PLINK2_PCA_BFILE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0':
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.eigenvec"), emit: eigenvec
    tuple val(meta), path("*.eigenval"), emit: eigenval
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def npcs = task.ext.npcs ?: params.genotype_qc_pca_n ?: 10
    def bfile = bed.baseName
    """
    plink2 \\
        --bfile ${bfile} \\
        --pca ${npcs} \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eigenvec
    touch ${prefix}.eigenval
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: 2.00
    END_VERSIONS
    """
}
