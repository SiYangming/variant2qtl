process OMIGA_GWAS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/omiga:1.8.17' :
        'quay.io/bioinfortools/omiga:1.8.17' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(covariates)

    output:
    tuple val(meta), path("*.gwas*"), emit: gwas, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('omiga'), eval('omiga --version 2>&1 | head -1 | sed "s/[^0-9.]//g" || echo 1.8.17'), emit: versions_omiga, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--qtl-map-model a+A --verbose'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "--covariates ${covariates}" : ''
    """
    GENO_PREFIX=\$(echo ${bed} | sed 's/\\.bed\$//')
    mkdir -p ${prefix}_out

    omiga \\
        --mode gwas \\
        --genotype \${GENO_PREFIX} \\
        --phenotype ${phenotype} \\
        ${cov_arg} \\
        --prefix ${prefix} \\
        --output-dir ${prefix}_out \\
        --threads ${task.cpus} \\
        ${args}

    find ${prefix}_out -type f -name '*gwas*' -exec cp -t . {} + 2>/dev/null || true
    ls *gwas* >/dev/null 2>&1 || touch ${prefix}.gwas.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    touch ${prefix}.gwas.txt
    touch ${prefix}_out/${prefix}.gwas.txt
    """
}
