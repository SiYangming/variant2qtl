process OMIGA_CIS {
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
    tuple val(meta), path("*.cis_qtl*.txt.gz"), emit: cis_qtl, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('omiga'), eval('omiga --version 2>&1 | head -1 | sed "s/[^0-9.]//g" || echo 1.8.17'), emit: versions_omiga, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--qtl-map-model a+A --cis-window 1000000 --verbose'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "--covariates ${covariates}" : ''
    """
    GENO_PREFIX=\$(echo ${bed} | sed 's/\\.bed\$//')
    mkdir -p ${prefix}_out

    omiga \\
        --mode cis \\
        --genotype \${GENO_PREFIX} \\
        --phenotype ${phenotype} \\
        ${cov_arg} \\
        --prefix ${prefix} \\
        --output-dir ${prefix}_out \\
        --threads ${task.cpus} \\
        ${args}

    find ${prefix}_out -type f -name '*.cis_qtl*.txt.gz' -exec cp -t . {} + 2>/dev/null || true
    ls *.cis_qtl*.txt.gz >/dev/null 2>&1 || touch ${prefix}.cis_qtl.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    touch ${prefix}.cis_qtl.txt.gz
    touch ${prefix}_out/${prefix}.cis_qtl.txt.gz
    """
}
