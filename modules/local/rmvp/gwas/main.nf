process RMVP_GWAS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/rmvp:1.4.6' :
        'quay.io/bioinfortools/rmvp:1.4.6' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(phenotype)

    output:
    tuple val(meta), path("*.csv"), emit: results, optional: true
    tuple val(meta), path("*.png"), emit: plots, optional: true
    tuple val(meta), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('r-rmvp'), eval('Rscript -e \'cat(as.character(packageVersion("rMVP")))\''), emit: versions_rmvp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--method FarmCPU'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_r = "${moduleDir}/templates/rmvp_run.R"
    """
    cp '${template_r}' rmvp_run.R

    Rscript rmvp_run.R \\
        --vcf ${vcf} \\
        --phenotype ${phenotype} \\
        --out ${prefix} \\
        --ncpus ${task.cpus} \\
        ${args} \\
        > ${prefix}.log 2>&1

    # Ensure channel-friendly placeholders if MVP naming differs
    ls *.csv >/dev/null 2>&1 || touch ${prefix}.FarmCPU.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.FarmCPU.csv
    touch ${prefix}.log
    """
}
