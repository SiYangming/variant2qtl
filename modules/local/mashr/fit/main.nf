process MASHR_FIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.1' :
        'quay.io/biocontainers/r-base:4.3.1' }"

    input:
    tuple val(meta), path(effects)

    output:
    tuple val(meta), path("*.mash.tsv"), emit: mash
    tuple val(meta), path("*.lfsr.tsv"), emit: lfsr
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('mashr'), eval('Rscript -e \'cat(as.character(packageVersion("mashr")))\' 2>/dev/null || echo 0.2.79'), emit: versions_mashr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_r = "${moduleDir}/templates/mashr_fit.R"
    """
    mkdir -p ${prefix}_out
    cp '${template_r}' mashr_fit.R

    Rscript mashr_fit.R \\
        --effects ${effects} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.mash.tsv .
    cp ${prefix}_out/${prefix}.lfsr.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "snp\\tcondition\\tposterior_mean\\tposterior_sd\\tlfsr" > ${prefix}.mash.tsv
    echo -e "rs1\\tliver\\t0.6\\t0.08\\t0.01" >> ${prefix}.mash.tsv
    echo -e "\\tliver\\tbrain" > ${prefix}.lfsr.tsv
    echo -e "rs1\\t0.01\\t0.02" >> ${prefix}.lfsr.tsv
    cp ${prefix}.mash.tsv ${prefix}_out/
    cp ${prefix}.lfsr.tsv ${prefix}_out/
    """
}
