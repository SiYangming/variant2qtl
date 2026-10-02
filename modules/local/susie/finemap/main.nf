process SUSIE_FINEMAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.1' :
        'quay.io/biocontainers/r-base:4.3.1' }"

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta2), path(ld)

    output:
    tuple val(meta), path("*.pip.tsv"), emit: pip
    tuple val(meta), path("*.cs.tsv"), emit: credible_sets
    tuple val(meta), path("*.susie.log"), emit: log, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('susieR'), eval('Rscript -e \'cat(as.character(packageVersion("susieR")))\' 2>/dev/null || echo 0.14.2'), emit: versions_susier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ld_arg = ld ? "--ld ${ld}" : ''
    def template_r = "${moduleDir}/templates/susie_finemap.R"
    """
    mkdir -p ${prefix}_out
    cp '${template_r}' susie_finemap.R

    Rscript susie_finemap.R \\
        --sumstats ${sumstats} \\
        ${ld_arg} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.pip.tsv .
    cp ${prefix}_out/${prefix}.cs.tsv .
    cp ${prefix}_out/${prefix}.susie.log . 2>/dev/null || touch ${prefix}.susie.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "variant_id\\tz\\tpip" > ${prefix}.pip.tsv
    echo -e "snp1\\t5.0\\t0.9" >> ${prefix}.pip.tsv
    echo -e "credible_set\\tvariant_id\\tpip" > ${prefix}.cs.tsv
    echo -e "L1\\tsnp1\\t0.9" >> ${prefix}.cs.tsv
    echo "stub susieR" > ${prefix}.susie.log
    cp ${prefix}.pip.tsv ${prefix}_out/
    cp ${prefix}.cs.tsv ${prefix}_out/
    cp ${prefix}.susie.log ${prefix}_out/
    """
}
