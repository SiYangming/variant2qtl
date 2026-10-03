process COLOC_ABF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.1' :
        'quay.io/biocontainers/r-base:4.3.1' }"

    input:
    tuple val(meta), path(qtl_sumstats)
    tuple val(meta2), path(gwas_sumstats)

    output:
    tuple val(meta), path("*.summary.tsv"), emit: summary
    tuple val(meta), path("*.abf.tsv"), emit: abf
    tuple val(meta), path("*.coloc.log"), emit: log, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('coloc'), eval('Rscript -e \'cat(as.character(packageVersion("coloc")))\' 2>/dev/null || echo 5.2.3'), emit: versions_coloc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_r = "${moduleDir}/templates/coloc_abf.R"
    """
    mkdir -p ${prefix}_out
    cp '${template_r}' coloc_abf.R

    Rscript coloc_abf.R \\
        --qtl ${qtl_sumstats} \\
        --gwas ${gwas_sumstats} \\
        --out ${prefix}_out/${prefix} \\
        ${args}

    cp ${prefix}_out/${prefix}.summary.tsv .
    cp ${prefix}_out/${prefix}.abf.tsv .
    cp ${prefix}_out/${prefix}.coloc.log . 2>/dev/null || touch ${prefix}.coloc.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "nsnps\\tPP.H0.abf\\tPP.H1.abf\\tPP.H2.abf\\tPP.H3.abf\\tPP.H4.abf" > ${prefix}.summary.tsv
    echo -e "10\\t0.01\\t0.02\\t0.02\\t0.10\\t0.85" >> ${prefix}.summary.tsv
    echo -e "snp\\tSNP.PP.H4" > ${prefix}.abf.tsv
    echo -e "snp1\\t0.4" >> ${prefix}.abf.tsv
    echo "stub coloc" > ${prefix}.coloc.log
    cp ${prefix}.summary.tsv ${prefix}_out/
    cp ${prefix}.abf.tsv ${prefix}_out/
    cp ${prefix}.coloc.log ${prefix}_out/
    """
}
