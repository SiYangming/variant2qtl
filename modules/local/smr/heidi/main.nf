process SMR_HEIDI {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(qtl_sumstats)
    tuple val(meta2), path(gwas_sumstats)

    output:
    tuple val(meta), path("*.smr.tsv"), emit: smr
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def template_py = "${moduleDir}/templates/smr_heidi.py"
    """
    mkdir -p ${prefix}_out
    cp '${template_py}' smr_heidi.py

    python3 smr_heidi.py \\
        --qtl ${qtl_sumstats} \\
        --gwas ${gwas_sumstats} \\
        --out ${prefix}_out/${prefix}

    cp ${prefix}_out/${prefix}.smr.tsv .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "gene\\ttop_snp\\tb_smr\\tse_smr\\tt_smr\\tp_smr\\theidi_Q\\tp_heidi\\tn_snps" > ${prefix}.smr.tsv
    echo -e "GENE_A\\trs1\\t0.5\\t0.1\\t25\\t5e-7\\t1.2\\t0.3\\t8" >> ${prefix}.smr.tsv
    cp ${prefix}.smr.tsv ${prefix}_out/
    """
}
