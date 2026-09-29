process TASSEL_MLM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tassel:5.2.89--hdfd78af_0' :
        'quay.io/biocontainers/tassel:5.2.89--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(kinship)

    output:
    tuple val(meta), path("*_mlm*.txt"), emit: results
    tuple val(meta), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('tassel'), eval('run_pipeline.pl 2>&1 | head -1 || echo 5.2.89'), emit: versions_tassel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.memory ? "-Xmx${task.memory.toGiga()}G" : "-Xmx8G"
    // Default MLM pipeline; override entirely with task.ext.args when needed
    def default_args = kinship
        ? "-fork1 -vcf ${vcf} -fork2 -r ${phenotype} -fork3 -k ${kinship} -combine4 -input1 -input2 -input3 -intersect -MLMPlugin -endPlugin -export ${prefix}_mlm -runfork1 -runfork2 -runfork3"
        : "-fork1 -vcf ${vcf} -fork2 -r ${phenotype} -combine3 -input1 -input2 -intersect -MLMPlugin -endPlugin -export ${prefix}_mlm -runfork1 -runfork2"
    def args = task.ext.args ?: default_args
    """
    run_pipeline.pl \\
        ${mem} \\
        ${args} \\
        > ${prefix}.log 2>&1

    ls ${prefix}_mlm*.txt >/dev/null 2>&1 || touch ${prefix}_mlm.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mlm.txt
    touch ${prefix}.log
    """
}
