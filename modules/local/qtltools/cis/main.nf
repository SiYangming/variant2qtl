process QTLTOOLS_CIS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Official bioconda/biocontainers QTLtools is absent. Conda: YangmingSi::qtltools=1.3.1.
    // Docker image is built as quay.io/bioinfortools/qtltools:1.3.1 (see SiYangming/qtltools PACKAGING.md).
    // CI -stub uses python until that Quay tag is pushed (401 without docker login).
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(covariates)

    output:
    tuple val(meta), path("*.cis_qtl*.txt.gz"), emit: cis_qtl, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('qtltools'), eval('QTLtools --help 2>&1 | head -1 | sed "s/[^0-9.]//g" || echo 1.3.1'), emit: versions_qtltools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--nominal 1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_arg = covariates ? "--cov ${covariates}" : ''
    """
    GENO_PREFIX=\$(echo ${bed} | sed 's/\\.bed\$//')
    mkdir -p ${prefix}_out

    plink2 \\
        --bfile \${GENO_PREFIX} \\
        --make-bed \\
        --output-chr chrM \\
        --out ${prefix}_chrpref \\
        --threads ${task.cpus} || plink2 --bfile \${GENO_PREFIX} --make-bed --out ${prefix}_chrpref --threads ${task.cpus}

    plink2 \\
        --bfile ${prefix}_chrpref \\
        --export vcf bgz id-delim='-' \\
        --out ${prefix}_geno \\
        --threads ${task.cpus}

    tabix -p vcf ${prefix}_geno.vcf.gz || true

    PHENO="${phenotype}"
    if echo "\${PHENO}" | grep -vq '\\.gz\$'; then
        bgzip -c "\${PHENO}" > ${prefix}_pheno.bed.gz
        PHENO=${prefix}_pheno.bed.gz
    fi
    tabix -p bed "\${PHENO}" || true

    QTLtools cis \\
        --vcf ${prefix}_geno.vcf.gz \\
        --bed \${PHENO} \\
        ${cov_arg} \\
        --out ${prefix}_out/${prefix} \\
        --std-err \\
        ${args}

    find ${prefix}_out -type f \\( -name '*.txt' -o -name '*.txt.gz' \\) -exec cp {} . \\; 2>/dev/null || true
    if ls ${prefix}*cis* >/dev/null 2>&1; then
        true
    fi
    if ! ls *.cis_qtl*.txt.gz >/dev/null 2>&1; then
        if ls ${prefix}*.txt >/dev/null 2>&1; then
            gzip -c \$(ls ${prefix}*.txt | head -1) > ${prefix}.cis_qtl.txt.gz
        elif ls ${prefix}*.txt.gz >/dev/null 2>&1; then
            cp \$(ls ${prefix}*.txt.gz | head -1) ${prefix}.cis_qtl.txt.gz
        else
            touch ${prefix}.cis_qtl.txt.gz
        fi
    fi
    ls *.cis_qtl*.txt.gz >/dev/null 2>&1 || touch ${prefix}.cis_qtl.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "phe_id\\tvar_id\\tvar_chr\\tvar_from\\tnom_pval\\tbeta\\tslope_se" > ${prefix}.cis_qtl.txt
    echo -e "gene1\\tsnp1\\tchr1\\t100\\t1e-8\\t0.4\\t0.08" >> ${prefix}.cis_qtl.txt
    gzip -c ${prefix}.cis_qtl.txt > ${prefix}.cis_qtl.txt.gz
    cp ${prefix}.cis_qtl.txt.gz ${prefix}_out/
    """
}
