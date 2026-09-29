process PHENOCOVAR_ADAPT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(phenotype)
    tuple val(meta2), path(covariates)

    output:
    tuple val(meta), path("*.gemma.pheno.txt")  , emit: gemma_pheno
    tuple val(meta), path("*.gemma.covar.txt")  , emit: gemma_covar , optional: true
    tuple val(meta), path("*.emmax.pheno.txt")  , emit: emmax_pheno
    tuple val(meta), path("*.emmax.covar.txt")  , emit: emmax_covar , optional: true
    tuple val(meta), path("*.tassel.pheno.txt") , emit: tassel_pheno
    tuple val(meta), path("*.rmvp.pheno.txt")   , emit: rmvp_pheno
    tuple val(meta), path("*.omiga.pheno.txt")  , emit: omiga_pheno
    tuple val(meta), path("*.omiga.covar.txt")  , emit: omiga_covar , optional: true
    tuple val("${task.process}"), val('gawk'), eval('gawk --version | head -n1 | sed "s/GNU Awk //; s/,.*//"'), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def has_cov = covariates ? true : false
    """
    # Phenotype: FID IID trait... (plink header style)
    # Covariates (optional): FID IID cov1... (plink header style)

    # GEMMA: no header; one phenotype column (sample order = input order)
    # GEMMA: covariates = numeric columns only, no header (tool adds intercept)
    awk 'NR==1 { next } { print \$3 }' ${phenotype} > ${prefix}.gemma.pheno.txt

    # EMMAX: no header; FAMID INDID pheno
    awk 'NR==1 { next } { print \$1, \$2, \$3 }' ${phenotype} > ${prefix}.emmax.pheno.txt

    # TASSEL: Taxa + trait (header); Taxa = IID
    TRAIT=\$(awk 'NR==1 { print \$3; exit }' ${phenotype})
    {
        echo -e "Taxa\t\${TRAIT}"
        awk 'NR==1 { next } { print \$2 "\t" \$3 }' ${phenotype}
    } > ${prefix}.tassel.pheno.txt

    # rMVP: header Taxa + trait; Taxa = IID (matches VCF sample IDs)
    {
        echo -e "Taxa\t\${TRAIT}"
        awk 'NR==1 { next } { print \$2 "\t" \$3 }' ${phenotype}
    } > ${prefix}.rmvp.pheno.txt

    # OmiGA GWAS: keep FID IID + trait with header
    awk 'NR==1 { print \$1 "\t" \$2 "\t" \$3; next } { print \$1 "\t" \$2 "\t" \$3 }' ${phenotype} > ${prefix}.omiga.pheno.txt

    if ${has_cov}; then
        # GEMMA covar: drop FID/IID, no header
        awk 'NR==1 { next } { for (i=3; i<=NF; i++) printf "%s%s", \$i, (i<NF ? OFS : ORS) }' ${covariates} > ${prefix}.gemma.covar.txt

        # EMMAX covar: FID IID intercept cov... (intercept in column 3)
        awk 'NR==1 { next } { printf "%s %s 1", \$1, \$2; for (i=3; i<=NF; i++) printf " %s", \$i; printf ORS }' ${covariates} > ${prefix}.emmax.covar.txt

        # OmiGA covar: keep header FID IID + covariates
        awk '{
            printf "%s\t%s", \$1, \$2
            for (i=3; i<=NF; i++) printf "\t%s", \$i
            printf ORS
        }' ${covariates} > ${prefix}.omiga.covar.txt
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def has_cov = covariates ? true : false
    """
    echo -e "1.0\\n2.0" > ${prefix}.gemma.pheno.txt
    echo -e "fam1 ind1 1.0\\nfam2 ind2 2.0" > ${prefix}.emmax.pheno.txt
    echo -e "Taxa\tQuantitativeTrait\\nind1\t1.0\\nind2\t2.0" > ${prefix}.tassel.pheno.txt
    echo -e "Taxa\tQuantitativeTrait\\nind1\t1.0\\nind2\t2.0" > ${prefix}.rmvp.pheno.txt
    echo -e "FID\tIID\tQuantitativeTrait\\nfam1\tind1\t1.0\\nfam2\tind2\t2.0" > ${prefix}.omiga.pheno.txt
    if ${has_cov}; then
        echo -e "0 1\\n1 0" > ${prefix}.gemma.covar.txt
        echo -e "fam1 ind1 1 0 1\\nfam2 ind2 1 1 0" > ${prefix}.emmax.covar.txt
        echo -e "FID\tIID\tSex\tAge\\nfam1\tind1\t1\t40\\nfam2\tind2\t2\t41" > ${prefix}.omiga.covar.txt
    fi
    """
}
