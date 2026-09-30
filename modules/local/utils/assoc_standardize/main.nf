process ASSOC_STANDARDIZE {
    tag "${meta.id}-${tool}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(assoc), val(tool)

    output:
    tuple val(meta), path("*.std.tsv"), emit: standardized
    tuple val("${task.process}"), val('gawk'), eval('gawk --version | head -n1 | sed "s/GNU Awk //; s/,.*//"'), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${tool}"
    // Use awk string matches (not /regex/) so Groovy does not parse slashy strings inside """
    """
    TOOL="${tool}"
    ASSOC="${assoc}"
    OUT="${prefix}.std.tsv"

    echo -e "CHR\tPOS\tSNP\tP" > "\${OUT}"

    case "\${TOOL}" in
        gemma)
            # GEMMA assoc.txt header: chr rs ps ... p_wald|p_lrt|p_score
            awk -v OFS='\\t' '
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        h=\$i; gsub("\\r", "", h)
                        if (h=="chr") c=i
                        if (h=="rs") s=i
                        if (h=="ps") p=i
                        if (h=="p_lrt" || h=="p_wald" || h=="p_score") pv=i
                    }
                    next
                }
                c && s && p && pv { print \$c, \$p, \$s, \$pv }
            ' "\${ASSOC}" >> "\${OUT}"
            ;;
        emmax)
            # EMMAX .ps: SNP beta p (no header; no chr/pos)
            awk -v OFS='\\t' 'NF>=3 { print "NA", "NA", \$1, \$3 }' "\${ASSOC}" >> "\${OUT}"
            ;;
        tassel)
            # TASSEL MLM export: Marker / Chr / Pos / p (name variants)
            awk -v OFS='\\t' '
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        h=tolower(\$i); gsub("\\r", "", h)
                        if (h=="marker" || h=="snp") s=i
                        if (h=="chr" || h=="chromosome") c=i
                        if (h=="pos" || h=="position" || h=="site") p=i
                        if (h=="p" || h=="p-value" || h=="pvalue" || h=="marker_p") pv=i
                    }
                    next
                }
                {
                    chr = c ? \$c : "NA"
                    pos = p ? \$p : "NA"
                    snp = s ? \$s : "NA"
                    pval = pv ? \$pv : "NA"
                    if (snp != "NA" || pval != "NA") print chr, pos, snp, pval
                }
            ' "\${ASSOC}" >> "\${OUT}"
            ;;
        rmvp)
            # rMVP CSV: SNP, Chrs, Pos, *.P.value / pvalue
            awk -v FS=',' -v OFS='\\t' '
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        h=\$i; gsub("\\r", "", h); gsub("\\"", "", h); hl=tolower(h)
                        if (hl=="snp") s=i
                        if (hl=="chrs" || hl=="chr" || hl=="chrom") c=i
                        if (hl=="pos" || hl=="position") p=i
                        if (hl ~ "p\\\\.value\$" || hl=="pvalue" || hl=="p" || hl ~ "\\\\.p\\\\.value\$") pv=i
                    }
                    next
                }
                s && pv {
                    chr = c ? \$c : "NA"
                    pos = p ? \$p : "NA"
                    gsub("\\"", "", \$s); gsub("\\"", "", chr); gsub("\\"", "", pos); gsub("\\"", "", \$pv)
                    print chr, pos, \$s, \$pv
                }
            ' "\${ASSOC}" >> "\${OUT}"
            ;;
        omiga)
            # OmiGA GWAS: best-effort common colnames
            awk -v OFS='\\t' '
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        h=\$i; gsub("\\r", "", h); hl=tolower(h)
                        if (hl=="snp" || hl=="rs" || hl=="id" || hl=="variant") s=i
                        if (hl=="chr" || hl=="chrom" || hl=="chromosome") c=i
                        if (hl=="pos" || hl=="bp" || hl=="position") p=i
                        if (hl=="p" || hl=="pvalue" || hl=="p_value" || hl=="pval") pv=i
                    }
                    next
                }
                s && pv {
                    chr = c ? \$c : "NA"
                    pos = p ? \$p : "NA"
                    print chr, pos, \$s, \$pv
                }
            ' "\${ASSOC}" >> "\${OUT}"
            ;;
        *)
            # Generic: try common header names
            awk -v OFS='\\t' '
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        h=\$i; gsub("\\r", "", h); hl=tolower(h)
                        if (hl=="snp" || hl=="rs" || hl=="marker" || hl=="id") s=i
                        if (hl=="chr" || hl=="chrom" || hl=="chromosome" || hl=="chrs") c=i
                        if (hl=="pos" || hl=="bp" || hl=="ps" || hl=="position") p=i
                        if (hl=="p" || hl=="pvalue" || hl=="p_value" || hl=="p_lrt" || hl=="p_wald" || hl ~ "p\\\\.value") pv=i
                    }
                    next
                }
                s && pv {
                    chr = c ? \$c : "NA"
                    pos = p ? \$p : "NA"
                    print chr, pos, \$s, \$pv
                }
            ' "\${ASSOC}" >> "\${OUT}"
            ;;
    esac
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${tool}"
    """
    echo -e "CHR\tPOS\tSNP\tP" > ${prefix}.std.tsv
    echo -e "1\t100\trsStub1\t0.001" >> ${prefix}.std.tsv
    echo -e "2\t200\trsStub2\t0.05" >> ${prefix}.std.tsv
    """
}
