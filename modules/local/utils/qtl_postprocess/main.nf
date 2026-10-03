process QTL_POSTPROCESS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(cis_qtl)

    output:
    tuple val(meta), path("*.cis_std.tsv"), emit: standardized
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def engine = meta.engine ?: 'unknown'
    """
    python3 - <<'PY'
import csv, gzip, math, sys
from pathlib import Path

src = Path("${cis_qtl}")
out = Path("${prefix}.cis_std.tsv")
engine = "${engine}"

def open_maybe_gz(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")

def bh_qvalues(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        p = pvals[i]
        adj = min(prev, p * n / (n - rank + 1))
        q[i] = adj
        prev = adj
    return [min(1.0, x) for x in q]

rows = []
with open_maybe_gz(src) as fh:
    sample = fh.read(4096)
    fh.seek(0)
    dialect = csv.Sniffer().sniff(sample, delimiters="\\t, ")
    reader = csv.DictReader(fh, dialect=dialect)
    if reader.fieldnames is None:
        sys.exit(0)
    fields = [f.strip().lstrip("#") for f in reader.fieldnames]
    remap = dict(zip(reader.fieldnames, fields))
    def get(row, names):
        lower = {remap.get(k, k).lower(): v for k, v in row.items()}
        for n in names:
            if n in lower and lower[n] not in (None, "", "NA", "nan"):
                return lower[n]
        return ""
    for row in reader:
        gene = get(row, ["gene", "gene_id", "phenotype_id", "phe_id", "pheno_id", "molecular_trait_id"])
        snp = get(row, ["snp", "variant_id", "rsid", "var_id", "rs", "id"])
        chrom = get(row, ["chr", "chrom", "chromosome", "phe_chr", "var_chr", "chrom_var"])
        pos = get(row, ["pos", "position", "ps", "var_from", "start", "bp"])
        p = get(row, ["p", "pval", "p_value", "pvalue", "nom_pval", "p_nominal", "pval_nominal", "p_wald"])
        beta = get(row, ["beta", "slope", "b", "effect"])
        se = get(row, ["se", "slope_se", "std_error", "stderr"])
        try:
            pnum = float(p)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(pnum) or pnum < 0:
            continue
        rows.append((gene, snp, chrom, pos, pnum, beta, se))

pvals = [r[4] for r in rows]
qs = bh_qvalues(pvals)
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\\t")
    w.writerow(["gene", "variant_id", "chr", "pos", "p", "q", "beta", "se", "engine"])
    for (gene, snp, chrom, pos, pnum, beta, se), q in zip(rows, qs):
        w.writerow([gene, snp, chrom, pos, f"{pnum:.8g}", f"{q:.8g}", beta, se, engine])
PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def engine = meta.engine ?: 'unknown'
    """
    echo -e "gene\\tvariant_id\\tchr\\tpos\\tp\\tq\\tbeta\\tse\\tengine" > ${prefix}.cis_std.tsv
    echo -e "gene1\\tsnp1\\tchr1\\t100\\t1e-8\\t1e-8\\t0.4\\t0.08\\t${engine}" >> ${prefix}.cis_std.tsv
    """
}
