process MATRIXEQTL_CIS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Official bioconda MatrixEQTL packaging is unused; pin YangmingSi + bioinfortools 2.4.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/matrixeqtl:2.4' :
        'quay.io/bioinfortools/matrixeqtl:2.4' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(covariates)

    output:
    tuple val(meta), path("*.cis_qtl*.txt.gz"), emit: cis_qtl, optional: true
    tuple val(meta), path("*_out"), emit: outdir
    tuple val("${task.process}"), val('matrixeqtl'), eval('Rscript -e \'cat(as.character(packageVersion("MatrixEQTL")))\' 2>/dev/null || echo 2.4'), emit: versions_matrixeqtl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--cis_window 1000000 --pv_cis 1 --pv_trans 0 --model linear'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cov_path = covariates ? "${covariates}" : ''
    """
    GENO_PREFIX=\$(echo ${bed} | sed 's/\\.bed\$//')
    mkdir -p ${prefix}_out ${prefix}_meqtl

    plink2 \\
        --bfile \${GENO_PREFIX} \\
        --export A-transpose \\
        --out ${prefix}_meqtl/geno \\
        --threads ${task.cpus}

    python3 - <<'PY'
import csv, gzip, sys
from pathlib import Path

prefix = "${prefix}"
pheno = Path("${phenotype}")
cov_in = Path("${cov_path}") if "${cov_path}" else None
outdir = Path(f"{prefix}_meqtl")
bim = Path("${bim}")

# --- genotype: plink2 .traw → MatrixEQTL SNP matrix + snpsloc ---
traw = outdir / "geno.traw"
with traw.open() as fh:
    header = fh.readline().rstrip("\\n").split("\\t")
# CHR SNP (C)M POS COUNTED ALT <samples...>
sample_cols = header[6:]
samples = [c.split("_")[-1] if "_" in c else c for c in sample_cols]
# plink2 A-transpose sample headers are often FID_IID
samples = []
for c in sample_cols:
    parts = c.split("_")
    samples.append(parts[-1] if len(parts) >= 2 else c)

snp_path = outdir / "SNP.txt"
snploc_path = outdir / "snpsloc.txt"
with traw.open() as fh, snp_path.open("w", newline="") as snp_out, snploc_path.open("w", newline="") as loc_out:
    reader = csv.DictReader(fh, delimiter="\\t")
    snp_w = csv.writer(snp_out, delimiter="\\t", lineterminator="\\n")
    loc_w = csv.writer(loc_out, delimiter="\\t", lineterminator="\\n")
    snp_w.writerow(["id"] + samples)
    loc_w.writerow(["snp", "chr", "pos"])
    for row in reader:
        snp = row["SNP"]
        chrom = row["CHR"].replace("chr", "")
        pos = row["POS"]
        dosages = []
        for col, sid in zip(sample_cols, samples):
            val = row.get(col, "NA")
            if val in ("", "NA", None):
                dosages.append("NA")
            else:
                try:
                    dosages.append(str(float(val)))
                except ValueError:
                    dosages.append("NA")
        snp_w.writerow([snp] + dosages)
        loc_w.writerow([snp, chrom, pos])

# --- phenotype BED → GE + geneloc ---
def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()

with open_maybe_gz(pheno) as fh:
    rows = list(csv.reader(fh, delimiter="\\t"))
if not rows:
    sys.exit("Empty phenotype BED")
header = rows[0]
# FastQTL: #chr start end gene_id ... samples  OR chr start end gene_id strand ...
if header[0].startswith("#"):
    header[0] = header[0].lstrip("#")
# detect whether first row is header
hdr0 = header[0].lower().lstrip("#")
has_header = hdr0 in ("chr", "chrom") or (len(header) > 3 and header[3].lower() in ("gene_id", "gene", "id", "phenotype_id", "pheno_id"))
gene_idx = 3
sample_start = 4
if has_header:
    data = rows[1:]
    # FastQTL: chrom start end gene_id [strand] samples...
    if len(header) > 4 and header[4].lower() == "strand":
        sample_start = 5
    pheno_samples = header[sample_start:]
else:
    data = rows
    n_samp = len(rows[0]) - 4
    pheno_samples = samples[:n_samp] if n_samp == len(samples) else [f"S{i}" for i in range(n_samp)]

# Align to genotype sample order when possible
samp_index = {s: i for i, s in enumerate(pheno_samples)}
ordered = [s for s in samples if s in samp_index]
if not ordered:
    ordered = pheno_samples
    samples = ordered

ge_path = outdir / "GE.txt"
geneloc_path = outdir / "geneloc.txt"
with ge_path.open("w", newline="") as ge_out, geneloc_path.open("w", newline="") as gl_out:
    ge_w = csv.writer(ge_out, delimiter="\\t", lineterminator="\\n")
    gl_w = csv.writer(gl_out, delimiter="\\t", lineterminator="\\n")
    ge_w.writerow(["id"] + ordered)
    gl_w.writerow(["geneid", "chr", "s1", "s2"])
    for row in data:
        if len(row) <= sample_start:
            continue
        chrom, start, end, gene = row[0], row[1], row[2], row[gene_idx]
        vals = row[sample_start:]
        if len(vals) < len(pheno_samples):
            continue
        mapped = []
        for s in ordered:
            i = samp_index.get(s)
            mapped.append(vals[i] if i is not None else "NA")
        ge_w.writerow([gene] + mapped)
        gl_w.writerow([gene, chrom.replace("chr", ""), start, end])

# --- covariates → MatrixEQTL (cov x sample) ---
cov_path = outdir / "Covariates.txt"
have_cov = False
if cov_in and cov_in.exists() and cov_in.stat().st_size > 0:
    with cov_in.open() as fh:
        cov_rows = list(csv.reader(fh, delimiter="\\t"))
    if cov_rows:
        ch = cov_rows[0]
        # sample×cov (tensorqtl): row0 = sample id; cov×sample (omiga): col headers = samples
        body = cov_rows[1:]
        sample_x_cov = bool(body) and body[0][0] in set(samples)
        if sample_x_cov:
            cov_names = ch[1:]
            sample_to_vals = {r[0]: r[1:] for r in body if r}
            with cov_path.open("w", newline="") as out:
                w = csv.writer(out, delimiter="\\t", lineterminator="\\n")
                w.writerow(["id"] + ordered)
                for j, cname in enumerate(cov_names):
                    row = [cname]
                    for s in ordered:
                        vals = sample_to_vals.get(s, [])
                        row.append(vals[j] if j < len(vals) else "NA")
                    w.writerow(row)
            have_cov = True
        else:
            cov_samples = ch[1:]
            idx = {s: i for i, s in enumerate(cov_samples)}
            with cov_path.open("w", newline="") as out:
                w = csv.writer(out, delimiter="\\t", lineterminator="\\n")
                w.writerow(["id"] + ordered)
                for r in body:
                    if not r:
                        continue
                    row = [r[0]]
                    for s in ordered:
                        i = idx.get(s)
                        row.append(r[i + 1] if i is not None and i + 1 < len(r) else "NA")
                    w.writerow(row)
            have_cov = True

Path("have_cov.txt").write_text("1" if have_cov else "0")
Path("samples.txt").write_text("\\n".join(ordered) + "\\n")
PY

    COV_ARG=""
    if [ "\$(cat have_cov.txt)" = "1" ]; then
        COV_ARG="--covariates_file ${prefix}_meqtl/Covariates.txt"
    fi

    matrixeqtl \\
        --SNP_file ${prefix}_meqtl/SNP.txt \\
        --exp_file ${prefix}_meqtl/GE.txt \\
        --snps_loc ${prefix}_meqtl/snpsloc.txt \\
        --gene_loc ${prefix}_meqtl/geneloc.txt \\
        --output_prefix ${prefix}_out/${prefix} \\
        --threads ${task.cpus} \\
        \${COV_ARG} \\
        ${args}

    if [ -f ${prefix}_out/${prefix}_eQTL_cis.txt.gz ]; then
        cp ${prefix}_out/${prefix}_eQTL_cis.txt.gz ${prefix}.cis_qtl.txt.gz
    elif [ -f ${prefix}_out/${prefix}_eQTL_cis.txt ]; then
        gzip -c ${prefix}_out/${prefix}_eQTL_cis.txt > ${prefix}.cis_qtl.txt.gz
    else
        # Standardize header for postprocess even if empty
        echo -e "gene\\tsnp\\tbeta\\tt-stat\\tpvalue\\tFDR" | gzip -c > ${prefix}.cis_qtl.txt.gz
    fi
    # MatrixEQTL cis columns: snps gene beta t-stat p-value FDR → rename for postprocess
    python3 - <<'PY'
import gzip
from pathlib import Path
src = Path("${prefix}.cis_qtl.txt.gz")
lines = gzip.open(src, "rt").read().splitlines()
if not lines:
    raise SystemExit(0)
header = lines[0].lstrip().split("\\t")
lower = [h.lower() for h in header]
# Map MatrixEQTL names
rename = {
    "snps": "snp",
    "gene": "gene",
    "beta": "beta",
    "t-stat": "tstat",
    "p-value": "pvalue",
    "fdr": "FDR",
}
out_header = [rename.get(h.lower(), h) for h in header]
with gzip.open(src, "wt") as fh:
    fh.write("\\t".join(out_header) + "\\n")
    for line in lines[1:]:
        fh.write(line + "\\n")
PY
    cp ${prefix}.cis_qtl.txt.gz ${prefix}_out/
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out
    echo -e "gene\\tsnp\\tbeta\\ttstat\\tpvalue\\tFDR" > ${prefix}.cis_qtl.txt
    echo -e "gene1\\tsnp1\\t0.4\\t5.0\\t1e-8\\t1e-6" >> ${prefix}.cis_qtl.txt
    gzip -c ${prefix}.cis_qtl.txt > ${prefix}.cis_qtl.txt.gz
    cp ${prefix}.cis_qtl.txt.gz ${prefix}_out/
    """
}
