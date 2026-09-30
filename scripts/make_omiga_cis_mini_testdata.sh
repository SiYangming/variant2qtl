#!/usr/bin/env bash
# Regenerate assets/testdata/omiga_cis_mini from nf-core plink_simulated.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/assets/testdata/omiga_cis_mini"
BASE_URL="${MODULES_TESTDATA_BASE_PATH:-https://raw.githubusercontent.com/nf-core/test-datasets/modules/data}"
POPGEN="${BASE_URL}/genomics/homo_sapiens/popgen"
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/omiga_cis_mini.XXXXXX")"
trap 'rm -rf "${WORKDIR}"' EXIT

mkdir -p "${OUT}"
curl -fsSL "${POPGEN}/plink_simulated.bed" -o "${OUT}/geno.bed"
curl -fsSL "${POPGEN}/plink_simulated.bim" -o "${WORKDIR}/plink_simulated.bim"
curl -fsSL "${POPGEN}/plink_simulated.fam" -o "${OUT}/geno.fam"
curl -fsSL "${POPGEN}/plink_simulated.fam" -o "${WORKDIR}/plink_simulated.fam"
curl -fsSL "${POPGEN}/plink_simulated_covariates.txt" -o "${WORKDIR}/plink_simulated_covariates.txt"
# plink_simulated alleles are D/d; OmiGA expects ACGT — recode labels only (bed dosages unchanged)
awk 'BEGIN{OFS="\t"} {
  a5=$5; a6=$6
  gsub(/D/,"A",a5); gsub(/d/,"T",a5)
  gsub(/D/,"A",a6); gsub(/d/,"T",a6)
  print $1,$2,$3,$4,a5,a6
}' "${WORKDIR}/plink_simulated.bim" > "${OUT}/geno.bim"

python3 - "${WORKDIR}" "${OUT}" <<'PY'
import gzip
import random
import sys
from pathlib import Path

workdir, outdir = Path(sys.argv[1]), Path(sys.argv[2])
fam = (workdir / "plink_simulated.fam").read_text().strip().splitlines()
cov_lines = (workdir / "plink_simulated_covariates.txt").read_text().strip().splitlines()
iids = [ln.split()[1] for ln in fam]
cov_header = cov_lines[0].split()
cov_by_iid = {r[1]: r for r in (ln.split() for ln in cov_lines[1:])}
missing = [i for i in iids if i not in cov_by_iid]
if missing:
    raise SystemExit(f"covariates missing IIDs: {missing[:5]}")

# BED phenotype: #chr start end pheno_id + samples (start = 0-based TSS near plink_simulated SNP coords)
genes = [
    ("1", 49, 50, "GENE_A"),
    ("1", 99, 100, "GENE_B"),
    ("1", 149, 150, "GENE_C"),
]
rng = random.Random(42)
pheno = outdir / "phenotype.bed.gz"
with gzip.open(pheno, "wt") as fh:
    fh.write("\t".join(["#chr", "start", "end", "pheno_id"] + iids) + "\n")
    for gi, (chrom, start, end, gid) in enumerate(genes):
        vals = []
        for iid in iids:
            sex = float(cov_by_iid[iid][2])
            vals.append(f"{rng.gauss(0.0, 1.0) + 0.15 * sex + 0.05 * gi:.6f}")
        fh.write("\t".join([chrom, str(start), str(end), gid] + vals) + "\n")

cov_names = ["Sex", "Age", "PC1", "PC2"]
cov_idx = [cov_header.index(n) for n in cov_names]
cov_path = outdir / "covariates.txt"
with cov_path.open("w") as fh:
    fh.write("\t".join(["id"] + iids) + "\n")
    for name, idx in zip(cov_names, cov_idx):
        fh.write("\t".join([name] + [cov_by_iid[iid][idx] for iid in iids]) + "\n")

(outdir / "README.md").write_text(
    """# OmiGA cis mini testdata

Tiny molecular phenotype + covariates + genotype derived from nf-core `plink_simulated`
(200 samples, chr1 SNPs at positions 1–220).

| File | Description |
|------|-------------|
| `geno.{bed,bim,fam}` | Copy of `plink_simulated` PLINK binary set |
| `phenotype.bed.gz` | FastQTL-style BED: `#chr start end pheno_id` + per-sample values (3 fake genes) |
| `covariates.txt` | OmiGA default orientation (covariate × sample); Sex/Age/PC1/PC2 |

Regenerate:

```bash
bash scripts/make_omiga_cis_mini_testdata.sh
```
"""
)
print(f"Wrote {outdir}/geno.bed, geno.bim, geno.fam")
print(f"Wrote {pheno} ({pheno.stat().st_size} bytes)")
print(f"Wrote {cov_path} ({cov_path.stat().st_size} bytes)")
PY

echo "Done: ${OUT}"
