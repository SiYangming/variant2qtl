#!/usr/bin/env bash
# Download official tensorQTL GEUVADIS example data, convert pgen→bed, and
# subset phenotypes for a feasible CPU real-run smoke test.
#
# Source: https://github.com/broadinstitute/tensorqtl/tree/master/example/data
# (~80MB download; binaries are gitignored — do not commit).
#
# Usage:
#   bash scripts/fetch_tensorqtl_geuvadis_testdata.sh
#   N_PHENO=5 bash scripts/fetch_tensorqtl_geuvadis_testdata.sh
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${OUT_DIR:-${ROOT}/assets/testdata/tensorqtl_geuvadis}"
N_PHENO="${N_PHENO:-5}"
BASE_URL="${BASE_URL:-https://raw.githubusercontent.com/broadinstitute/tensorqtl/master/example/data}"
PFILE_PREFIX="GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18"

mkdir -p "${OUT}/raw"
cd "${OUT}"

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $1" >&2
        exit 1
    }
}

need_cmd curl
need_cmd plink2
need_cmd python3

download() {
    local name="$1"
    local dest="raw/${name}"
    if [[ -s "${dest}" ]]; then
        echo "keep ${dest}"
        return 0
    fi
    echo "download ${name}"
    curl -fsSL -o "${dest}.partial" "${BASE_URL}/${name}"
    mv "${dest}.partial" "${dest}"
}

download "${PFILE_PREFIX}.pgen"
download "${PFILE_PREFIX}.pvar"
download "${PFILE_PREFIX}.psam"
download "GEUVADIS.445_samples.expression.bed.gz"
download "GEUVADIS.445_samples.covariates.txt"

# Always (re)build bed with chr-prefixed contig names to match expression BED (chr18).
# Official notebook converts with --output-chr chrM; default plink2 bed uses bare "18".
if [[ "${FORCE_REBUILD:-0}" == "1" || ! -s geno.bed || ! -s geno.bim || ! -s geno.fam ]] ||
    ! head -1 geno.bim | awk '{exit !($1 ~ /^chr/)}'; then
    echo "plink2 pgen → bed (chr-prefixed contigs)"
    plink2 \
        --pfile "raw/${PFILE_PREFIX}" \
        --make-bed \
        --output-chr chrM \
        --out geno \
        --threads "${THREADS:-4}"
fi


# Subset expression BED to first N_PHENO gene rows (keep header).
if [[ ! -s phenotype.bed.gz ]]; then
    echo "subset expression to ${N_PHENO} phenotypes"
    python3 - <<PY
import gzip
from pathlib import Path

n = int("${N_PHENO}")
src = Path("raw/GEUVADIS.445_samples.expression.bed.gz")
dst = Path("phenotype.bed.gz")
with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
    header = fin.readline()
    if not header.startswith("#"):
        raise SystemExit("unexpected phenotype BED header")
    fout.write(header)
    written = 0
    for line in fin:
        # Prefer chr18 rows when present
        if not line.startswith("chr18") and not line.startswith("18"):
            continue
        fout.write(line)
        written += 1
        if written >= n:
            break
    if written == 0:
        fin.seek(0)
        _ = fin.readline()
        for line in fin:
            fout.write(line)
            written += 1
            if written >= n:
                break
print(f"wrote {written} phenotype rows → {dst}")
PY
fi

cp -f raw/GEUVADIS.445_samples.covariates.txt covariates.txt

cat > SOURCE.txt <<EOF
source: broadinstitute/tensorqtl example/data (GEUVADIS chr18)
url: ${BASE_URL}
pfile: ${PFILE_PREFIX}
n_pheno_subset: ${N_PHENO}
note: Full official example recommends GPU/~50GB; this subset is for CPU smoke.
EOF

cat > README.md <<'EOF'
# tensorQTL GEUVADIS official example (local cache)

Fetched by `scripts/fetch_tensorqtl_geuvadis_testdata.sh` from
[broadinstitute/tensorqtl example/data](https://github.com/broadinstitute/tensorqtl/tree/master/example/data).

| File | Role |
| --- | --- |
| `geno.{bed,bim,fam}` | PLINK 1 bed converted from official PLINK2 pgen |
| `phenotype.bed.gz` | Subset of official expression BED (default 5 chr18 genes) |
| `covariates.txt` | Official covariates (covariate × sample) |
| `raw/` | Original downloaded pgen/pvar/psam + full expression |

Large binaries are gitignored. Re-run the fetch script before tagged real nf-tests.
EOF

echo "Ready under ${OUT}"
ls -lh geno.bed geno.bim geno.fam phenotype.bed.gz covariates.txt
