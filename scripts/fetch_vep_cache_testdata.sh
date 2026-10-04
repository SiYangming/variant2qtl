#!/usr/bin/env bash
# Fetch a tiny Ensembl VEP cache for real (non-stub) smoke tests.
# Source: Seqera annotation-cache (same as nf-core ensemblvep tests).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${ROOT}/assets/testdata/vep_cache_real"
CACHE_NAME="115_WBcel235"
mkdir -p "${DEST}"
if [[ -d "${DEST}/${CACHE_NAME}" ]]; then
  echo "VEP cache already present: ${DEST}/${CACHE_NAME}"
  exit 0
fi
echo "Downloading ${CACHE_NAME} from s3://annotation-cache (anonymous)..."
# Prefer AWS CLI if present; else curl via https gateway if available.
if command -v aws >/dev/null 2>&1; then
  aws s3 sync --no-sign-request "s3://annotation-cache/vep_cache/${CACHE_NAME}/" "${DEST}/${CACHE_NAME}/"
else
  echo "aws CLI not found; trying nextflow/s3 via docker pull of a tiny helper is not used."
  echo "Install aws CLI or manually sync:"
  echo "  aws s3 sync --no-sign-request s3://annotation-cache/vep_cache/${CACHE_NAME}/ ${DEST}/${CACHE_NAME}/"
  exit 1
fi
echo "Done: ${DEST}/${CACHE_NAME}"
