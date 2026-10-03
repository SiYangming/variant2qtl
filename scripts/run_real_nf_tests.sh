#!/usr/bin/env bash
# Run optional real (non-stub) module nf-tests that are ignored by default CI.
# Swaps in nf-test.real.config for the duration of the run.
#
# Usage:
#   bash scripts/run_real_nf_tests.sh susie
#   bash scripts/run_real_nf_tests.sh tensorqtl
#   bash scripts/run_real_nf_tests.sh both
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

TARGET="${1:-both}"
PROFILE="${PROFILE:-test,conda}"

if [[ ! -f nf-test.config || ! -f nf-test.real.config ]]; then
    echo "ERROR: nf-test.config / nf-test.real.config missing" >&2
    exit 1
fi

cleanup() {
    if [[ -f nf-test.config.ci_backup ]]; then
        mv -f nf-test.config.ci_backup nf-test.config
    fi
}
trap cleanup EXIT

cp -f nf-test.config nf-test.config.ci_backup
cp -f nf-test.real.config nf-test.config

run_susie() {
    echo "==> SuSiE real (official-style sumstats)"
    nf-test test modules/local/susie/finemap/tests/main.susie_real.nf.test \
        --tag susie_real \
        --profile "${PROFILE}"
}

run_tensorqtl() {
    echo "==> tensorQTL real (official GEUVADIS example)"
    bash scripts/fetch_tensorqtl_geuvadis_testdata.sh
    nf-test test modules/local/tensorqtl/cis/tests/main.tensorqtl_real.nf.test \
        --tag tensorqtl_real \
        --profile "${PROFILE}"
}

case "${TARGET}" in
    susie) run_susie ;;
    tensorqtl) run_tensorqtl ;;
    both)
        run_susie
        run_tensorqtl
        ;;
    *)
        echo "Usage: $0 [susie|tensorqtl|both]" >&2
        exit 2
        ;;
esac
