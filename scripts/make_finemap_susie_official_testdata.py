#!/usr/bin/env python3
"""Generate SuSiE sumstats + LD matching the susieR vignette simulation recipe.

Prefer scripts/make_finemap_susie_official_testdata.R when susieR is installed
(uses N3finemapping when available). This Python fallback keeps the repo
regenerable without R.

Usage:
  python3 scripts/make_finemap_susie_official_testdata.py
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--out-dir",
        default="assets/testdata/finemap_susie_official",
        help="Output directory",
    )
    parser.add_argument("--n", type=int, default=200, help="Sample size (vignette uses 200)")
    parser.add_argument("--p", type=int, default=100, help="Variant count (trimmed from vignette 1000)")
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(1)
    n, p = args.n, args.p
    beta = np.zeros(p)
    beta[:4] = 1.0
    X = rng.standard_normal((n, p))
    X = X - X.mean(axis=0, keepdims=True)
    y = X @ beta + rng.standard_normal(n)
    X = (X - X.mean(axis=0, keepdims=True)) / X.std(axis=0, keepdims=True, ddof=1)
    y = y - y.mean()

    # Univariate OLS z-scores
    beta_hat = np.empty(p)
    se_hat = np.empty(p)
    z = np.empty(p)
    for j in range(p):
        x = X[:, j]
        b = float(np.dot(x, y) / np.dot(x, x))
        resid = y - b * x
        sigma2 = float(np.dot(resid, resid) / (n - 2))
        se = float(np.sqrt(sigma2 / np.dot(x, x)))
        beta_hat[j] = b
        se_hat[j] = se
        z[j] = b / se

    R = np.corrcoef(X, rowvar=False)
    np.fill_diagonal(R, 1.0)
    R = np.clip(R, -1.0, 1.0)

    ids = [f"var{i+1}" for i in range(p)]
    with (out / "sumstats.tsv").open("w", encoding="utf-8") as fh:
        fh.write("variant_id\tbeta\tse\tz\tn\n")
        for i in range(p):
            fh.write(f"{ids[i]}\t{beta_hat[i]:.8g}\t{se_hat[i]:.8g}\t{z[i]:.8g}\t{n}\n")
    np.savetxt(out / "ld.txt", R, fmt="%.8f", delimiter="\t")

    source = (
        f"susieR vignette simulation recipe (python fallback; n={n},p={p},causal=1:4)\n"
        f"generated: {datetime.now(timezone.utc).isoformat()}\n"
    )
    (out / "SOURCE.txt").write_text(source, encoding="utf-8")
    print(f"Wrote {out}/sumstats.tsv and ld.txt ({p} variants)")


if __name__ == "__main__":
    main()
