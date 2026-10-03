#!/usr/bin/env python3
"""Annotation enrichment prior (TORUS-style) from SNP z-scores.

Official TORUS is not on bioconda. Input: snp, z, plus binary/numeric annot columns.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def irls_logit(x: list[list[float]], y: list[float], max_iter: int = 50) -> list[float]:
    n = len(y)
    p = len(x[0])
    beta = [0.0] * p
    for _ in range(max_iter):
        w = []
        z = []
        mu = []
        for i in range(n):
            eta = sum(beta[j] * x[i][j] for j in range(p))
            m = sigmoid(eta)
            m = min(1.0 - 1e-12, max(1e-12, m))
            mu.append(m)
            wt = m * (1.0 - m)
            w.append(wt)
            z.append(eta + (y[i] - m) / wt)
        # solve (X'WX) b = X'Wz
        xtwx = [[0.0] * p for _ in range(p)]
        xtwz = [0.0] * p
        for i in range(n):
            for j in range(p):
                xtwz[j] += x[i][j] * w[i] * z[i]
                for k in range(p):
                    xtwx[j][k] += x[i][j] * w[i] * x[i][k]
        # Gaussian elimination with intercept ridge
        for j in range(p):
            xtwx[j][j] += 1e-8
        a = [row[:] + [xtwz[i]] for i, row in enumerate(xtwx)]
        for col in range(p):
            pivot = max(range(col, p), key=lambda r: abs(a[r][col]))
            a[col], a[pivot] = a[pivot], a[col]
            if abs(a[col][col]) < 1e-18:
                continue
            div = a[col][col]
            for k in range(col, p + 1):
                a[col][k] /= div
            for r in range(p):
                if r == col:
                    continue
                fac = a[r][col]
                for k in range(col, p + 1):
                    a[r][k] -= fac * a[col][k]
        new = [a[j][p] for j in range(p)]
        if max(abs(new[j] - beta[j]) for j in range(p)) < 1e-8:
            return new
        beta = new
    return beta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annot", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--z-threshold", type=float, default=2.0)
    args = parser.parse_args()

    path = Path(args.annot)
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    reserved = {"snp", "variant_id", "rsid", "id", "z", "beta", "se", "p", "pval"}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit("empty input")
        raw_names = list(reader.fieldnames)
        lower = {c.lower().strip(): c for c in raw_names}
        snp_key = next((lower[k] for k in ("snp", "variant_id", "rsid", "id") if k in lower), None)
        z_key = lower.get("z")
        if not snp_key or not z_key:
            raise SystemExit("need snp and z")
        annot_names = [c for c in raw_names if c.lower().strip() not in reserved]
        snps: list[str] = []
        zs: list[float] = []
        anns: list[list[float]] = []
        for row in reader:
            snp = (row.get(snp_key) or "").strip()
            try:
                z = float(row[z_key])
            except (TypeError, ValueError, KeyError):
                continue
            if not snp:
                continue
            vec = []
            ok = True
            for name in annot_names:
                try:
                    vec.append(float(row.get(name) or 0.0))
                except (TypeError, ValueError):
                    ok = False
                    break
            if not ok:
                continue
            snps.append(snp)
            zs.append(z)
            anns.append(vec)

    if len(snps) < 4:
        raise SystemExit("need >=4 SNPs")
    if not annot_names:
        annot_names = ["intercept_only"]
        anns = [[] for _ in snps]

    y = [1.0 if abs(z) >= args.z_threshold else 0.0 for z in zs]
    x = [[1.0] + row for row in anns]
    beta = irls_logit(x, y)
    names = ["intercept"] + annot_names

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    prior_path = out_prefix.with_suffix(".prior.tsv")
    with prior_path.open("w") as fh:
        fh.write("snp\tz\tprior\n")
        for snp, z, row in zip(snps, zs, x):
            eta = sum(beta[j] * row[j] for j in range(len(beta)))
            fh.write(f"{snp}\t{z:.8g}\t{sigmoid(eta):.8g}\n")

    enr_path = out_prefix.with_suffix(".enrich.tsv")
    with enr_path.open("w") as fh:
        fh.write("term\tlog_odds\n")
        for name, b in zip(names, beta):
            fh.write(f"{name}\t{b:.8g}\n")
    log = out_prefix.with_suffix(".torus.log")
    log.write_text(f"nsnps\t{len(snps)}\nn_annot\t{len(annot_names)}\nz_threshold\t{args.z_threshold}\n")


if __name__ == "__main__":
    main()
