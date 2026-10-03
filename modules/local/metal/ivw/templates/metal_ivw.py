#!/usr/bin/env python3
"""Inverse-variance meta-analysis (METAL-style) for long-format cohort sumstats.

Official METAL is not on bioconda. Columns: snp, cohort|study, beta, se.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def p_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def chi2_sf(x: float, df: int) -> float:
    if df <= 0 or x < 0:
        return 1.0
    # survival of chi^2 with integer df via incomplete-gamma regularized
    k = df / 2.0
    s = x / 2.0
    # Q(k,s) = 1 - P(k,s); series for lower gamma / Gamma when k is half-integer or integer
    if s == 0:
        return 1.0
    # continued fraction for upper incomplete / Gamma (Lentz)
    # start with series for lower P when s < k+1 else CF
    if s < k + 1:
        term = 1.0 / k
        acc = term
        for n in range(1, 200):
            term *= s / (k + n)
            acc += term
            if abs(term) < 1e-14 * acc:
                break
        p = acc * math.exp(-s + k * math.log(s) - math.lgamma(k))
        return max(0.0, min(1.0, 1.0 - p))
    a = 1.0 / s
    b = s + 1.0 - k
    c = 1e30
    d = 1.0 / b
    h = d
    for i in range(1, 200):
        an = -i * (i - k)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-30:
            d = 1e-30
        c = b + an / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    q = h * math.exp(-s + k * math.log(s) - math.lgamma(k))
    return max(0.0, min(1.0, q))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohorts", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    path = Path(args.cohorts)
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    by_snp: dict[str, list[tuple[float, float]]] = defaultdict(list)
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit("empty input")
        cols = {c.lower().strip(): c for c in reader.fieldnames}
        snp_key = next((cols[k] for k in ("snp", "variant_id", "rsid", "id") if k in cols), None)
        beta_key = cols.get("beta")
        se_key = cols.get("se")
        if not snp_key or not beta_key or not se_key:
            raise SystemExit("need snp, beta, se")
        for row in reader:
            snp = (row.get(snp_key) or "").strip()
            try:
                beta = float(row[beta_key])
                se = float(row[se_key])
            except (TypeError, ValueError, KeyError):
                continue
            if not snp or se <= 0:
                continue
            by_snp[snp].append((beta, se))

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for snp, pairs in sorted(by_snp.items()):
        wsum = 0.0
        bw = 0.0
        for beta, se in pairs:
            w = 1.0 / (se * se)
            wsum += w
            bw += w * beta
        if wsum <= 0:
            continue
        beta_m = bw / wsum
        se_m = math.sqrt(1.0 / wsum)
        z = beta_m / se_m
        q = 0.0
        for beta, se in pairs:
            q += ((beta - beta_m) / se) ** 2
        df = len(pairs) - 1
        p_het = chi2_sf(q, df) if df > 0 else 1.0
        rows.append(
            {
                "snp": snp,
                "beta": f"{beta_m:.8g}",
                "se": f"{se_m:.8g}",
                "z": f"{z:.8g}",
                "p": f"{p_two_sided(z):.8g}",
                "n_studies": str(len(pairs)),
                "q_het": f"{q:.8g}",
                "p_het": f"{p_het:.8g}",
            }
        )

    tsv = out_prefix.with_suffix(".metal.tsv")
    fields = ["snp", "beta", "se", "z", "p", "n_studies", "q_het", "p_het"]
    with tsv.open("w") as fh:
        fh.write("\t".join(fields) + "\n")
        for row in rows:
            fh.write("\t".join(row[k] for k in fields) + "\n")
    log = out_prefix.with_suffix(".metal.log")
    log.write_text(f"nsnps\t{len(rows)}\n")


if __name__ == "__main__":
    main()
