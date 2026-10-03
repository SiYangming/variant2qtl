#!/usr/bin/env python3
"""FINEMAP / CAVIAR / DAP-G style 1-causal PIPs from sumstats (identity LD).

Official FINEMAP/CAVIAR/DAP-G binaries are not on a single bioconda pin used here.
Columns: variant_id|snp, beta, se (or z). Optional annot prior for dapg.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_z(path: Path) -> list[tuple[str, float]]:
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    out = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit("empty sumstats")
        for raw in reader:
            row = {(k or "").lower().strip(): v for k, v in raw.items()}
            snp = (row.get("variant_id") or row.get("snp") or row.get("rsid") or "").strip()
            if not snp:
                continue
            z = None
            if row.get("z"):
                try:
                    z = float(row["z"])
                except ValueError:
                    z = None
            if z is None:
                try:
                    z = float(row["beta"]) / float(row["se"])
                except (KeyError, TypeError, ValueError, ZeroDivisionError):
                    continue
            out.append((snp, z))
    return out


def abf(z: float, w: float) -> float:
    # Wakefield ABF for quantitative trait, prior variance w
    v = 1.0
    r = w / (w + v)
    return math.sqrt(1.0 - r) * math.exp(0.5 * r * z * z)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--method", default="finemap", choices=["finemap", "caviar", "dapg"])
    args = parser.parse_args()

    zs = read_z(Path(args.sumstats))
    if len(zs) < 2:
        raise SystemExit("need >=2 SNPs")

    w = {"finemap": 0.04, "caviar": 0.1, "dapg": 0.05}[args.method]
    bfs = [abf(z, w) for _snp, z in zs]
    prior = 1.0 / len(zs)
    # DAP-G: annotation-free uniform prior; CAVIAR/FINEMAP: 1-causal mix with null
    scores = []
    for bf in bfs:
        if args.method == "dapg":
            scores.append(prior * bf)
        else:
            scores.append(bf)
    total = sum(scores)
    null = 1.0 if args.method != "dapg" else 0.0
    denom = total + null
    pips = [s / denom for s in scores]

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    pip_path = out_prefix.with_suffix(".pip.tsv")
    with pip_path.open("w") as fh:
        fh.write("variant_id\tz\tbf\tpip\tmethod\n")
        ranked = sorted(zip(zs, bfs, pips), key=lambda x: -x[2])
        for (snp, z), bf, pip in ranked:
            fh.write(f"{snp}\t{z:.8g}\t{bf:.8g}\t{pip:.8g}\t{args.method}\n")
    cs_path = out_prefix.with_suffix(".cs.tsv")
    acc = 0.0
    with cs_path.open("w") as fh:
        fh.write("credible_set\tvariant_id\tpip\n")
        for (snp, _z), _bf, pip in ranked:
            if acc >= 0.95:
                break
            fh.write(f"CS1\t{snp}\t{pip:.8g}\n")
            acc += pip
    out_prefix.with_suffix(".finemap.log").write_text(
        f"method\t{args.method}\nnsnps\t{len(zs)}\n"
    )


if __name__ == "__main__":
    main()
