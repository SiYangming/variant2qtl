#!/usr/bin/env python3
"""FUSION/PrediXcan-style TWAS: gene weights x GWAS z (independent SNPs).

Official PrediXcan/FUSION binaries are not required. Columns:
  weights: gene, snp, weight
  gwas: snp, beta, se (or z)
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def p_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def read_table(path: Path) -> list[dict[str, str]]:
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit(f"empty table: {path}")
        rows = []
        for raw in reader:
            rows.append({(k or "").lower().strip(): v for k, v in raw.items()})
        return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    gwas_z: dict[str, float] = {}
    for row in read_table(Path(args.gwas)):
        snp = (row.get("snp") or row.get("variant_id") or row.get("rsid") or "").strip()
        if not snp:
            continue
        if row.get("z"):
            try:
                gwas_z[snp] = float(row["z"])
                continue
            except ValueError:
                pass
        try:
            beta = float(row["beta"])
            se = float(row["se"])
        except (KeyError, TypeError, ValueError):
            continue
        if se <= 0:
            continue
        gwas_z[snp] = beta / se

    by_gene: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for row in read_table(Path(args.weights)):
        gene = (row.get("gene") or row.get("feature") or "").strip()
        snp = (row.get("snp") or row.get("variant_id") or row.get("rsid") or "").strip()
        try:
            w = float(row.get("weight") or row.get("w") or "")
        except ValueError:
            continue
        if gene and snp:
            by_gene[gene].append((snp, w))

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = out_prefix.with_suffix(".twas.tsv")
    with tsv.open("w") as fh:
        fh.write("gene\tn_snps\tz\tp\n")
        for gene in sorted(by_gene):
            num = 0.0
            den = 0.0
            n = 0
            for snp, w in by_gene[gene]:
                if snp not in gwas_z:
                    continue
                num += w * gwas_z[snp]
                den += w * w
                n += 1
            if n < 1 or den <= 0:
                continue
            z = num / math.sqrt(den)
            fh.write(f"{gene}\t{n}\t{z:.8g}\t{p_two_sided(z):.8g}\n")
    out_prefix.with_suffix(".twas.log").write_text(f"ngenes\t{len(by_gene)}\nngwas\t{len(gwas_z)}\n")


if __name__ == "__main__":
    main()
