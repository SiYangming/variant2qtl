#!/usr/bin/env python3
"""FUSION/PrediXcan/MetaXcan-style TWAS: gene weights x GWAS z.

Supports:
  - long weights: gene|GENE, snp|rsid|variant_id, weight|WEIGHT
  - optional LD matrix TSV (first column snp ids, remaining columns matching SNPs)
    for FUSION-style variance w'R w; without LD, independent-SNP w'w.
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


def read_ld(path: Path) -> dict[str, dict[str, float]]:
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    mat: dict[str, dict[str, float]] = {}
    with path.open() as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader, None)
        if not header:
            raise SystemExit("empty LD matrix")
        cols = [h.strip() for h in header[1:]]
        for row in reader:
            if not row:
                continue
            rid = row[0].strip()
            mat[rid] = {}
            for j, col in enumerate(cols):
                try:
                    mat[rid][col] = float(row[j + 1])
                except (IndexError, ValueError):
                    mat[rid][col] = 0.0
    return mat


def gene_z(weights: list[tuple[str, float]], gwas_z: dict[str, float], ld: dict[str, dict[str, float]] | None) -> tuple[float, int] | None:
    pairs = [(snp, w) for snp, w in weights if snp in gwas_z]
    if not pairs:
        return None
    num = sum(w * gwas_z[snp] for snp, w in pairs)
    if ld is None:
        den = sum(w * w for _snp, w in pairs)
    else:
        den = 0.0
        for snp_i, w_i in pairs:
            for snp_j, w_j in pairs:
                r = 1.0 if snp_i == snp_j else ld.get(snp_i, {}).get(snp_j, ld.get(snp_j, {}).get(snp_i, 0.0))
                den += w_i * w_j * r
    if den <= 0:
        return None
    return num / math.sqrt(den), len(pairs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--ld", default=None)
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
        gene = (row.get("gene") or row.get("feature") or row.get("gene_name") or "").strip()
        snp = (row.get("snp") or row.get("variant_id") or row.get("rsid") or row.get("snpid") or "").strip()
        w_raw = row.get("weight") or row.get("w") or row.get("effect") or ""
        try:
            w = float(w_raw)
        except ValueError:
            continue
        if gene and snp:
            by_gene[gene].append((snp, w))

    ld = read_ld(Path(args.ld)) if args.ld else None

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = out_prefix.with_suffix(".twas.tsv")
    with tsv.open("w") as fh:
        fh.write("gene\tn_snps\tz\tp\tld\n")
        for gene in sorted(by_gene):
            res = gene_z(by_gene[gene], gwas_z, ld)
            if res is None:
                continue
            z, n = res
            fh.write(f"{gene}\t{n}\t{z:.8g}\t{p_two_sided(z):.8g}\t{'yes' if ld else 'no'}\n")
    out_prefix.with_suffix(".twas.log").write_text(
        f"ngenes\t{len(by_gene)}\nngwas\t{len(gwas_z)}\nld\t{'yes' if ld else 'no'}\n"
    )


if __name__ == "__main__":
    main()
