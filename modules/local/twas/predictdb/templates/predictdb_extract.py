#!/usr/bin/env python3
"""Extract PredictDB / PrediXcan SQLite weights into long TSV (gene, snp, weight).

Expects PredictDB-style tables:
  weights(rsid|varID, gene, weight, ...)
  optional extra(gene, genename, ...) for gene symbol aliases
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    db = Path(args.db)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "weights" not in tables:
        raise SystemExit(f"PredictDB missing weights table: {db}")

    gene_alias: dict[str, str] = {}
    if "extra" in tables:
        cols = {r[1].lower() for r in con.execute("PRAGMA table_info(extra)")}
        name_col = "genename" if "genename" in cols else ("gene_name" if "gene_name" in cols else None)
        if name_col:
            for row in con.execute(f"SELECT gene, {name_col} AS gname FROM extra"):
                gene = (row["gene"] or "").strip()
                gname = (row["gname"] or "").strip()
                if gene and gname:
                    gene_alias[gene] = gname

    wcols = {r[1].lower(): r[1] for r in con.execute("PRAGMA table_info(weights)")}
    snp_col = None
    for cand in ("rsid", "varid", "snp", "variant_id"):
        if cand in wcols:
            snp_col = wcols[cand]
            break
    if snp_col is None:
        raise SystemExit("weights table needs rsid/varID/snp column")
    if "gene" not in wcols or "weight" not in wcols:
        raise SystemExit("weights table needs gene and weight columns")

    gene_col = wcols["gene"]
    weight_col = wcols["weight"]

    n = 0
    with out.open("w") as fh:
        fh.write("gene\tsnp\tweight\tgene_id\n")
        for row in con.execute(f"SELECT {gene_col} AS gene, {snp_col} AS snp, {weight_col} AS weight FROM weights"):
            gene_id = (row["gene"] or "").strip()
            snp = (row["snp"] or "").strip()
            if not gene_id or not snp:
                continue
            try:
                w = float(row["weight"])
            except (TypeError, ValueError):
                continue
            gene = gene_alias.get(gene_id, gene_id)
            fh.write(f"{gene}\t{snp}\t{w:.8g}\t{gene_id}\n")
            n += 1

    con.close()
    if n == 0:
        raise SystemExit(f"no weights extracted from {db}")
    out.with_suffix(out.suffix + ".log").write_text(f"nweights\t{n}\n")


if __name__ == "__main__":
    main()
