#!/usr/bin/env python3
"""Intersect samples, drop high-missing genes, INV/quantile, write FastQTL BED."""
import argparse
import csv
import gzip
import math
import statistics
import sys
from pathlib import Path


def open_maybe_gz(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def ppf(p):
    """Standard-normal quantile (Acklam approximation)."""
    p = min(max(float(p), 1e-12), 1.0 - 1e-12)
    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577509590705e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068071649571e01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464858e00,
        2.938163982698783e00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    ]
    plow = 0.02425
    phigh = 1.0 - plow
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )
    q = p - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])
        * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    )


def invnorm(vals):
    n = len(vals)
    order = sorted(range(n), key=lambda i: vals[i])
    out = [0.0] * n
    for rank, i in enumerate(order, start=1):
        # Blom
        p = (rank - 0.375) / (n + 0.25)
        out[i] = ppf(p)
    return out


def quantile_normalize(matrix):
    """Bolstad QN treating columns as samples (rows = genes)."""
    if not matrix or not matrix[0]:
        return matrix
    n_genes = len(matrix)
    n_samp = len(matrix[0])
    sorted_cols = []
    orders = []
    for j in range(n_samp):
        col = [matrix[i][j] for i in range(n_genes)]
        order = sorted(range(n_genes), key=lambda i: col[i])
        orders.append(order)
        sorted_cols.append([col[i] for i in order])
    mean_sorted = [
        statistics.fmean(sorted_cols[j][r] for j in range(n_samp)) for r in range(n_genes)
    ]
    out = [[0.0] * n_samp for _ in range(n_genes)]
    for j in range(n_samp):
        for r, i in enumerate(orders[j]):
            out[i][j] = mean_sorted[r]
    return out


def parse_float(s):
    if s is None:
        return None
    t = str(s).strip()
    if t == "" or t.upper() in {"NA", "NAN", "."}:
        return None
    try:
        v = float(t)
    except ValueError:
        return None
    if not math.isfinite(v):
        return None
    return v


def load_table(path):
    with open_maybe_gz(path) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t, ")
        reader = csv.reader(fh, dialect)
        rows = list(reader)
    if not rows:
        sys.exit("empty phenotype matrix")
    header = [c.strip().lstrip("#") for c in rows[0]]
    return header, rows[1:]


def is_bed_header(header):
    h0 = header[0].lower()
    return h0 in {"chr", "chrom", "chromosome"} or (
        len(header) >= 4 and header[1].lower() in {"start", "chromstart"}
    )


def load_sample_keep(path):
    if not path:
        return None
    ids = []
    lines = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            lines.append(line.split())
    if not lines:
        return ids
    fam_like = all(len(p) >= 6 for p in lines[: min(5, len(lines))])
    for parts in lines:
        if parts[0].lower() in {"fid", "id", "iid"}:
            continue
        if fam_like:
            ids.append(parts[1])
        else:
            ids.append(parts[0])
    return ids


def load_gene_map(path):
    mapping = {}
    if not path:
        return mapping
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                parts = line.split()
            if len(parts) < 4:
                continue
            chrom, start, end, gid = parts[0], parts[1], parts[2], parts[3]
            mapping[gid] = (chrom, start, end)
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--gene-map", default="")
    ap.add_argument("--samples", default="")
    ap.add_argument("--max-missing", type=float, default=0.2)
    ap.add_argument("--transform", default="invnorm", choices=["invnorm", "quantile", "none"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    header, body = load_table(args.matrix)
    gene_map = load_gene_map(args.gene_map or None)
    keep = load_sample_keep(args.samples or None)

    if is_bed_header(header):
        # FastQTL BED: chr start end pheno_id samples...
        if len(header) < 5:
            sys.exit("BED phenotype needs chr/start/end/id + samples")
        sample_names = header[4:]
        gene_ids = []
        coords = []
        raw = []
        for row in body:
            if len(row) < 5:
                continue
            chrom, start, end, gid = row[0], row[1], row[2], row[3]
            gene_ids.append(gid)
            coords.append((chrom, start, end))
            raw.append([parse_float(x) for x in row[4 : 4 + len(sample_names)]])
    else:
        sample_names = header[1:]
        gene_ids = []
        coords = []
        raw = []
        for idx, row in enumerate(body):
            if not row:
                continue
            gid = row[0]
            gene_ids.append(gid)
            if gid in gene_map:
                coords.append(gene_map[gid])
            else:
                coords.append(("1", str(idx), str(idx + 1)))
            vals = [parse_float(x) for x in row[1 : 1 + len(sample_names)]]
            while len(vals) < len(sample_names):
                vals.append(None)
            raw.append(vals[: len(sample_names)])

    if keep is not None:
        keep_set = set(keep)
        idx = [i for i, s in enumerate(sample_names) if s in keep_set]
        sample_names = [sample_names[i] for i in idx]
        raw = [[row[i] for i in idx] for row in raw]

    n_samp = len(sample_names)
    if n_samp == 0:
        sys.exit("no samples left after intersection")

    kept_genes = []
    kept_coords = []
    kept_mat = []
    for gid, xy, row in zip(gene_ids, coords, raw):
        miss = sum(1 for v in row if v is None) / float(n_samp)
        if miss > args.max_missing:
            continue
        finite = [v for v in row if v is not None]
        fill = statistics.median(finite) if finite else 0.0
        kept_genes.append(gid)
        kept_coords.append(xy)
        kept_mat.append([fill if v is None else v for v in row])

    if not kept_genes:
        sys.exit("no genes left after missingness filter")

    if args.transform == "invnorm":
        kept_mat = [invnorm(row) for row in kept_mat]
    elif args.transform == "quantile":
        kept_mat = quantile_normalize(kept_mat)

    out = Path(args.out)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["#chr", "start", "end", "pheno_id"] + sample_names)
        for gid, (chrom, start, end), row in zip(kept_genes, kept_coords, kept_mat):
            w.writerow([chrom, start, end, gid] + [f"{v:.6g}" for v in row])


if __name__ == "__main__":
    main()
