#!/usr/bin/env python3
"""LeafCutter perind counts → intron ratios, INV/quantile, FastQTL BED + phenotype_group."""
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
        p = (rank - 0.375) / (n + 0.25)
        out[i] = ppf(p)
    return out


def parse_ratio(cell):
    t = str(cell).strip()
    if t == "" or t.upper() in {"NA", "NAN", "."}:
        return None
    if "/" in t:
        num, den = t.split("/", 1)
        try:
            n = float(num)
            d = float(den)
        except ValueError:
            return None
        if d <= 0:
            return None
        return n / d
    try:
        v = float(t)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def parse_intron(iid):
    parts = iid.replace("::", ":").split(":")
    if len(parts) >= 4:
        chrom, start, end, cluster = parts[0], parts[1], parts[2], ":".join(parts[3:])
        return chrom, start, end, cluster, iid
    bits = iid.split()
    if len(bits) >= 4:
        return bits[0], bits[1], bits[2], bits[3], iid
    return "1", "0", "1", "clu", iid


def load_genes(path):
    genes = []
    if not path:
        return genes
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                parts = line.split()
            if len(parts) < 4:
                continue
            chrom, start, end, gid = parts[0], int(parts[1]), int(parts[2]), parts[3]
            genes.append((chrom, start, end, gid))
    return genes


def annotate(chrom, start, end, genes):
    try:
        s, e = int(start), int(end)
    except ValueError:
        return "NA"
    c = chrom if chrom.startswith("chr") else chrom
    hits = []
    for gchrom, gs, ge, gid in genes:
        gc = gchrom if str(gchrom).startswith("chr") else str(gchrom)
        cc = c if c.startswith("chr") else "chr" + c
        gcn = gc if gc.startswith("chr") else "chr" + gc
        if cc != gcn:
            continue
        if e < gs or s > ge:
            continue
        hits.append(gid)
    return hits[0] if hits else "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--genes", default="")
    ap.add_argument("--min-ratio", type=float, default=0.01)
    ap.add_argument("--transform", default="invnorm", choices=["invnorm", "none"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open_maybe_gz(args.counts) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t, ")
        reader = csv.reader(fh, dialect)
        rows = list(reader)
    if not rows:
        sys.exit("empty counts")
    header = [c.strip().lstrip("#") for c in rows[0]]
    samples = header[1:]
    genes = load_genes(args.genes or None)

    kept = []
    for row in rows[1:]:
        if not row:
            continue
        iid = row[0]
        ratios = [parse_ratio(x) for x in row[1 : 1 + len(samples)]]
        while len(ratios) < len(samples):
            ratios.append(None)
        n_miss = sum(1 for v in ratios if v is None)
        if n_miss / float(len(samples)) > 0.5:
            continue
        finite = [v for v in ratios if v is not None]
        fill = statistics.median(finite) if finite else 0.0
        ratios = [fill if v is None else v for v in ratios]
        if statistics.fmean(ratios) < args.min_ratio:
            continue
        chrom, start, end, cluster, pheno_id = parse_intron(iid)
        gene = annotate(chrom, start, end, genes)
        kept.append((chrom, start, end, pheno_id, gene, cluster, ratios))

    if not kept:
        sys.exit("no introns left after filtering")

    if args.transform == "invnorm":
        for i, rec in enumerate(kept):
            chrom, start, end, pheno_id, gene, cluster, ratios = rec
            kept[i] = (chrom, start, end, pheno_id, gene, cluster, invnorm(ratios))

    out_prefix = Path(args.out)
    bed = Path(str(out_prefix) + ".phenotype.bed")
    grp = Path(str(out_prefix) + ".phenotype_group.tsv")
    with bed.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["#chr", "start", "end", "pheno_id"] + samples)
        for chrom, start, end, pheno_id, _gene, _clu, ratios in kept:
            w.writerow([chrom, start, end, pheno_id] + [f"{v:.6g}" for v in ratios])
    with grp.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["phenotype_id", "group_id"])
        for _c, _s, _e, pheno_id, gene, cluster, _r in kept:
            gid = gene if gene != "NA" else cluster
            w.writerow([pheno_id, gid])


if __name__ == "__main__":
    main()
