#!/usr/bin/env python3
"""LD-score regression style heritability (Bulik-Sullivan) on summary stats.

Official ldsc is not pinned here. Columns:
  sumstats: snp, z|beta+se, n|n_samples; optional l2 (LD score)
  ldscores (optional): reference-panel snp + l2 / L2 (merged by SNP id)
  annot (optional): snp + binary/numeric annotation columns for partitioned h2
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit(f"empty table: {path}")
        return [{(k or "").lower().strip(): v for k, v in raw.items()} for raw in reader]


def ols(x: list[float], y: list[float]) -> tuple[float, float, float]:
    n = len(x)
    if n < 3:
        raise SystemExit("need >=3 SNPs for LDSC")
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    if sxx <= 0:
        raise SystemExit("LD scores have zero variance")
    slope = sxy / sxx
    intercept = my - slope * mx
    resid = [(yi - intercept - slope * xi) for xi, yi in zip(x, y)]
    sse = sum(r * r for r in resid)
    se = math.sqrt(sse / (n - 2) / sxx) if n > 2 else float("nan")
    return intercept, slope, se


def load_ldscores(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in read_table(path):
        snp = (row.get("snp") or row.get("variant_id") or row.get("rsid") or "").strip()
        if not snp:
            continue
        raw = row.get("l2") or row.get("ldscore") or row.get("ld_score") or ""
        try:
            out[snp] = float(raw)
        except ValueError:
            continue
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--annot", default=None)
    parser.add_argument("--ldscores", default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ref_l2 = load_ldscores(Path(args.ldscores)) if args.ldscores else {}

    snps = []
    zs = []
    ns = []
    l2s = []
    for row in read_table(Path(args.sumstats)):
        snp = (row.get("snp") or row.get("variant_id") or row.get("rsid") or "").strip()
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
        try:
            n = float(row.get("n") or row.get("n_samples") or row.get("sample_size") or 10000)
        except ValueError:
            n = 10000.0

        l2 = None
        if snp in ref_l2:
            l2 = ref_l2[snp]
        elif row.get("l2"):
            try:
                l2 = float(row["l2"])
            except ValueError:
                l2 = None
        if l2 is None:
            # proxy when neither reference nor inline l2 is present
            l2 = 1.0

        snps.append(snp)
        zs.append(z)
        ns.append(n)
        l2s.append(max(l2, 1e-6))

    if not snps:
        raise SystemExit("no usable sumstats rows")

    used_ref = bool(ref_l2) and any(s in ref_l2 for s in snps)
    if (not used_ref) and all(v == 1.0 for v in l2s):
        # synthetic LD scores from chi2 ranks so slope is identifiable
        order = sorted(range(len(zs)), key=lambda i: zs[i] * zs[i])
        for rank, i in enumerate(order, start=1):
            l2s[i] = 1.0 + rank / len(zs)

    chi2 = [z * z for z in zs]
    # E[chi2] ≈ 1 + (N h2 / M) l2  → regress chi2 ~ l2, h2 ≈ slope * M / mean(N)
    _intercept, slope, slope_se = ols(l2s, chi2)
    m = float(len(snps))
    n_bar = sum(ns) / len(ns)
    h2 = max(0.0, slope * m / n_bar)
    h2_se = abs(slope_se * m / n_bar)

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    h2_path = out_prefix.with_suffix(".h2.tsv")
    with h2_path.open("w") as fh:
        fh.write("trait\tnsnp\tn_bar\th2\th2_se\tslope\tldscores\n")
        fh.write(
            f"trait1\t{int(m)}\t{n_bar:.6g}\t{h2:.8g}\t{h2_se:.8g}\t{slope:.8g}\t"
            f"{'ref' if used_ref else 'inline'}\n"
        )

    part_path = out_prefix.with_suffix(".part.tsv")
    with part_path.open("w") as fh:
        fh.write("annotation\tnsnp\th2\tenrichment\n")
        fh.write(f"base\t{int(m)}\t{h2:.8g}\t1\n")
        if args.annot:
            reserved = {"snp", "variant_id", "rsid", "id", "z", "beta", "se", "n", "n_samples", "l2"}
            annot_rows = read_table(Path(args.annot))
            if annot_rows:
                cols = [c for c in annot_rows[0] if c not in reserved]
                by_snp = {
                    (r.get("snp") or r.get("variant_id") or r.get("rsid") or "").strip(): r for r in annot_rows
                }
                for col in cols:
                    idx = [i for i, s in enumerate(snps) if s in by_snp]
                    if len(idx) < 3:
                        continue
                    try:
                        mask = []
                        for i in idx:
                            mask.append(float(by_snp[snps[i]].get(col) or 0.0) != 0.0)
                    except ValueError:
                        continue
                    if not any(mask) or all(mask):
                        continue
                    x = [l2s[i] for i, keep in zip(idx, mask) if keep]
                    y = [chi2[i] for i, keep in zip(idx, mask) if keep]
                    if len(x) < 3:
                        continue
                    _b0, s_ann, _se = ols(x, y)
                    h2_ann = max(0.0, s_ann * len(x) / n_bar)
                    prop = len(x) / m
                    enrich = (h2_ann / h2 / prop) if h2 > 0 and prop > 0 else float("nan")
                    fh.write(f"{col}\t{len(x)}\t{h2_ann:.8g}\t{enrich:.8g}\n")

    out_prefix.with_suffix(".ldsc.log").write_text(
        f"nsnp\t{int(m)}\nh2\t{h2:.8g}\nldscores\t{'ref' if used_ref else 'inline'}\n"
    )


if __name__ == "__main__":
    main()
