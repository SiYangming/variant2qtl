#!/usr/bin/env python3
"""SMR & HEIDI from overlapping QTL and GWAS summary statistics (Zhu et al. 2016)."""
import argparse
import csv
import math
import sys


def chi2_sf(x, df=1):
    """Survival function of chi^2 (df=1 uses erfc)."""
    if x <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))
    # Wilson-Hilferty-ish fallback for small df=k: not used
    return math.erfc(math.sqrt(x / 2.0))


def parse_float(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def read_ss(path):
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t, ")
        reader = csv.DictReader(fh, dialect=dialect)
        if reader.fieldnames is None:
            sys.exit(f"empty header: {path}")
        fields = {f.strip().lstrip("#").lower(): f for f in reader.fieldnames}

        def col(*names):
            for n in names:
                if n in fields:
                    return fields[n]
            return None

        snp_c = col("snp", "variant_id", "rsid", "id")
        gene_c = col("gene", "gene_id", "phenotype_id", "probe", "pheno_id")
        beta_c = col("beta", "b", "effect", "slope")
        se_c = col("se", "stderr", "std_error", "slope_se")
        p_c = col("p", "pval", "p_value", "pvalue")
        if not snp_c or not beta_c or not se_c:
            sys.exit(f"{path} needs snp, beta, se")
        out = {}
        for row in reader:
            snp = str(row[snp_c]).strip()
            beta = parse_float(row[beta_c])
            se = parse_float(row[se_c])
            if not snp or beta is None or se is None or se <= 0:
                continue
            gene = str(row[gene_c]).strip() if gene_c else "NA"
            p = parse_float(row[p_c]) if p_c else None
            out.setdefault(gene, {})[snp] = {"beta": beta, "se": se, "p": p}
    return out


def smr_one(bzx, sezx, bzy, sezy):
    # bxy = bzy / bzx
    bxy = bzy / bzx
    var = (sezy ** 2) / (bzx ** 2) + (bzy ** 2) * (sezx ** 2) / (bzx ** 4)
    if var <= 0:
        return None
    t = (bxy ** 2) / var
    p = chi2_sf(t, 1)
    return bxy, math.sqrt(var), t, p


def heidi(estimates):
    """Cochran Q on SMR effect sizes vs top eQTL."""
    if len(estimates) < 2:
        return "NA", "NA"
    top = estimates[0]
    q = 0.0
    df = 0
    for bxy, se, _t, _p in estimates[1:]:
        if se <= 0:
            continue
        q += ((bxy - top[0]) ** 2) / (se ** 2)
        df += 1
    if df < 1:
        return "NA", "NA"
    # chi2_sf for df>1: Wilson–Hilferty
    if df == 1:
        p = chi2_sf(q, 1)
    else:
        z = ((q / df) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df))) / math.sqrt(2.0 / (9.0 * df))
        p = 0.5 * math.erfc(z / math.sqrt(2.0))
    return f"{q:.6g}", f"{p:.6g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qtl", required=True)
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    qtl = read_ss(args.qtl)
    gwas_all = read_ss(args.gwas)
    # GWAS often a single dummy gene key
    gwas_snps = {}
    for _g, snps in gwas_all.items():
        gwas_snps.update(snps)

    rows_out = []
    for gene, snps in qtl.items():
        shared = []
        for snp, qx in snps.items():
            if snp not in gwas_snps:
                continue
            gy = gwas_snps[snp]
            est = smr_one(qx["beta"], qx["se"], gy["beta"], gy["se"])
            if est is None:
                continue
            p_eqtl = qx["p"]
            if p_eqtl is None:
                z = abs(qx["beta"] / qx["se"])
                p_eqtl = chi2_sf(z * z, 1)
            shared.append((snp, p_eqtl, est))
        if not shared:
            continue
        shared.sort(key=lambda x: x[1])
        top_snp, _p, top_est = shared[0]
        heidi_q, heidi_p = heidi([e for _s, _p, e in shared])
        bxy, se, t, psmr = top_est
        rows_out.append((gene, top_snp, bxy, se, t, psmr, heidi_q, heidi_p, len(shared)))

    with open(args.out + ".smr.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "top_snp", "b_smr", "se_smr", "t_smr", "p_smr", "heidi_Q", "p_heidi", "n_snps"])
        for rec in rows_out:
            gene, snp, bxy, se, t, psmr, hq, hp, n = rec
            w.writerow([gene, snp, f"{bxy:.6g}", f"{se:.6g}", f"{t:.6g}", f"{psmr:.6g}", hq, hp, n])


if __name__ == "__main__":
    main()
