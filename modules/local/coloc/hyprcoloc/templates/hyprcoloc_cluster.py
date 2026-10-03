#!/usr/bin/env python3
"""HyPrColoc-style multi-trait clustering via pairwise Wakefield ABF (no official bioconda hyprcoloc)."""
import argparse
import csv
import math
import sys
from collections import defaultdict


def ppf_lbf(beta, se, W=0.15 ** 2):
    if se is None or se <= 0 or beta is None:
        return None
    v = se * se
    z = beta / se
    # Wakefield approximate Bayes factor (log)
    return 0.5 * math.log(v / (v + W)) + 0.5 * z * z * (W / (v + W))


def coloc_h4(rows_a, rows_b):
    """Approximate PP.H4 for two traits on shared SNPs (product of ABFs)."""
    by_a = {r["snp"]: r for r in rows_a}
    by_b = {r["snp"]: r for r in rows_b}
    snps = sorted(set(by_a) & set(by_b))
    if len(snps) < 2:
        return 0.0, snps[:1]
    lbf1 = []
    lbf2 = []
    for s in snps:
        a = ppf_lbf(by_a[s]["beta"], by_a[s]["se"])
        b = ppf_lbf(by_b[s]["beta"], by_b[s]["se"])
        if a is None or b is None:
            continue
        lbf1.append(a)
        lbf2.append(b)
    if len(lbf1) < 2:
        return 0.0, snps[:1]
    # unnormalized H4 ~ sum exp(lbf1+lbf2); H0=1
    m = max(x + y for x, y in zip(lbf1, lbf2))
    h4 = sum(math.exp(x + y - m) for x, y in zip(lbf1, lbf2))
    h4 *= math.exp(m)
    h0 = 1.0
    h1 = sum(math.exp(x) for x in lbf1)
    h2 = sum(math.exp(y) for y in lbf2)
    tot = h0 + h1 + h2 + h4
    pp4 = h4 / tot if tot > 0 else 0.0
    best = snps[max(range(len(lbf1)), key=lambda i: lbf1[i] + lbf2[i])]
    return pp4, [best]


def parse_float(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traits", required=True)
    ap.add_argument("--pp4-threshold", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    by_trait = defaultdict(list)
    with open(args.traits, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t, ")
        reader = csv.DictReader(fh, dialect=dialect)
        if reader.fieldnames is None:
            sys.exit("traits table needs a header")
        fields = {f.strip().lstrip("#").lower(): f for f in reader.fieldnames}

        def col(*names):
            for n in names:
                if n in fields:
                    return fields[n]
            return None

        snp_c = col("snp", "variant_id", "rsid", "id")
        trait_c = col("trait", "phenotype", "name")
        beta_c = col("beta", "b", "effect")
        se_c = col("se", "stderr", "std_error")
        if not all([snp_c, trait_c, beta_c, se_c]):
            sys.exit("need snp, trait, beta, se columns")
        for row in reader:
            snp = str(row[snp_c]).strip()
            trait = str(row[trait_c]).strip()
            beta = parse_float(row[beta_c])
            se = parse_float(row[se_c])
            if not snp or not trait or beta is None or se is None:
                continue
            by_trait[trait].append({"snp": snp, "beta": beta, "se": se})

    traits = sorted(by_trait)
    if len(traits) < 2:
        sys.exit("need >=2 traits")

    parent = {t: t for t in traits}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    pairs = []
    for i, a in enumerate(traits):
        for b in traits[i + 1 :]:
            pp4, snps = coloc_h4(by_trait[a], by_trait[b])
            pairs.append((a, b, pp4, snps[0] if snps else ""))
            if pp4 >= args.pp4_threshold:
                pa, pb = find(a), find(b)
                if pa != pb:
                    parent[pb] = pa

    clusters = defaultdict(list)
    for t in traits:
        clusters[find(t)].append(t)

    with open(args.out + ".clusters.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster", "traits", "n_traits"])
        for i, (_root, members) in enumerate(sorted(clusters.items()), start=1):
            w.writerow([f"C{i}", ",".join(sorted(members)), len(members)])
    with open(args.out + ".pairs.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait1", "trait2", "PP.H4", "candidate_snp"])
        for a, b, pp4, snp in pairs:
            w.writerow([a, b, f"{pp4:.6g}", snp])


if __name__ == "__main__":
    main()
