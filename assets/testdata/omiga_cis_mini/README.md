# OmiGA cis mini testdata

Tiny molecular phenotype + covariates + genotype derived from nf-core `plink_simulated`
(200 samples, chr1 SNPs at positions 1–220).

| File                 | Description                                                                     |
| -------------------- | ------------------------------------------------------------------------------- |
| `geno.{bed,bim,fam}` | Copy of `plink_simulated`; bim alleles recoded `D/d`→`A/T` for OmiGA            |
| `phenotype.bed.gz`   | FastQTL-style BED: `#chr start end pheno_id` + per-sample values (3 fake genes) |
| `covariates.txt`     | OmiGA default orientation (covariate × sample); Sex/Age/PC1/PC2                 |

Regenerate:

```bash
bash scripts/make_omiga_cis_mini_testdata.sh
```
