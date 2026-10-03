# tensorQTL GEUVADIS official example (local cache)

Fetched by `scripts/fetch_tensorqtl_geuvadis_testdata.sh` from
[broadinstitute/tensorqtl example/data](https://github.com/broadinstitute/tensorqtl/tree/master/example/data).

| File | Role |
| --- | --- |
| `geno.{bed,bim,fam}` | PLINK 1 bed converted from official PLINK2 pgen |
| `phenotype.bed.gz` | Subset of official expression BED (default 5 chr18 genes) |
| `covariates.txt` | Official covariates (covariate × sample) |
| `raw/` | Original downloaded pgen/pvar/psam + full expression |

Large binaries are gitignored. Re-run the fetch script before tagged real nf-tests.
