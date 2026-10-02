# SiYangming/variant2qtl: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - [unreleased]

### `Added`

- Installed nf-core modules/subworkflows for GWAS/QTL scaffolding (plink/plink2, regenie, gcta, VCF annotate/phase/impute, etc.).
- Local engines: GEMMA, EMMAX, TASSEL, rMVP, OmiGA (cis/GWAS) with packaging pins.
- Subworkflows (default OFF): `genotype_ingest_harmonize`, `genotype_qc`, `genotype_to_gwas_formats`, `gwas_benchmark_parallel`, `molqtl_map_omiga`.
- Adapters `phenocovar_adapt` / `assoc_standardize`; OmiGA cis mini testdata under `assets/testdata/omiga_cis_mini/`.
- Docs: `docs/modules_prebuild.md`, `docs/subworkflows_inventory.md`.
- Profile `test_gwas` (`conf/test_gwas.config`) + pipeline nf-test `tests/gwas.nf.test` for QC → GEMMA/EMMAX → standardize; usage docs for genotype/GWAS/OmiGA flags.
- PLINK_RECODE stub emits tped/tfam/vcf.gz when `ext.args` request transpose/vcf so GWAS stub runs can exercise EMMAX.
- Synced pipeline template to nf-core/tools v4.1.0.
- Optional `--genotype_input` cohort samplesheet (`assets/schema_genotype_input.json`): VCF ingest and/or bed+bim+fam + phenotype/covariates; `test_gwas` now exercises ingest→QC→GWAS via samplesheet.
- `--genotype_input` `molqtl_phenotype` / `molqtl_covariates` columns drive OmiGA cis; profile `test_omiga` + nf-test for samplesheet → QC → cis.

### `Fixed`

### `Dependencies`

### `Deprecated`
