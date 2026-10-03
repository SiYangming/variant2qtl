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
- Local **tensorQTL** cis (`modules/local/tensorqtl/cis`, `molqtl_map_tensorqtl`) behind `params.run_tensorqtl_cis` (default OFF); PyPI pin `tensorqtl==1.0.10` (no bioconda); profile `test_tensorqtl` + stub nf-test; reuses `omiga_cis_mini` geno/pheno + `covariates_tensorqtl.txt`.
- Local **SuSiE** fine-mapping (`modules/local/susie/finemap`, `qtl_finemap_susie`) behind `params.run_finemap_susie` (default OFF); conda-forge `r-susier=0.14.2`; mini sumstats under `assets/testdata/finemap_susie_mini/`; profile `test_finemap` + stub nf-test.
- Optional **real** module nf-tests (ignored by default CI): tensorQTL on official GEUVADIS example (`scripts/fetch_tensorqtl_geuvadis_testdata.sh`, `--tag tensorqtl_real`); SuSiE on susieR vignette/N3finemapping-derived sumstats (`assets/testdata/finemap_susie_official/`, `--tag susie_real`).
- Local **coloc.abf** (`modules/local/coloc/abf`, `qtl_coloc`) behind `params.run_coloc` (default OFF); conda-forge `r-coloc=5.2.3`; mini overlapping sumstats under `assets/testdata/coloc_mini/`; profile `test_coloc` + stub nf-test. hyprcoloc not included.
- Local **QTLtools** cis (`modules/local/qtltools/cis`, `molqtl_map_qtltools`) behind `params.run_qtltools_cis` (default OFF); conda pin `YangmingSi::qtltools=1.3.1` and Docker `quay.io/bioinfortools/qtltools:1.3.1` (fork https://github.com/SiYangming/qtltools; no official bioconda/biocontainers); profile `test_qtltools` + stub nf-test.
- Local **cis-QTL postprocess** (`modules/local/utils/qtl_postprocess`, `qtl_postprocess_cis`) behind `params.run_qtl_postprocess` (default OFF); BH q-values + standard table; profile `test_postprocess` + stub nf-test.
- Local **phenotype_prepare** (`modules/local/utils/phenotype_prepare`) behind `params.run_phenotype_prepare` (default OFF); sample intersect, missingness filter, INV/quantile, FastQTL BED; profile `test_pheno` + stub nf-test.
- Local **PEER** covariates (`modules/local/peer/factors`, `covariate_peer`) behind `params.run_peer` (default OFF); bioconda `r-peer=1.3` / `quay.io/biocontainers/peer:1.3--h503566f_1`; feeds OmiGA/QTLtools and tensorQTL `--cov` when those branches have no explicit covariates; profile `test_peer` + stub nf-test.

### `Fixed`

### `Dependencies`

### `Deprecated`
