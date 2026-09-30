# nf-core/variant2qtl: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - 2026-09-29

### `Added`

- Installed nf-core modules/subworkflows for GWAS/QTL scaffolding (plink/plink2, regenie, gcta, VCF annotate/phase/impute, etc.).
- Local engines: GEMMA, EMMAX, TASSEL, rMVP, OmiGA (cis/GWAS) with packaging pins.
- Subworkflows (default OFF): `genotype_ingest_harmonize`, `genotype_qc`, `genotype_to_gwas_formats`, `gwas_benchmark_parallel`, `molqtl_map_omiga`.
- Adapters `phenocovar_adapt` / `assoc_standardize`; OmiGA cis mini testdata under `assets/testdata/omiga_cis_mini/`.
- Docs: `docs/modules_prebuild.md`, `docs/subworkflows_inventory.md`.

### `Fixed`

### `Dependencies`

### `Deprecated`

## v0.0.1 - 2026-09-24

Initial release of **variant2qtl**, created with the [nf-core](https://nf-co.re/) template.

An end-to-end, multi-entry Nextflow pipeline for comprehensive QTL (eQTL/sQTL/pQTL) analysis, integrating diverse genetic variants (SNPs, Indels, STRs, and SVs) with automated functional annotation.

一款基于 Nextflow DSL2 构建的、支持多源变异（SNP/Indel/SV/STR）输入的全自动 QTL（eQTL/sQTL/pQTL/gQTL）高通量关联分析流程。

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
