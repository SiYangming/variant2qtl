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
- Local **LeafCutter sQTL prepare** (`modules/local/leafcutter/prepare`, `sqtl_leafcutter`) behind `params.run_sqtl_leafcutter` (default OFF); intron ratios + INV + gene annotation → FastQTL BED and `phenotype_group`; profile `test_sqtl` + stub nf-test.
- Local **HyPrColoc-style** clustering (`modules/local/coloc/hyprcoloc`) behind `params.run_hyprcoloc` (default OFF); pairwise Wakefield ABF (official `r-hyprcoloc` is not on bioconda); profile `test_hyprcoloc` + stub nf-test.
- Local **SMR/HEIDI** (`modules/local/smr/heidi`) behind `params.run_smr` (default OFF); Zhu et al. 2016 on overlapping sumstats (official SMR binary is not on bioconda); profile `test_smr` + stub nf-test.
- Local **mashr** (`modules/local/mashr/fit`, `qtl_mashr`) behind `params.run_mashr` (default OFF); conda-forge `r-mashr=0.2.79` / r-base stub container; profile `test_mashr` + stub nf-test.
- Local **METAL-style IVW** (`modules/local/metal/ivw`, `qtl_metal`) behind `params.run_metal` (default OFF); official METAL is not on bioconda; profile `test_metal` + stub nf-test.
- Local **TORUS-style enrichment** (`modules/local/torus/enrich`, `qtl_torus`) behind `params.run_torus` (default OFF); official TORUS is not on bioconda; profile `test_torus` + stub nf-test.
- Local **TWAS** (`modules/local/twas/fusion`, `qtl_twas`) behind `params.run_twas` (default OFF); FUSION/PrediXcan-style independent-SNP gene score; profile `test_twas` + stub nf-test.
- Extra fine-mapping (`modules/local/finemap/extra`) behind `params.run_finemap_extra` (default OFF); FINEMAP/CAVIAR/DAP-G style 1-causal PIPs; profile `test_finemap_extra` + stub nf-test.
- P3 **SNP/Indel umbrella** `params.run_snp_indel` (default OFF) runs ingest (if needed) + QC + OmiGA cis; profile `test_snp_indel`.
- P3 **SV** (`variant_sv`) behind `params.run_sv` (default OFF); nf-core smoove/manta/delly; profile `test_sv`.
- P3 **STR** (`variant_str`) behind `params.run_str` (default OFF); nf-core ExpansionHunter/GangSTR/HipSTR/TRGT; profile `test_str`.
- QTL modality umbrellas `params.run_eqtl` / `run_sqtl` / `run_pqtl` (default OFF) + `qtl_modality_engine` (omiga/tensorqtl/qtltools); profiles `test_eqtl` / `test_modality_sqtl` / `test_pqtl`.
- Local **LDSC-style h2** (`modules/local/ldsc/h2`, `qtl_ldsc`) behind `params.run_ldsc` (default OFF); optional partitioned enrichment; profile `test_ldsc` + stub nf-test.
- TWAS deepen: MetaXcan `GENE/RSID/WEIGHT` columns and optional `--twas_ld` for FUSION-style `w'Rw` variance.
- gQTL modality umbrella `params.run_gqtl` (default OFF): ingest + QC + GWAS benchmark; profile `test_gqtl`.
- PredictDB / PrediXcan SQLite import (`modules/local/twas/predictdb`) when `--twas_weights` ends with `.db`; profile `test_twas_predictdb`.
- LDSC reference-panel LD scores via `--ldsc_ldscores` (snp + l2), merged by SNP id.
- Cap nf-test process resources in `tests/nextflow.config` so stub GWAS/gQTL paths stay under local memory limits.
- Wire nf-core **annotate / phase / impute**: `run_annotate` (VEP/snpEff), `run_phase` (SHAPEIT5), `run_impute` (beagle5/minimac4/glimpse); profiles `test_annotate` / `test_phase` / `test_impute`.
- Somalier `run_relate`, `run_vcf_prep` umbrella, optional `genotype_qc_use_somalier` hook; CI-ignored real smokes for LDSC/PredictDB panels and VEP/SHAPEIT5/Beagle.
- Sequential `run_vcf_prep` chain (annotate VCF feeds phase, phased VCF feeds impute) instead of OR-ing the three independent flags.
- `run_cache` (Ensembl VEP / snpEff download), `run_impute_bam` (GLIMPSE2), `run_fasta_index` (bgzip + samtools faidx/dict); profiles `test_cache` / `test_impute_bam` / `test_fasta_index`.
- Example params YAML (`assets/params_eqtl.yml`, `params_vcf_prep.yml`, `params_cache_annotate.yml`) and output-directory docs.
- Feed `run_fasta_index` / `run_cache` into relate, VEP, SV/STR, and GLIMPSE2 when per-step FASTA/cache params are unset.
- Optional `annotate_filter` via `vcf_filter_bcftools_ensemblvep`; `run_impute_bam` uses `bam_vcf_impute_glimpse2` (chunk/phase/ligate).
- `--sv_feed_ingest` / `--str_feed_ingest` mix called VCFs into `genotype_ingest` with a biallelic `bcftools view` (complex SV/STR alleles are dropped before PLINK).
- Optional `--phase_scatter_bed` / `--impute_scatter_bed` (`bed_scatter_bedtools` + `bed_to_region` + `vcf_gather_bcftools`).
- README shortest command uses `-params-file assets/params_eqtl.yml`. **No 0.0.1 tag.**

### `Fixed`

### `Dependencies`

### `Deprecated`
