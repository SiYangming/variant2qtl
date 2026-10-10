# SiYangming/variant2qtl: Output

## Introduction

Directories below are relative to `--outdir`. Most optional branches are **off by default**; only enabled flags create their folders. FastQC/MultiQC always run from `--input`.

## Pipeline overview

- [FastQC](#fastqc) — raw read QC
- [MultiQC](#multiqc) — aggregated QC
- [Pipeline information](#pipeline-information) — Nextflow reports and software versions
- [Genotype / QTL](#genotype--qtl) — ingest, QC, GWAS, cis mapping, coloc, TWAS, LDSC
- [VCF prep](#vcf-prep) — annotate, filter, phase, impute, relate, FASTA index, caches

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) collates pipeline QC from supported tools.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Nextflow reports: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`, `pipeline_dag.*`
  - `variant2qtl_software_mqc_versions.yml`
  - `params.json`

</details>

### Genotype / QTL

<details markdown="1">
<summary>Output files</summary>

- `genotype_ingest/` — VCF→PLINK ingest (includes SV/STR feed when `--sv_feed_ingest` / `--str_feed_ingest`; biallelic filter drops complex alleles)
- `genotype_qc/` — PLINK2 filter, optional het/relatedness/PCA, `plink_recode_vcf` after QC, Somalier outlier remove
- `genotype_analysis/` — optional BGEN (`plink2_vcf2bgen`) when `--genotype_to_bgen` (uses QC-rebuilt VCF when QC ran)
- `gwas_benchmark/` — GEMMA / EMMAX / TASSEL / rMVP / OmiGA GWAS
- `phenotype_prepare/` — FastQTL BED from expression matrix
- `covariate_peer/` — PEER factors
- `sqtl_leafcutter/` — intron ratios / sQTL BED; optional `regtools_junctionsextract` / `leafcutter_clusterregtools`
- `molqtl_omiga/` — `omiga_cis` / `omiga_trans` / `omiga_independent`
- `molqtl_matrixeqtl/` — MatrixEQTL cis
- `molqtl_tensorqtl/` / `molqtl_qtltools/` — cis QTL
- `qtl_postprocess/` — BH q-values
- `qtl_finemap_susie/` / `qtl_finemap_extra/` — PIPs
- `qtl_coloc/` / `qtl_hyprcoloc/` / `qtl_smr/` / `qtl_mashr/` / `qtl_metal/` / `qtl_torus/`
- `qtl_twas/` — FUSION / PredictDB
- `qtl_ldsc/` — h2 / partitioned enrichment
- `variant_sv/` / `variant_str/` — SV / STR callers; optional `survivor_merge` / `trtools_mergestr` when 2+ engines

</details>

### VCF prep

<details markdown="1">
<summary>Output files</summary>

- `annotation_cache/` — VEP / snpEff cache download
- `reference_fasta/` — bgzip FASTA, faidx, dict
- `variant_annotate/` — optional `vcfanno`; VEP / snpEff (scatter helpers when `annotate_sites_per_chunk`); optional `bcftools_view` / `ensemblvep_filtervep` when `annotate_filter`
- `variant_phase/` — SHAPEIT5
- `variant_impute/` — Beagle5 / Minimac4 / GLIMPSE
- `genotype_ingest/` — also receives `--vcf_prep_feed_ingest` prepared VCFs
- `scatter/` — optional `--phase_scatter_bed` / `--impute_scatter_bed` (BED split, region strings, concat/index)
- `variant_impute_bam/` — GLIMPSE2 chunk / phase / ligate
- `variant_relate/` — Somalier extract / relate

</details>
