# SiYangming/variant2qtl: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- NOTE: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Genotype ingest, QC, GWAS benchmark, molQTL, and fine-mapping

These branches are **off by default**. Enable them with boolean params; they do not replace the FASTQC/MultiQC template path (still driven by `--input`).

### Preferred input: `--genotype_input`

Cohort CSV (schema: `assets/schema_genotype_input.json`). One row per analysis set:

| Column              | Required | Notes                                                                                               |
| ------------------- | -------- | --------------------------------------------------------------------------------------------------- |
| `id`                | yes      | Cohort / meta id                                                                                    |
| `vcf`               | xor bed  | When set without bed → ingest (VCF→PLINK)                                                           |
| `bed`/`bim`/`fam`   | xor vcf  | All three required together when skipping ingest                                                    |
| `phenotype`         | no\*     | PLINK-style `FID IID Trait` (needed for `--run_gwas_benchmark`)                                     |
| `covariates`        | no       | Optional GWAS `FID IID cov...`                                                                      |
| `molqtl_phenotype`  | no\*     | FastQTL-style BED/OPF (needed for `--run_omiga_cis` / `--run_tensorqtl_cis` / `--run_qtltools_cis`) |
| `molqtl_covariates` | no       | Optional molQTL covariates (OmiGA: covariate×sample; tensorQTL: sample×covar)                       |

\*Provide the phenotype column matching the branch you enable (GWAS and/or molQTL).

Examples: `assets/genotype_samplesheet.csv` (VCF ingest + GWAS), `assets/genotype_samplesheet_bed.csv` (bed-only GWAS), `assets/genotype_samplesheet_omiga.csv` (bed + OmiGA cis), `assets/genotype_samplesheet_tensorqtl.csv` (bed + tensorQTL cis), `assets/genotype_samplesheet_qtltools.csv` (bed + QTLtools cis).

When `--genotype_input` is set, it **overrides** scattered `--gwas_benchmark_bed/bim/fam/vcf/phenotype/covariates`, `--genotype_ingest_vcf`, and (for molQTL) `--omiga_cis_*` / `--tensorqtl_cis_*` file params when the matching samplesheet columns are present.

### Feature flags

| Param                 | Default | What it does                                                                                                                               |
| --------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `run_genotype_ingest` | `false` | Legacy single-VCF ingest via `--genotype_ingest_vcf` when **not** using `--genotype_input`. VCF-only samplesheet rows always ingest.       |
| `run_genotype_qc`     | `false` | MAF/HWE/geno filter (`PLINK2_FILTER`); optional het / relatedness / PCA via `genotype_qc_run_*`. Filtered bed feeds GWAS when both are on. |
| `run_gwas_benchmark`  | `false` | Parallel engines via `gwas_benchmark_parallel` (pheno adapt → formats → engines → optional `assoc_standardize`).                           |
| `run_omiga_cis`       | `false` | OmiGA cis-molQTL (`molqtl_map_omiga`). Use dedicated bed/bim/fam or reuse QC bed.                                                          |
| `run_tensorqtl_cis`   | `false` | tensorQTL cis-molQTL (`molqtl_map_tensorqtl`). Prefer samplesheet molQTL columns or `--tensorqtl_use_qc_bed`.                              |
| `run_finemap_susie`   | `false` | SuSiE fine-mapping (`qtl_finemap_susie`) from `--finemap_susie_sumstats` or cis QTL outputs from OmiGA/tensorQTL/QTLtools.                 |
| `run_coloc`           | `false` | QTL–GWAS `coloc.abf` (`qtl_coloc`) from `--coloc_*_sumstats` or cis QTL + GWAS standardized tables. hyprcoloc not wired.                   |
| `run_qtltools_cis`    | `false` | QTLtools cis-molQTL (`molqtl_map_qtltools`). Prefer samplesheet molQTL columns or `--qtltools_use_qc_bed`.                                 |
| `run_qtl_postprocess` | `false` | Harmonise cis tables from any engine and add BH q-values (`qtl_postprocess_cis`).                                                          |

### Related file / option params

**GWAS benchmark** (legacy single-cohort; ignored for genotype files when `--genotype_input` is set)

- `gwas_benchmark_id` — meta id (default `gwas_benchmark`)
- `gwas_benchmark_bed` / `_bim` / `_fam` — PLINK inputs when ingest is off
- `gwas_benchmark_vcf` — optional VCF (else rebuilt from bed)
- `gwas_benchmark_phenotype` — PLINK-style `FID IID Trait` (header required)
- `gwas_benchmark_covariates` — optional `FID IID cov...`
- `gwas_benchmark_engines` — comma list, e.g. `gemma,emmax` (default includes tassel/rmvp/omiga)
- `gwas_benchmark_standardize` — run `assoc_standardize` on engine outputs (default `true`)

**Genotype ingest**

- `genotype_ingest_vcf` — required when ingest is on **without** `--genotype_input`
- `genotype_ingest_args`, `genotype_ingest_fasta`, `genotype_ingest_norm_args`
- `genotype_ingest_samples`, `genotype_ingest_view_args`
- `genotype_ingest_liftover_chain` / `_fasta` / `_dict`

**Genotype QC**

- `genotype_qc_args` — PLINK2 filter args
- `genotype_qc_run_het` / `genotype_qc_het_sd`
- `genotype_qc_run_relatedness` / `genotype_qc_pi_hat`
- `genotype_qc_run_pca` / `genotype_qc_pca_n`

**OmiGA cis**

- Prefer samplesheet `molqtl_phenotype` / `molqtl_covariates` with `--genotype_input`
- Legacy: `omiga_cis_id`, `omiga_cis_bed` / `_bim` / `_fam`, `omiga_cis_vcf`, `omiga_cis_phenotype`, `omiga_cis_covariates`
- `omiga_cis_use_qc_bed` — reuse shared QC bed from the genotype QC branch instead of samplesheet/legacy bed

**tensorQTL cis**

- Prefer samplesheet `molqtl_phenotype` / `molqtl_covariates` (sample × covariate orientation; see `covariates_tensorqtl.txt`)
- Legacy: `tensorqtl_cis_id`, `tensorqtl_cis_bed` / `_bim` / `_fam`, `tensorqtl_cis_phenotype`, `tensorqtl_cis_covariates`
- `tensorqtl_use_qc_bed` — reuse shared QC bed
- Packaging: not on bioconda; pin `tensorqtl==1.0.10` via conda `pip` / prefer `conda` or `wave` for real runs; `-stub` for CI

**SuSiE fine-mapping**

- `finemap_susie_sumstats` — TSV/CSV with `variant_id`/`snp`, `beta`+`se` and/or `z`
- `finemap_susie_ld` — optional square LD matrix (identity used when omitted)
- When sumstats unset, reuses OmiGA/tensorQTL/QTLtools `cis_qtl` outputs if those branches ran

**QTLtools cis**

- Prefer samplesheet `molqtl_phenotype` / `molqtl_covariates`
- Legacy: `qtltools_cis_id`, `qtltools_cis_bed` / `_bim` / `_fam`, `qtltools_cis_phenotype`, `qtltools_cis_covariates`
- `qtltools_use_qc_bed` — reuse shared QC bed
- Converts PLINK bed → VCF inside the process. Official bioconda/biocontainers package is absent; conda pin `YangmingSi::qtltools=1.3.1`; Docker `quay.io/bioinfortools/qtltools:1.3.1`; `-stub` for CI

**QTL–GWAS coloc**

- `coloc_qtl_sumstats` / `coloc_gwas_sumstats` — `snp`/`variant_id`, `beta`, `se` (optional `maf`, `n`)
- When unset, reuses cis QTL outputs and GWAS `assoc_standardize` tables if those branches ran
- `coloc.abf` pairwise; HyPrColoc-style clustering via `--run_hyprcoloc --hyprcoloc_sumstats` (long `snp/trait/beta/se`; official r-hyprcoloc is not on bioconda)

**SMR / HEIDI**

- `smr_qtl_sumstats` / `smr_gwas_sumstats` — else cis QTL + GWAS standardised tables
- Summary-stat SMR and HEIDI (Zhu et al. 2016); official SMR binary is not on bioconda

**LeafCutter sQTL**

- `sqtl_counts` — LeafCutter perind `count/total`; `sqtl_genes` optional gene BED
- Writes FastQTL BED + `phenotype_group`; feeds cis engines like phenotype_prepare
- Cluster BAM junctions with nf-core `leafcutter/clusterregtools` upstream if needed

**cis-QTL postprocess**

- `qtl_postprocess_input` — optional cis table; else uses OmiGA/tensorQTL/QTLtools `cis_qtl`
- Writes `*.cis_std.tsv` with BH q-values (`gene variant_id chr pos p q beta se engine`)

**Phenotype prepare**

- `phenotype_matrix` — genes×samples TSV/CSV or FastQTL BED
- `phenotype_gene_bed` — optional `chr start end gene_id`
- `phenotype_samples` — optional keep list or FAM
- `phenotype_transform` — `invnorm` (default) / `quantile` / `none`; `phenotype_max_missing` (default 0.2)
- Writes `*.phenotype.bed.gz` for OmiGA/tensorQTL/QTLtools when those branches have no explicit phenotype

**PEER covariates**

- `peer_phenotype` — FastQTL BED (else uses phenotype_prepare output)
- `peer_covariates` — optional known covariates (either orientation)
- `peer_nk` — hidden factor count (default 10)
- Writes sample×factor (`*.peer.tensor.tsv`) for tensorQTL and factor×sample (`*.peer.omiga.tsv`) for OmiGA/QTLtools when those branches have no explicit covariates
- Pin: `bioconda::r-peer=1.3` / `quay.io/biocontainers/peer:1.3--h503566f_1`

### Recommended combinations

1. **Samplesheet VCF → ingest → QC → GWAS**: `--genotype_input assets/genotype_samplesheet.csv --run_genotype_qc --run_gwas_benchmark --gwas_benchmark_engines gemma,emmax`
2. **Samplesheet bed → QC → GWAS**: `--genotype_input assets/genotype_samplesheet_bed.csv --run_genotype_qc --run_gwas_benchmark`
3. **Legacy params bed → QC → GWAS**: `--run_genotype_qc --run_gwas_benchmark` + `--gwas_benchmark_{bed,bim,fam,phenotype}`
4. **Samplesheet bed → QC → OmiGA cis**: `--genotype_input assets/genotype_samplesheet_omiga.csv --run_genotype_qc --run_omiga_cis --omiga_cis_use_qc_bed`
5. **Legacy QC bed → OmiGA cis**: `--run_omiga_cis --omiga_cis_use_qc_bed --omiga_cis_phenotype ...`
6. **Samplesheet bed → QC → tensorQTL cis**: `--genotype_input assets/genotype_samplesheet_tensorqtl.csv --run_genotype_qc --run_tensorqtl_cis --tensorqtl_use_qc_bed`
7. **SuSiE from sumstats**: `--run_finemap_susie --finemap_susie_sumstats assets/testdata/finemap_susie_mini/sumstats.tsv`
8. **Samplesheet bed → QC → QTLtools cis**: `--genotype_input assets/genotype_samplesheet_qtltools.csv --run_genotype_qc --run_qtltools_cis --qtltools_use_qc_bed`
9. **coloc from mini sumstats**: `--run_coloc --coloc_qtl_sumstats assets/testdata/coloc_mini/qtl.tsv --coloc_gwas_sumstats assets/testdata/coloc_mini/gwas.tsv`
10. **cis postprocess from a table**: `--run_qtl_postprocess --qtl_postprocess_input assets/testdata/qtl_postprocess_mini/cis_qtl.tsv`
11. **Phenotype matrix → FastQTL BED**: `--run_phenotype_prepare --phenotype_matrix assets/testdata/phenotype_prepare_mini/expr.tsv --phenotype_gene_bed assets/testdata/phenotype_prepare_mini/genes.bed --phenotype_samples assets/testdata/phenotype_prepare_mini/samples.txt`
12. **PEER factors from a BED**: `--run_peer --peer_phenotype assets/testdata/peer_mini/phenotype.bed --peer_nk 2`
13. **LeafCutter sQTL BED**: `--run_sqtl_leafcutter --sqtl_counts assets/testdata/sqtl_leafcutter_mini/perind.counts.tsv --sqtl_genes assets/testdata/sqtl_leafcutter_mini/genes.bed`
14. **HyPrColoc-style clustering**: `--run_hyprcoloc --hyprcoloc_sumstats assets/testdata/hyprcoloc_mini/traits.tsv`
15. **SMR/HEIDI from mini sumstats**: `--run_smr --smr_qtl_sumstats assets/testdata/coloc_mini/qtl.tsv --smr_gwas_sumstats assets/testdata/coloc_mini/gwas.tsv`
16. **Phenotype matrix → FastQTL BED**: `--run_phenotype_prepare --phenotype_matrix assets/testdata/phenotype_prepare_mini/expr.tsv --phenotype_gene_bed assets/testdata/phenotype_prepare_mini/genes.bed --phenotype_samples assets/testdata/phenotype_prepare_mini/samples.txt`
17. **PEER factors from a BED**: `--run_peer --peer_phenotype assets/testdata/peer_mini/phenotype.bed --peer_nk 2`

Notes: `plink_simulated` alleles are `D`/`d`; EMMAX paths recode to numeric (`12 transpose`). Real EMMAX binaries are amd64-oriented — on arm64 prefer engines without `emmax`, or use `-stub`. OmiGA/tensorQTL mini testdata uses A/T-recoded alleles under `assets/testdata/omiga_cis_mini/`.

### Example: `test_gwas` profile

Smoke-tests samplesheet **VCF ingest → QC → gemma/emmax → standardize** (`assets/genotype_samplesheet.csv`). Prefer `-stub` for CI:

```bash
nextflow run . -profile test_gwas,docker -stub --outdir results_test_gwas
```

Without stub (amd64 docker recommended if including EMMAX):

```bash
nextflow run . -profile test_gwas,docker --outdir results_test_gwas
```

### Example: `test_omiga` profile

Smoke-tests samplesheet **bed → QC → OmiGA cis** (`assets/genotype_samplesheet_omiga.csv`):

```bash
nextflow run . -profile test_omiga,docker -stub --outdir results_test_omiga
```

### Example: `test_tensorqtl` profile

Smoke-tests samplesheet **bed → QC → tensorQTL cis** (`assets/genotype_samplesheet_tensorqtl.csv`):

```bash
nextflow run . -profile test_tensorqtl,docker -stub --outdir results_test_tensorqtl
```

Optional **real** module test on official GEUVADIS example (downloads ~80MB; not in default CI):

```bash
bash scripts/run_real_nf_tests.sh tensorqtl
```

### Example: `test_finemap` profile

Smoke-tests **SuSiE** on mini sumstats:

```bash
nextflow run . -profile test_finemap,docker -stub --outdir results_test_finemap
```

Optional **real** module test on susieR vignette / `N3finemapping`-style sumstats (not in default CI):

```bash
bash scripts/run_real_nf_tests.sh susie
```

### Example: `test_qtltools` profile

Smoke-tests samplesheet **bed → QC → QTLtools cis** (`assets/genotype_samplesheet_qtltools.csv`):

```bash
nextflow run . -profile test_qtltools,docker -stub --outdir results_test_qtltools
```

### Example: `test_coloc` profile

Smoke-tests **coloc.abf** on mini overlapping sumstats:

```bash
nextflow run . -profile test_coloc,docker -stub --outdir results_test_coloc
```

### Example: `test_postprocess` profile

Smoke-tests **BH q-value** standardisation of a cis table:

```bash
nextflow run . -profile test_postprocess,docker -stub --outdir results_test_postprocess
```

### Example: `test_pheno` profile

Smoke-tests **phenotype_prepare** (intersect / missing / INV → FastQTL BED):

```bash
nextflow run . -profile test_pheno,docker -stub --outdir results_test_pheno
```

### Example: `test_peer` profile

Smoke-tests **PEER** factors from a mini FastQTL BED:

```bash
nextflow run . -profile test_peer,docker -stub --outdir results_test_peer
```

### Example: `test_sqtl` profile

```bash
nextflow run . -profile test_sqtl,docker -stub --outdir results_test_sqtl
```

### Example: `test_hyprcoloc` profile

```bash
nextflow run . -profile test_hyprcoloc,docker -stub --outdir results_test_hyprcoloc
```

### Example: `test_smr` profile

```bash
nextflow run . -profile test_smr,docker -stub --outdir results_test_smr
```

The existing `-profile test,docker` FASTQC-only path is unchanged.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run SiYangming/variant2qtl --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run SiYangming/variant2qtl -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull SiYangming/variant2qtl
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [SiYangming/variant2qtl releases page](https://github.com/SiYangming/variant2qtl/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_gwas`
- `test_omiga`
- `test_tensorqtl`
- `test_finemap`
- `test_qtltools`
- `test_coloc`
- `test_postprocess`
- `test_pheno`
- `test_peer`
- `test_sqtl`
- `test_hyprcoloc`
- `test_smr`
  - Genotype QC + GWAS / molQTL / SuSiE smoke profiles; prefer with `-stub` in CI; see [Genotype ingest, QC, GWAS benchmark, molQTL, and fine-mapping](#genotype-ingest-qc-gwas-benchmark-molqtl-and-fine-mapping)
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
