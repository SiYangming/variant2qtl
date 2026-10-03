# variant2qtl — 预构建模块名单（官方 nf-core 缺失）

> 表达定量软件（STAR/Salmon/RSEM/Kallisto/featureCounts/DESeq2 等）**不纳入**本流程；表型矩阵由上游准备好后直接输入。

官网检索基准：nf-core/modules remote（约 2000 个）。下列工具当前均无对应官方 module，需在 `modules/local/` 自行封装。

## P0 — QTL 核心引擎

| 建议模块名           | 说明                                                                                                                                            | 参考                                   |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| `omiga`              | **OmiGA**：超高效 molQTL（cis / independent-cis / interaction-cis / trans）与 GWAS；支持 LM/LMM、加性/非加性、GRM、自动协变量；适合复杂亲缘群体 | https://omiga.bio/#/ ；Nat Commun 2026 |
| `tensorqtl`          | GTEx 系主流 eQTL/sQTL/pQTL 引擎（GPU/CPU）                                                                                                      | —                                      |
| `qtltools`           | cis/trans QTL、quan、correct、pca 等                                                                                                            | —                                      |
| `fastqtl`            | 经典 cis-eQTL                                                                                                                                   | —                                      |
| `matrixeqtl`         | R Matrix eQTL                                                                                                                                   | —                                      |
| `saige-qtl` / `apex` | 大规模 / 混合模型 QTL                                                                                                                           | —                                      |
| `mmqtl`              | 多样本 / 多组织 QTL                                                                                                                             | —                                      |

### P0 — OmiGA 模块状态（已建 / 未接入 workflow）

| 模块路径                   | 状态                                                                 | Conda / 容器 pin                                                  |
| -------------------------- | -------------------------------------------------------------------- | ----------------------------------------------------------------- |
| `modules/local/omiga/cis`  | **已建，已挂入** `molqtl_map_omiga`（`params.run_omiga_cis` 默认关） | `YangmingSi::omiga=1.8.17` / `quay.io/bioinfortools/omiga:1.8.17` |
| `modules/local/omiga/gwas` | **已建，已挂入** `gwas_benchmark_parallel`                           | 同上                                                              |

- 官方无 bioconda/biocontainers；自建包装见 fork https://github.com/SiYangming/OmiGA （[PR #1](https://github.com/SiYangming/OmiGA/pull/1)，[release 1.8.17](https://github.com/SiYangming/OmiGA/releases/tag/1.8.17)）
- conda version 与 Docker tag 均为 `1.8.17`
- 上游二进制：`OmiGA-build-v1.8.17-x86_64.tar.xz`（官网 latest.json）

#### OmiGA cis — mini testdata & smoke

Official nf-core testdata has **no** molecular phenotype matrix. Mini files live in [`assets/testdata/omiga_cis_mini/`](../assets/testdata/omiga_cis_mini/) (~30 KB total):

| File                 | Role                                                                                             |
| -------------------- | ------------------------------------------------------------------------------------------------ |
| `geno.{bed,bim,fam}` | Copy of nf-core `plink_simulated` (200 samples, chr1); bim alleles recoded `D/d`→`A/T` for OmiGA |
| `phenotype.bed.gz`   | FastQTL-style BED (`#chr start end pheno_id` + matching IIDs; 3 fake genes)                      |
| `covariates.txt`     | OmiGA default covariate×sample matrix (Sex/Age/PC1/PC2)                                          |

Regenerate: `bash scripts/make_omiga_cis_mini_testdata.sh`

**nf-test (module):**

```bash
# stub (CI-friendly)
nf-test test modules/local/omiga/cis/tests/main.nf.test --tag stub --profile test,docker

# real OmiGA run (needs docker image; tags gwas_real / omiga_cis_real)
nf-test test modules/local/omiga/cis/tests/main.nf.test --tag omiga_cis_real --profile test,docker
```

**Pipeline smoke** (booleans via params-file; default `run_omiga_cis=false` does not touch the main path):

```bash
# example params JSON (paths absolute or relative to launch dir)
# { "input": "...", "outdir": "...", "run_omiga_cis": true,
#   "omiga_cis_bed": "assets/testdata/omiga_cis_mini/geno.bed",
#   "omiga_cis_bim": "assets/testdata/omiga_cis_mini/geno.bim",
#   "omiga_cis_fam": "assets/testdata/omiga_cis_mini/geno.fam",
#   "omiga_cis_phenotype": "assets/testdata/omiga_cis_mini/phenotype.bed.gz",
#   "omiga_cis_covariates": "assets/testdata/omiga_cis_mini/covariates.txt" }
nextflow run . -profile test,docker -params-file params_omiga_cis.json
```

### P0 — tensorQTL 模块状态（已建 / 已接入 workflow）

| 模块路径                      | 状态                                                                         | Conda / 容器 pin                                                                           |
| ----------------------------- | ---------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `modules/local/tensorqtl/cis` | **已建，已挂入** `molqtl_map_tensorqtl`（`params.run_tensorqtl_cis` 默认关） | PyPI `tensorqtl==1.0.10` via conda `pip`（无 bioconda/biocontainers）；stub / wave / conda |

- Official bioconda/biocontainers package **absent** — pin strategy mirrors other local engines: document version in `environment.yml`; CI uses `-stub`; real runs prefer `-profile conda` or `wave`.
- Stub phenotype: FastQTL BED (reuse `assets/testdata/omiga_cis_mini/`). Covariates: **sample × covariate** (`covariates_tensorqtl.txt`).
- **Official real testdata:** GEUVADIS chr18 example from [broadinstitute/tensorqtl `example/data`](https://github.com/broadinstitute/tensorqtl/tree/master/example/data). Fetch (gitignored binaries): `bash scripts/fetch_tensorqtl_geuvadis_testdata.sh` → `assets/testdata/tensorqtl_geuvadis/`.
- GPU: full official notebook targets GPU/~50GB; tagged real nf-test uses phenotype subset + `--mode cis_nominal` for CPU smoke. Pass CUDA flags via `ext.args` when available.

**nf-test / smoke:**

```bash
nf-test test modules/local/tensorqtl/cis/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_tensorqtl,docker -stub --outdir results_test_tensorqtl

# Real (conda; downloads ~80MB GEUVADIS example on first run; ignored by default CI)
bash scripts/run_real_nf_tests.sh tensorqtl
```

### P0 — QTLtools 模块状态（已建 / 已接入 workflow）

| 模块路径                     | 状态                                                                       | Conda / 容器 pin                                                      |
| ---------------------------- | -------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `modules/local/qtltools/cis` | **已建，已挂入** `molqtl_map_qtltools`（`params.run_qtltools_cis` 默认关） | `YangmingSi::qtltools=1.3.1` / `quay.io/bioinfortools/qtltools:1.3.1` |

- Converts PLINK bed → bgzipped VCF inside the process (`plink2 --export vcf bgz`).
- Phenotype: FastQTL BED (reuse `assets/testdata/omiga_cis_mini/`).
- Packaging fork: https://github.com/SiYangming/qtltools (`PACKAGING.md`). Pins: `YangmingSi::qtltools=1.3.1` and `quay.io/bioinfortools/qtltools:1.3.1`.
- CI: `-stub`; real runs prefer `-profile conda` or `wave`.

```bash
nf-test test modules/local/qtltools/cis/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_qtltools,docker -stub --outdir results_test_qtltools
```

## P1 — 协变量与 sQTL 表型准备

| 建议模块名                        | 说明                                                                                   |
| --------------------------------- | -------------------------------------------------------------------------------------- |
| `peer`                            | **已建**：隐因子协变量（`params.run_peer` 默认关）                                     |
| `leafcutter/normalize`            | intron ratio 过滤 + 分位数归一化（sQTL）                                               |
| `leafcutter/annotate`             | intron→基因注释，供 phenotype group                                                    |
| `leafcutter/differentialsplicing` | 官方有开放 PR，尚未合入；可先 local 或等上游                                           |
| `phenotype_prepare`               | **已建**：样本交集 / 缺失 / INV / FastQTL BED（`params.run_phenotype_prepare` 默认关） |

> 已安装官方：`regtools/junctionsextract`、`leafcutter/clusterregtools`（聚类上游可用）。

### P1 — phenotype_prepare（已建 / 已接入 workflow）

| 模块路径                                | 状态                                                                          | Conda / 容器 pin              |
| --------------------------------------- | ----------------------------------------------------------------------------- | ----------------------------- |
| `modules/local/utils/phenotype_prepare` | **已建，已挂入** `phenotype_prepare`（`params.run_phenotype_prepare` 默认关） | biocontainers `python:3.9--1` |

Mini testdata: [`assets/testdata/phenotype_prepare_mini/`](../assets/testdata/phenotype_prepare_mini/).

```bash
nf-test test modules/local/utils/phenotype_prepare/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_pheno,docker -stub --outdir results_test_pheno
```

### P1 — PEER（已建 / 已接入 workflow）

| 模块路径                     | 状态                                                          | Conda / 容器 pin                                                      |
| ---------------------------- | ------------------------------------------------------------- | --------------------------------------------------------------------- |
| `modules/local/peer/factors` | **已建，已挂入** `covariate_peer`（`params.run_peer` 默认关） | `bioconda::r-peer=1.3` / `quay.io/biocontainers/peer:1.3--h503566f_1` |

- Input: FastQTL BED (`--peer_phenotype` or phenotype_prepare output); optional known covariates.
- Emits sample×factor (`*.peer.tensor.tsv`) for tensorQTL and factor×sample (`*.peer.omiga.tsv`) for OmiGA/QTLtools.
- Mini testdata: [`assets/testdata/peer_mini/`](../assets/testdata/peer_mini/).

```bash
nf-test test modules/local/peer/factors/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_peer,docker -stub --outdir results_test_peer
```

## P1b — GWAS/QTL 基准对照（已建模块 / 已挂入 `gwas_benchmark_parallel`）

> 与 OmiGA / regenie / gcta 并列作工具基准；conda 版本号与 Docker tag 同一字符串（自建包强制）。
> 开关：`params.run_gwas_benchmark`（**默认 false**）。启用：`--run_gwas_benchmark` + `--gwas_benchmark_{bed,bim,fam,phenotype}`（可选 vcf/covariates）。

| 模块路径                          | 状态             | Conda / 容器 pin                                                              |
| --------------------------------- | ---------------- | ----------------------------------------------------------------------------- |
| `modules/local/gemma/relatedness` | **已建，已挂入** | `bioconda::gemma=0.98.5` / `quay.io/biocontainers/gemma:0.98.5--h38cc83e_1`   |
| `modules/local/gemma/lmm`         | **已建，已挂入** | 同上                                                                          |
| `modules/local/tassel/mlm`        | **已建，已挂入** | `bioconda::tassel=5.2.89` / `quay.io/biocontainers/tassel:5.2.89--hdfd78af_0` |
| `modules/local/rmvp/gwas`         | 已建，已挂入     | `conda-forge::r-rmvp=1.4.6` / `quay.io/bioinfortools/rmvp:1.4.6`              |
| `modules/local/emmax/kin`         | 已建，已挂入     | `YangmingSi::emmax=0.0.20100307` / `quay.io/bioinfortools/emmax:0.0.20100307` |
| `modules/local/emmax/assoc`       | 已建，已挂入     | 同上                                                                          |

**Locked GWAS smoke testdata** (base = `params.modules_testdata_base_path` in [`tests/nextflow.config`](../tests/nextflow.config)):

| File             | Relative path                                                                                           |
| ---------------- | ------------------------------------------------------------------------------------------------------- |
| bed/bim/fam      | `genomics/homo_sapiens/popgen/plink_simulated.{bed,bim,fam}`                                            |
| VCF              | `genomics/homo_sapiens/popgen/plink_simulated.vcf.gz`                                                   |
| quantitative phe | `genomics/homo_sapiens/popgen/plink_simulated_quantitative_phenoname.phe` (`FID IID QuantitativeTrait`) |
| covariates       | `genomics/homo_sapiens/popgen/plink_simulated_covariates.txt` (`FID IID Sex Age PC1 PC2`)               |

Adapters: `modules/local/utils/phenocovar_adapt`, `modules/local/utils/assoc_standardize`; format fan-out: `subworkflows/local/genotype_to_gwas_formats` (EMMAX tped via `plink --recode 12 transpose`); benchmark orchestration: `subworkflows/local/gwas_benchmark_parallel`.

**Notes:** `plink_simulated` alleles are `D`/`d`; EMMAX requires numeric alleles (`12 transpose`). `EMMAX_ASSOC` real runs need amd64 (amd64 container on arm64 may hit Illegal instruction) — use `--gwas_benchmark_engines gemma,omiga,...` without `emmax`, or `--tag gwas_real_amd64` for the assoc nf-test. Minimal QC: `--run_genotype_qc` (PLINK2 MAF/HWE/geno) feeds filtered bed into the benchmark when both flags are on. Minimal ingest: `--run_genotype_ingest` + `--genotype_ingest_vcf` (VCF→bed via `PLINK_VCF`); when ingest is ON it supplies bed instead of `--gwas_benchmark_bed` (chain: ingest → optional QC → benchmark).

### 包装仓与发布账号

| 工具           | Git 策略                               | 仓库                                   | 发布                                            |
| -------------- | -------------------------------------- | -------------------------------------- | ----------------------------------------------- |
| GEMMA / TASSEL | 官方 bioconda + biocontainers          | —                                      | 直接引用                                        |
| OmiGA          | fork 官方仓 + Release 备份官方 tarball | https://github.com/SiYangming/OmiGA    | Conda `YangmingSi`；Quay `bioinfortools`        |
| rMVP           | fork 官方仓，Dockerfile 在 fork        | https://github.com/SiYangming/rMVP     | Quay `bioinfortools`；conda 用 conda-forge      |
| EMMAX          | 无官方 GitHub → 自建备份仓             | https://github.com/SiYangming/emmax    | Conda `YangmingSi`；Quay `bioinfortools`        |
| QTLtools       | fork 官方仓 + 官方 source tarball      | https://github.com/SiYangming/qtltools | Conda `YangmingSi`；Quay `bioinfortools`        |
| PEER           | 官方 bioconda + biocontainers          | —                                      | `bioconda::r-peer=1.3` / `peer:1.3--h503566f_1` |

构建环境：本地 `conda_build`。GitHub 操作：`gh` CLI。

## P2 — 精细定位 / 共定位 / 下游

| 建议模块名                     | 说明               |
| ------------------------------ | ------------------ |
| `susie` / `susier`             | SuSiE fine-mapping |
| `finemap` / `caviar` / `dap-g` | 其它 fine-mapping  |
| `coloc` / `hyprcoloc`          | QTL–GWAS 共定位    |
| `mashr`                        | 多组织 QTL 收缩    |
| `metal`                        | 跨队列 meta        |
| `smr`                          | SMR / HEIDI        |
| `torus`                        | QTL 富集先验       |

### P2 — SuSiE 模块状态（已建 / 已接入 workflow）

| 模块路径                      | 状态                                                                      | Conda / 容器 pin                                                      |
| ----------------------------- | ------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `modules/local/susie/finemap` | **已建，已挂入** `qtl_finemap_susie`（`params.run_finemap_susie` 默认关） | `conda-forge::r-susier=0.14.2` / biocontainers `r-base:4.3.1` + conda |

- Input: association/QTL sumstats (`variant_id`/`snp`, `beta`+`se` and/or `z`); optional LD matrix. Without LD, identity matrix is used (smoke / limited interpretation).
- Stub mini testdata: [`assets/testdata/finemap_susie_mini/`](../assets/testdata/finemap_susie_mini/).
- **Official-style real testdata:** [`assets/testdata/finemap_susie_official/`](../assets/testdata/finemap_susie_official/) from `susieR::N3finemapping` or the susieR vignette simulation recipe (`scripts/make_finemap_susie_official_testdata.R` / `.py`).
- Can consume OmiGA/tensorQTL/QTLtools `cis_qtl` outputs when `--finemap_susie_sumstats` is unset and those branches ran.

### P2 — coloc 模块状态（已建 / 已接入 workflow）

| 模块路径                  | 状态                                                      | Conda / 容器 pin                                                    |
| ------------------------- | --------------------------------------------------------- | ------------------------------------------------------------------- |
| `modules/local/coloc/abf` | **已建，已挂入** `qtl_coloc`（`params.run_coloc` 默认关） | `conda-forge::r-coloc=5.2.3` / biocontainers `r-base:4.3.1` + conda |

- Input: QTL + GWAS sumstats with overlapping `snp`/`variant_id`, `beta`, `se` (optional `maf`, `n`). Mini files: [`assets/testdata/coloc_mini/`](../assets/testdata/coloc_mini/).
- Can consume OmiGA/tensorQTL/QTLtools `cis_qtl` plus GWAS standardized tables when `--coloc_*_sumstats` unset.
- **hyprcoloc** remains a follow-up.

```bash
nf-test test modules/local/coloc/abf/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_coloc,docker -stub --outdir results_test_coloc
```

### P2 — cis postprocess（已建 / 已接入 workflow）

| 模块路径                              | 状态                                                                          | Conda / 容器 pin              |
| ------------------------------------- | ----------------------------------------------------------------------------- | ----------------------------- |
| `modules/local/utils/qtl_postprocess` | **已建，已挂入** `qtl_postprocess_cis`（`params.run_qtl_postprocess` 默认关） | biocontainers `python:3.9--1` |

```bash
nf-test test modules/local/utils/qtl_postprocess/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_postprocess,docker -stub --outdir results_test_postprocess
```

Other fine-mappers (FINEMAP/CAVIAR/DAP-G) remain follow-ups.

```bash
nf-test test modules/local/susie/finemap/tests/main.nf.test --tag stub --profile test,docker
nextflow run . -profile test_finemap,docker -stub --outdir results_test_finemap

# Real (conda; uses committed official-style sumstats+LD; ignored by default CI)
bash scripts/run_real_nf_tests.sh susie
```

## P3 — 可选互补

| 建议模块名                               | 说明                               | 现有替代                              |
| ---------------------------------------- | ---------------------------------- | ------------------------------------- |
| `crossmap`                               | 坐标升版本                         | `ucsc/liftover`、`picard/liftovervcf` |
| `annovar`                                | 变异注释                           | `ensemblvep`、`snpeff`                |
| `king`                                   | 亲缘/族系 QC                       | `somalier/relate`、`plink/genome`     |
| `bolt-lmm`                               | 大规模混合模型 GWAS                | `regenie`、`omiga`、见 P1b GEMMA      |
| `ldsc`                                   | LD score / 遗传力                  | 仅有 `gcta/calculateldscores`         |
| `predixcan` / `metaxcan` / TWAS-`fusion` | TWAS（勿与融合基因 fusion\* 混淆） | —                                     |

## 已从 nf-core 安装（摘要）

见仓库根目录 `modules_installed_list.txt`（当前约 95 个），含完整 **plink** + **plink2**、VCF/SV/STR、VEP/snpEff、相位填补、regenie/gcta 等。不含任何 RNA 定量模块。
