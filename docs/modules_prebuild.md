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

## P1 — 协变量与 sQTL 表型准备

| 建议模块名                        | 说明                                         |
| --------------------------------- | -------------------------------------------- |
| `peer`                            | 隐因子协变量估计（eQTL 常用）                |
| `leafcutter/normalize`            | intron ratio 过滤 + 分位数归一化（sQTL）     |
| `leafcutter/annotate`             | intron→基因注释，供 phenotype group          |
| `leafcutter/differentialsplicing` | 官方有开放 PR，尚未合入；可先 local 或等上游 |

> 已安装官方：`regtools/junctionsextract`、`leafcutter/clusterregtools`（聚类上游可用）。

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

| 工具           | Git 策略                               | 仓库                                | 发布                                       |
| -------------- | -------------------------------------- | ----------------------------------- | ------------------------------------------ |
| GEMMA / TASSEL | 官方 bioconda + biocontainers          | —                                   | 直接引用                                   |
| OmiGA          | fork 官方仓 + Release 备份官方 tarball | https://github.com/SiYangming/OmiGA | Conda `YangmingSi`；Quay `bioinfortools`   |
| rMVP           | fork 官方仓，Dockerfile 在 fork        | https://github.com/SiYangming/rMVP  | Quay `bioinfortools`；conda 用 conda-forge |
| EMMAX          | 无官方 GitHub → 自建备份仓             | https://github.com/SiYangming/emmax | Conda `YangmingSi`；Quay `bioinfortools`   |

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
