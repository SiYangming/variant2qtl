# variant2qtl — Subworkflow 盘点

> nf-core remote 合计 **126** 个 subworkflow。本流程**不跑** RNA 定量（STAR/Salmon 等），下列“无关”类仅作排除说明。

当前仓库已装：`utils_nextflow_pipeline`、`utils_nfcore_pipeline`、`utils_nfschema_plugin`，以及 local `utils_nfcore_variant2qtl_pipeline`。

---

## 一、官方现成、可直接安装复用

### A. 强烈推荐（变异注释 / 亲缘 / 相位填补）

| Subworkflow                        | 作用                               | 对应 variant2qtl 阶段 |
| ---------------------------------- | ---------------------------------- | --------------------- |
| `vcf_annotate_ensemblvep`          | VEP 注释 + bgzip/tabix             | 功能注释              |
| `vcf_annotate_snpeff`              | snpEff 注释                        | 功能注释              |
| `vcf_annotate_ensemblvep_snpeff`   | VEP/snpEff scatter-gather 并行注释 | 高通量注释（优先）    |
| `vcf_filter_bcftools_ensemblvep`   | bcftools + filter_vep              | 注释后过滤            |
| `vcf_gather_bcftools`              | 多分片 VCF concat（可排序）        | scatter 后汇总        |
| `vcf_extract_relate_somalier`      | Somalier extract + relate          | 样本亲缘 / 家系 QC    |
| `vcf_phase_shapeit5`               | SHAPEIT5 相位                      | 基因型相位            |
| `vcf_impute_beagle5`               | Beagle5 填补                       | 基因型填补            |
| `vcf_impute_minimac4`              | Minimac4 填补                      | 基因型填补            |
| `vcf_impute_glimpse`               | GLIMPSE 填补                       | 低深度填补            |
| `bam_vcf_impute_glimpse2`          | GLIMPSE2（BAM/VCF）                | 可选低深度路径        |
| `cache_download_ensemblvep_snpeff` | 下载 VEP/snpEff cache              | 注释前置              |
| `utils_annotation_cache`           | 注释缓存工具                       | 注释前置              |
| `utils_references`                 | 参考基因组资源准备                 | 参考准备              |
| `fasta_bgzip_index_dict_samtools`  | fasta bgzip + fai + dict           | 参考索引              |
| `fasta_index_dna`                  | DNA 参考索引                       | 参考索引              |
| `bed_scatter_bedtools`             | BED 区间分片                       | 按区间并行            |

安装示例：

```bash
nf-core subworkflows install vcf_annotate_ensemblvep_snpeff
nf-core subworkflows install vcf_extract_relate_somalier
nf-core subworkflows install vcf_phase_shapeit5
nf-core subworkflows install vcf_impute_minimac4
```

### B. 可选 / 边缘相关

| Subworkflow                              | 说明                                          |
| ---------------------------------------- | --------------------------------------------- |
| `bam_impute_quilt` / `quilt2` / `stitch` | BAM 级填补（本流程以 VCF 多源输入为主，次要） |
| `deepvariant`                            | 从 BAM 重新 call（用户已有变异时通常不需要）  |
| `bam_variant_calling_*`                  | 体细胞/种系 calling 链（非 QTL 主路径）       |

### C. 明确不需要（表达定量 / RNA / 其它组学）

`fastq_align_*`、`quantify_*`、`bam_qc_rnaseq`、`bam_rseqc`、`bam_stringtie_merge`、`differential_functional_enrichment`、`abundance_differential_filter`、`dia_proteomics_analysis` 等——与“表型矩阵上游已备好”的设计冲突，**不装**。

---

## 二、官方没有、需要本地构建（`subworkflows/local/`）

按流水线阶段排列（依赖 [`docs/modules_prebuild.md`](modules_prebuild.md) 中的 local modules）。

### P0 — 核心业务（必须自建）

| 建议名                        | 职责                                                                                | 依赖模块（多为 local）                                                                               |
| ----------------------------- | ----------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `genotype_ingest_harmonize`   | 多源变异入口（SNP/Indel）→ 标准化 VCF/PLINK；可选 liftover、bcftools norm、样本子集 | `picard/liftovervcf`、`bcftools/{index,norm,view}`、`plink/vcf`；**已深化**（SV/STR 仍在 P3 编排）   |
| `genotype_qc`                 | HWE / missing / MAF / 杂合度 / 亲缘异常 / PCA                                       | `plink2/{filter,het,remove}`、`plink/genome`、local het/relatedness/pca；**已深化**（extras 默认关） |
| `genotype_to_analysis_format` | VCF ↔ BED ↔ 可选 BGEN，供 OmiGA/tensorQTL                                           | **已实现**（bed passthrough + 可选 `plink2/vcf2bgen`）                                               |
| `molqtl_map_omiga`            | **OmiGA** cis（经 analysis-format）                                                 | **已实现**；可 `--omiga_cis_use_qc_bed` 复用 QC bed；pin `1.8.17`                                    |
| `molqtl_map_tensorqtl`        | tensorQTL 备用/对照引擎                                                             | local `tensorqtl`                                                                                    |
| `qtl_postprocess`             | 结果合并、FDR/q-value、按染色体汇总、导出标准表                                     | 轻量 R/Python local                                                                                  |

### P1 — 协变量与 sQTL 表型（表型矩阵侧）

| 建议名              | 职责                                                                     |
| ------------------- | ------------------------------------------------------------------------ |
| `covariate_peer`    | PEER 隐因子 → 协变量矩阵                                                 |
| `sqtl_leafcutter`   | junction → cluster → normalize → annotate（phenotype + phenotype_group） |
| `phenotype_prepare` | 通用表型 QC：样本交集、缺失过滤、quantile/INV、BED 化                    |

### P2 — 精细定位 / 共定位 / 富集

| 建议名           | 职责                                            |
| ---------------- | ----------------------------------------------- |
| `qtl_finemap`    | SuSiE / FINEMAP / CAVIAR / DAP-G                |
| `qtl_coloc`      | coloc / hyprcoloc（QTL–GWAS）                   |
| `qtl_meta_mashr` | 多组织 mashr / METAL                            |
| `qtl_enrichment` | TORUS / 富集（可与 OmiGA 自带富集二选一或串联） |

### P2b — GWAS 基准并行（已实现）

| 建议名                      | 职责                                       | 状态                                                                                |
| --------------------------- | ------------------------------------------ | ----------------------------------------------------------------------------------- |
| `gwas_benchmark_parallel`   | 并行 GEMMA / TASSEL / rMVP / EMMAX / OmiGA | **已实现**；`params.run_gwas_benchmark` 默认关；`--gwas_benchmark_engines` 可选引擎 |
| `genotype_qc`               | MAF/HWE/geno + 可选 het/relatedness/PCA    | **已深化**；`params.run_genotype_qc`；`--genotype_qc_run_{het,relatedness,pca}`     |
| `genotype_ingest_harmonize` | VCF→PLINK + 可选 liftover/norm/samples     | **已深化**；`params.run_genotype_ingest`；chain/fasta/samples 参数见 schema         |

包装与 pin 见 [`docs/modules_prebuild.md`](modules_prebuild.md) P1b。

### P3 — 编排层（pipeline 级）

| 建议名                                           | 职责                                                   |
| ------------------------------------------------ | ------------------------------------------------------ |
| `variant2qtl_snp_indel`                          | SNP/Indel 主路径：ingest → QC → annotate → map         |
| `variant2qtl_sv`                                 | SV 路径（manta/delly/smoove 已装 module，缺编排）      |
| `variant2qtl_str`                                | STR 路径（hipstr/gangstr/EH/trgt 已装 module，缺编排） |
| `variant2qtl_eqtl` / `_sqtl` / `_pqtl` / `_gqtl` | 按 QTL 模态挂表型 + 引擎                               |

---

## 三、结论（一句话）

- **现成可复用**：主要在 **VCF 注释、亲缘、相位/填补、参考与缓存**（约 15+ 个 subworkflow），没有现成的 “QTL mapping” 或 “PLINK QC 全流程” subworkflow。
- **必须自建**：多源变异整合、基因型 QC（plink1/2）、**OmiGA/tensorQTL 映射**、PEER、LeafCutter sQTL、fine-map/coloc，以及按变异类型/QTL 模态的编排层。

## 建议下一步：P1（PEER / LeafCutter / phenotype_prepare）或 P3 编排层（SV/STR）。

P0 SNP/Indel 核心路径（ingest → QC → analysis-format → OmiGA cis / GWAS benchmark）已可用。

## 安装状态

已安装到 `subworkflows/nf-core/`（共 20 个，含模板自带 utils）：

- `vcf_annotate_ensemblvep` / `vcf_annotate_snpeff` / `vcf_annotate_ensemblvep_snpeff`
- `vcf_filter_bcftools_ensemblvep` / `vcf_gather_bcftools`
- `vcf_extract_relate_somalier`
- `vcf_phase_shapeit5`
- `vcf_impute_beagle5` / `vcf_impute_minimac4` / `vcf_impute_glimpse` / `bam_vcf_impute_glimpse2`
- `cache_download_ensemblvep_snpeff` / `utils_annotation_cache` / `utils_references`
- `fasta_bgzip_index_dict_samtools` / `fasta_index_dna` / `bed_scatter_bedtools`
- （模板）`utils_nextflow_pipeline` / `utils_nfcore_pipeline` / `utils_nfschema_plugin`
