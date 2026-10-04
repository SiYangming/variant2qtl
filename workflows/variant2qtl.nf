/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { samplesheetToList      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variant2qtl_pipeline'
include { GWAS_BENCHMARK_PARALLEL   } from '../subworkflows/local/gwas_benchmark_parallel/main'
include { GENOTYPE_QC               } from '../subworkflows/local/genotype_qc/main'
include { GENOTYPE_INGEST_HARMONIZE } from '../subworkflows/local/genotype_ingest_harmonize/main'
include { MOLQTL_MAP_OMIGA          } from '../subworkflows/local/molqtl_map_omiga/main'
include { MOLQTL_MAP_TENSORQTL      } from '../subworkflows/local/molqtl_map_tensorqtl/main'
include { MOLQTL_MAP_QTLTOOLS       } from '../subworkflows/local/molqtl_map_qtltools/main'
include { QTL_FINEMAP_SUSIE         } from '../subworkflows/local/qtl_finemap_susie/main'
include { QTL_COLOC                 } from '../subworkflows/local/qtl_coloc/main'
include { QTL_POSTPROCESS_CIS       } from '../subworkflows/local/qtl_postprocess_cis/main'
include { PHENOTYPE_PREPARE_SWF     } from '../subworkflows/local/phenotype_prepare/main'
include { COVARIATE_PEER            } from '../subworkflows/local/covariate_peer/main'
include { SQTL_LEAFCUTTER           } from '../subworkflows/local/sqtl_leafcutter/main'
include { QTL_HYPRCOLOC             } from '../subworkflows/local/qtl_hyprcoloc/main'
include { QTL_SMR                   } from '../subworkflows/local/qtl_smr/main'
include { QTL_MASHR                 } from '../subworkflows/local/qtl_mashr/main'
include { QTL_METAL                 } from '../subworkflows/local/qtl_metal/main'
include { QTL_TORUS                 } from '../subworkflows/local/qtl_torus/main'
include { QTL_TWAS                  } from '../subworkflows/local/qtl_twas/main'
include { QTL_LDSC                  } from '../subworkflows/local/qtl_ldsc/main'
include { QTL_FINEMAP_EXTRA         } from '../subworkflows/local/qtl_finemap_extra/main'
include { VARIANT_SV                } from '../subworkflows/local/variant_sv/main'
include { VARIANT_STR               } from '../subworkflows/local/variant_str/main'
include { VARIANT_ANNOTATE          } from '../subworkflows/local/variant_annotate/main'
include { VARIANT_PHASE             } from '../subworkflows/local/variant_phase/main'
include { VARIANT_IMPUTE            } from '../subworkflows/local/variant_impute/main'
include { VARIANT_VCF_PREP          } from '../subworkflows/local/variant_vcf_prep/main'
include { ANNOTATION_CACHE          } from '../subworkflows/local/annotation_cache/main'
include { VARIANT_IMPUTE_BAM        } from '../subworkflows/local/variant_impute_bam/main'
include { REFERENCE_FASTA           } from '../subworkflows/local/reference_fasta/main'
include { VARIANT_RELATE            } from '../subworkflows/local/variant_relate/main'
include { SOMALIER_OUTLIERS         } from '../modules/local/utils/somalier_outliers/main'
include { PLINK2_REMOVE as PLINK2_REMOVE_SOMALIER } from '../modules/nf-core/plink2/remove/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT2QTL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def run_snp_indel = params.run_snp_indel
    def run_eqtl = params.run_eqtl
    def run_sqtl = params.run_sqtl
    def run_pqtl = params.run_pqtl
    def run_gqtl = params.run_gqtl
    def run_modality = run_eqtl || run_sqtl || run_pqtl
    def modality_engine = (params.qtl_modality_engine ?: 'omiga').toString().trim().toLowerCase()
    def run_pheno_prep = params.run_phenotype_prepare || run_eqtl || run_pqtl
    def run_leaf = params.run_sqtl_leafcutter || run_sqtl
    def run_peer_flag = params.run_peer || run_modality
    def run_gwas = params.run_gwas_benchmark || run_gqtl
    def run_ingest = params.run_genotype_ingest || run_snp_indel || run_gqtl
    def run_qc = params.run_genotype_qc || run_snp_indel || run_gqtl || (run_modality && params.genotype_input)
    def run_vcf_prep = params.run_vcf_prep
    def run_annotate_flag = params.run_annotate
    def run_phase_flag = params.run_phase
    def run_impute_flag = params.run_impute
    def vcf_prep_vcf = params.vcf_prep_vcf ?: params.annotate_vcf ?: params.phase_vcf ?: params.impute_vcf
    def vcf_prep_tbi = params.vcf_prep_vcf_tbi ?: params.annotate_vcf_tbi ?: params.phase_vcf_tbi ?: params.impute_vcf_tbi
    def run_omiga = params.run_omiga_cis || run_snp_indel || (run_modality && modality_engine == 'omiga')
    def run_tensor = params.run_tensorqtl_cis || (run_modality && modality_engine == 'tensorqtl')
    def run_qtltools = params.run_qtltools_cis || (run_modality && modality_engine == 'qtltools')
    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })
    // FastQC 4.1.0 publishes versions via topic("versions"), not FASTQC.out.versions

    //
    // Optional genotype ingest → QC → GWAS benchmark (all default OFF).
    // Prefer --genotype_input samplesheet (one row per cohort) over scattered
    // --gwas_benchmark_* / --genotype_ingest_vcf params. Rows with VCF (no bed)
    // are ingested; rows with bed/bim/fam skip ingest. Shared QC / geno bed and
    // molqtl_* columns can feed OmiGA / tensorQTL cis when those flags are set.
    //
    def ch_shared_qc_bed = channel.empty()
    def ch_shared_geno_plink = channel.empty()
    def ch_shared_geno_vcf = channel.empty()
    def ch_molqtl_pheno = channel.empty()
    def ch_molqtl_covar = channel.empty()
    def ch_qtl_cis_for_finemap = channel.empty()
    def ch_gwas_std = channel.empty()
    def ch_prepared_pheno = channel.empty()
    def ch_peer_cov_omiga = channel.empty()
    def ch_peer_cov_tensor = channel.empty()
    def ch_somalier_pairs = channel.empty()
    def ch_downloaded_vep = channel.empty()
    def ch_downloaded_snpeff = channel.empty()
    def ch_fasta_indexed = channel.empty()
    def fasta_index_ok = params.run_fasta_index && params.fasta_index_fasta

    if (run_ingest || run_qc || run_gwas || params.genotype_input) {
        def ch_gwas_plink_raw = channel.empty()
        def ch_gwas_vcf = channel.empty()
        def ch_gwas_pheno = channel.empty()
        def ch_gwas_covar = channel.empty()

        if (params.genotype_input) {
            def ch_geno_rows = channel
                .fromList(samplesheetToList(params.genotype_input, "${projectDir}/assets/schema_genotype_input.json"))
                .map { meta, vcf, bed, bim, fam, phenotype, covariates, molqtl_phenotype, molqtl_covariates ->
                    def paths = [vcf, bed, bim, fam, phenotype, covariates, molqtl_phenotype, molqtl_covariates].collect { pathish ->
                        // Prefer toUriString(): Path.toString() strips remote schemes
                        if (pathish == null || pathish instanceof List || pathish instanceof Collection) {
                            return null
                        }
                        def s = null
                        try {
                            s = pathish.toUriString()?.trim()
                        }
                        catch (Throwable _ignored) {
                        }
                        if (!s) {
                            s = pathish.toString()?.trim()
                        }
                        if (!s || s == 'null' || s == '[]') {
                            return null
                        }
                        if (s.contains('://')) {
                            return file(s, checkIfExists: true)
                        }
                        def candidates = [file(s), file("${projectDir}/${s}")]
                        def assetsIdx = s.indexOf('assets/')
                        if (assetsIdx >= 0) {
                            candidates << file("${projectDir}/${s.substring(assetsIdx)}")
                        }
                        def f = candidates.find { cand -> cand.exists() }
                        if (!f) {
                            error("genotype_input path not found: ${s}")
                        }
                        return f
                    }
                    def vcf_p = paths[0]
                    def bed_p = paths[1]
                    def bim_p = paths[2]
                    def fam_p = paths[3]
                    def phe_p = paths[4]
                    def cov_p = paths[5]
                    def mol_phe_p = paths[6]
                    def mol_cov_p = paths[7]
                    def has_bed = bed_p && bim_p && fam_p
                    def has_vcf = vcf_p != null
                    if (!has_bed && !has_vcf) {
                        error("genotype_input row '${meta.id}': provide vcf and/or bed+bim+fam")
                    }
                    if ((bed_p || bim_p || fam_p) && !has_bed) {
                        error("genotype_input row '${meta.id}': bed, bim, and fam must all be set together")
                    }
                    [meta, vcf_p, bed_p, bim_p, fam_p, phe_p, cov_p, mol_phe_p, mol_cov_p, has_bed, has_vcf]
                }

            def ch_geno_branched = ch_geno_rows.branch { row ->
                from_bed: row[9]
                from_vcf: row[10] && !row[9]
            }

            ch_gwas_pheno = ch_geno_rows
                .filter { row -> row[5] != null }
                .map { meta, _vcf, _bed, _bim, _fam, phe, _cov, _mp, _mc, _hb, _hv -> [meta, phe] }
            ch_gwas_covar = ch_geno_rows
                .filter { row -> row[6] != null }
                .map { meta, _vcf, _bed, _bim, _fam, _phe, cov, _mp, _mc, _hb, _hv -> [meta, cov] }
            ch_molqtl_pheno = ch_geno_rows
                .filter { row -> row[7] != null }
                .map { meta, _vcf, _bed, _bim, _fam, _phe, _cov, mp, _mc, _hb, _hv -> [meta, mp] }
            ch_molqtl_covar = ch_geno_rows
                .filter { row -> row[8] != null }
                .map { meta, _vcf, _bed, _bim, _fam, _phe, _cov, _mp, mc, _hb, _hv -> [meta, mc] }

            ch_gwas_plink_raw = ch_geno_branched.from_bed.map { meta, _vcf, bed, bim, fam, _phe, _cov, _mp, _mc, _hb, _hv ->
                [meta, bed, bim, fam]
            }
            ch_gwas_vcf = ch_geno_branched.from_bed
                .filter { row -> row[1] != null }
                .map { meta, vcf, _bed, _bim, _fam, _phe, _cov, _mp, _mc, _hb, _hv -> [meta, vcf] }
                .mix(
                    ch_geno_branched.from_vcf.map { meta, vcf, _bed, _bim, _fam, _phe, _cov, _mp, _mc, _hb, _hv -> [meta, vcf] }
                )

            // VCF-only rows always go through ingest; empty channel is a no-op.
            GENOTYPE_INGEST_HARMONIZE(
                ch_geno_branched.from_vcf.map { meta, vcf, _bed, _bim, _fam, _phe, _cov, _mp, _mc, _hb, _hv -> [meta, vcf] }
            )
            ch_versions = ch_versions.mix(GENOTYPE_INGEST_HARMONIZE.out.versions)
            ch_gwas_plink_raw = ch_gwas_plink_raw.mix(GENOTYPE_INGEST_HARMONIZE.out.bed)
            ch_gwas_vcf = ch_gwas_vcf.mix(GENOTYPE_INGEST_HARMONIZE.out.vcf)
        } else {
            def gwas_meta = [id: params.gwas_benchmark_id ?: 'gwas_benchmark']

            if (run_ingest) {
                if (!params.genotype_ingest_vcf) {
                    log.warn "run_genotype_ingest/run_snp_indel=true but missing --genotype_ingest_vcf; channels empty."
                } else {
                    def ch_ingest_vcf = channel.of([
                        gwas_meta,
                        file(params.genotype_ingest_vcf, checkIfExists: true)
                    ])
                    GENOTYPE_INGEST_HARMONIZE(ch_ingest_vcf)
                    ch_versions = ch_versions.mix(GENOTYPE_INGEST_HARMONIZE.out.versions)
                    ch_gwas_plink_raw = GENOTYPE_INGEST_HARMONIZE.out.bed
                    ch_gwas_vcf = GENOTYPE_INGEST_HARMONIZE.out.vcf
                }
            } else {
                ch_gwas_plink_raw = (params.gwas_benchmark_bed && params.gwas_benchmark_bim && params.gwas_benchmark_fam)
                    ? channel.of([
                        gwas_meta,
                        file(params.gwas_benchmark_bed, checkIfExists: true),
                        file(params.gwas_benchmark_bim, checkIfExists: true),
                        file(params.gwas_benchmark_fam, checkIfExists: true)
                    ])
                    : channel.empty()

                ch_gwas_vcf = params.gwas_benchmark_vcf
                    ? channel.of([gwas_meta, file(params.gwas_benchmark_vcf, checkIfExists: true)])
                    : channel.empty()

                if ((run_qc || run_gwas) &&
                    (!params.gwas_benchmark_bed || !params.gwas_benchmark_bim || !params.gwas_benchmark_fam)) {
                    log.warn "genotype_qc/gwas_benchmark/run_gqtl enabled but missing --gwas_benchmark_bed/bim/fam; channels empty."
                }
            }

            ch_gwas_pheno = params.gwas_benchmark_phenotype
                ? channel.of([gwas_meta, file(params.gwas_benchmark_phenotype, checkIfExists: true)])
                : channel.empty()

            ch_gwas_covar = params.gwas_benchmark_covariates
                ? channel.of([gwas_meta, file(params.gwas_benchmark_covariates, checkIfExists: true)])
                : channel.empty()
        }

        def ch_gwas_plink = ch_gwas_plink_raw
        if (run_qc) {
            GENOTYPE_QC(ch_gwas_plink_raw)
            ch_versions = ch_versions.mix(GENOTYPE_QC.out.versions)
            ch_gwas_plink = GENOTYPE_QC.out.bed
            ch_shared_qc_bed = GENOTYPE_QC.out.bed
            // QC changes SNP set — force VCF rebuild from filtered bed
            ch_gwas_vcf = channel.empty()
        }

        ch_shared_geno_plink = ch_gwas_plink
        ch_shared_geno_vcf = ch_gwas_vcf

        if (run_gwas) {
            if (!params.genotype_input && !params.gwas_benchmark_phenotype) {
                log.warn "run_gwas_benchmark/run_gqtl=true but missing --gwas_benchmark_phenotype (or --genotype_input)."
            }
            GWAS_BENCHMARK_PARALLEL(
                ch_gwas_plink,
                ch_gwas_vcf,
                ch_gwas_pheno,
                ch_gwas_covar
            )
            ch_versions = ch_versions.mix(GWAS_BENCHMARK_PARALLEL.out.versions)
            ch_gwas_std = GWAS_BENCHMARK_PARALLEL.out.standardized
        }
    }

    //
    // Optional phenotype_prepare → FastQTL BED (default OFF)
    //
    if (run_pheno_prep) {
        def pheno_meta = [id: params.phenotype_prepare_id ?: 'phenotype_prepare']
        if (!params.phenotype_matrix) {
            log.warn "run_phenotype_prepare/run_eqtl/run_pqtl=true but missing --phenotype_matrix."
        } else {
            def ch_pheno_matrix = channel.of([pheno_meta, file(params.phenotype_matrix, checkIfExists: true)])
            def ch_pheno_genes = params.phenotype_gene_bed
                ? channel.of([pheno_meta, file(params.phenotype_gene_bed, checkIfExists: true)])
                : channel.empty()
            def ch_pheno_samples = params.phenotype_samples
                ? channel.of([pheno_meta, file(params.phenotype_samples, checkIfExists: true)])
                : channel.empty()
            PHENOTYPE_PREPARE_SWF(ch_pheno_matrix, ch_pheno_genes, ch_pheno_samples)
            ch_versions = ch_versions.mix(PHENOTYPE_PREPARE_SWF.out.versions)
            ch_prepared_pheno = PHENOTYPE_PREPARE_SWF.out.bed
        }
    }

    //
    // Optional LeafCutter sQTL phenotype (default OFF)
    //
    if (run_leaf) {
        def sqtl_meta = [id: params.sqtl_id ?: 'sqtl_leafcutter']
        if (!params.sqtl_counts) {
            log.warn "run_sqtl_leafcutter/run_sqtl=true but missing --sqtl_counts."
        } else {
            def ch_sqtl_counts = channel.of([sqtl_meta, file(params.sqtl_counts, checkIfExists: true)])
            def ch_sqtl_genes = params.sqtl_genes
                ? channel.of([sqtl_meta, file(params.sqtl_genes, checkIfExists: true)])
                : channel.empty()
            SQTL_LEAFCUTTER(ch_sqtl_counts, ch_sqtl_genes)
            ch_versions = ch_versions.mix(SQTL_LEAFCUTTER.out.versions)
            ch_prepared_pheno = ch_prepared_pheno.mix(SQTL_LEAFCUTTER.out.bed)
        }
    }

    //
    // Optional PEER hidden covariates (default OFF)
    // Prefer --peer_phenotype; else reuse phenotype_prepare BED.
    //
    if (run_peer_flag) {
        def peer_meta = [id: params.peer_id ?: 'peer']
        def ch_peer_pheno = params.peer_phenotype
            ? channel.of([peer_meta, file(params.peer_phenotype, checkIfExists: true)])
            : ch_prepared_pheno.map { _meta, bed -> [peer_meta, bed] }
        def ch_peer_known = params.peer_covariates
            ? channel.of([peer_meta, file(params.peer_covariates, checkIfExists: true)])
            : channel.empty()
        if (!params.peer_phenotype && !run_pheno_prep && !run_leaf) {
            log.warn "run_peer/modality=true but missing --peer_phenotype (and no prepared phenotype)."
        }
        COVARIATE_PEER(ch_peer_pheno, ch_peer_known)
        ch_versions = ch_versions.mix(COVARIATE_PEER.out.versions)
        ch_peer_cov_omiga = COVARIATE_PEER.out.cov_omiga
        ch_peer_cov_tensor = COVARIATE_PEER.out.cov_tensorqtl
    }

    //
    // Optional OmiGA cis-molQTL (default OFF)
    // Prefer --genotype_input molqtl_phenotype[/covariates] (+ optional QC bed).
    // Legacy: --omiga_cis_bed/bim/fam + --omiga_cis_phenotype, or --omiga_cis_use_qc_bed.
    //
    if (run_omiga) {
        def omiga_meta = [id: params.omiga_cis_id ?: 'omiga_cis']
        def ch_omiga_plink = channel.empty()
        def ch_omiga_pheno = channel.empty()
        def ch_omiga_covar = channel.empty()
        def ch_omiga_vcf = channel.empty()

        if (params.genotype_input) {
            if (params.omiga_cis_use_qc_bed || run_snp_indel || run_modality) {
                if (run_qc) {
                    ch_omiga_plink = ch_shared_qc_bed
                } else {
                    log.warn "omiga_cis_use_qc_bed/run_snp_indel/modality=true but QC did not run; falling back to samplesheet genotype bed."
                    ch_omiga_plink = ch_shared_geno_plink
                }
            } else {
                ch_omiga_plink = ch_shared_geno_plink
            }

            def ch_param_molqtl_pheno = params.omiga_cis_phenotype
                ? ch_omiga_plink.map { meta, _bed, _bim, _fam ->
                    [meta, file(params.omiga_cis_phenotype, checkIfExists: true)]
                }
                : channel.empty()
            ch_omiga_pheno = ch_molqtl_pheno.ifEmpty(ch_param_molqtl_pheno)

            def ch_param_molqtl_covar = params.omiga_cis_covariates
                ? ch_omiga_plink.map { meta, _bed, _bim, _fam ->
                    [meta, file(params.omiga_cis_covariates, checkIfExists: true)]
                }
                : channel.empty()
            ch_omiga_covar = ch_molqtl_covar.ifEmpty(ch_param_molqtl_covar)

            ch_omiga_vcf = ch_shared_geno_vcf
            if (params.omiga_cis_vcf) {
                ch_omiga_vcf = ch_omiga_vcf.ifEmpty(
                    ch_omiga_plink.map { meta, _bed, _bim, _fam ->
                        [meta, file(params.omiga_cis_vcf, checkIfExists: true)]
                    }
                )
            }
        } else {
            if (params.omiga_cis_use_qc_bed) {
                ch_omiga_plink = ch_shared_qc_bed.map { _meta, bed, bim, fam ->
                    [omiga_meta, bed, bim, fam]
                }
            } else if (params.omiga_cis_bed && params.omiga_cis_bim && params.omiga_cis_fam) {
                ch_omiga_plink = channel.of([
                    omiga_meta,
                    file(params.omiga_cis_bed, checkIfExists: true),
                    file(params.omiga_cis_bim, checkIfExists: true),
                    file(params.omiga_cis_fam, checkIfExists: true)
                ])
            }

            ch_omiga_pheno = params.omiga_cis_phenotype
                ? channel.of([omiga_meta, file(params.omiga_cis_phenotype, checkIfExists: true)])
                : channel.empty()

            ch_omiga_covar = params.omiga_cis_covariates
                ? channel.of([omiga_meta, file(params.omiga_cis_covariates, checkIfExists: true)])
                : channel.empty()

            ch_omiga_vcf = params.omiga_cis_vcf
                ? channel.of([omiga_meta, file(params.omiga_cis_vcf, checkIfExists: true)])
                : channel.empty()

            if (!params.omiga_cis_use_qc_bed &&
                (!params.omiga_cis_bed || !params.omiga_cis_bim || !params.omiga_cis_fam)) {
                log.warn "run_omiga_cis=true but missing --omiga_cis_bed/bim/fam (or --omiga_cis_use_qc_bed); channels empty."
            }
            if (!params.omiga_cis_phenotype && !run_pheno_prep && !run_leaf) {
                log.warn "run_omiga_cis=true but missing --omiga_cis_phenotype (and no prepared phenotype)."
            }
        }

        ch_omiga_pheno = ch_omiga_pheno.ifEmpty(
            ch_prepared_pheno.map { _meta, bed -> [omiga_meta, bed] }
        )
        ch_omiga_covar = ch_omiga_covar.ifEmpty(
            ch_peer_cov_omiga.map { _meta, cov -> [omiga_meta, cov] }
        )

        MOLQTL_MAP_OMIGA(
            ch_omiga_plink,
            ch_omiga_pheno,
            ch_omiga_covar,
            ch_omiga_vcf
        )
        ch_versions = ch_versions.mix(MOLQTL_MAP_OMIGA.out.versions)
        ch_qtl_cis_for_finemap = ch_qtl_cis_for_finemap.mix(
            MOLQTL_MAP_OMIGA.out.cis_qtl.map { meta, cis -> [meta + [engine: 'omiga'], cis] }
        )
    }

    //
    // Optional tensorQTL cis-molQTL (default OFF)
    // Prefer --genotype_input molqtl_phenotype[/covariates] (+ optional QC bed).
    // Legacy: --tensorqtl_cis_bed/bim/fam + --tensorqtl_cis_phenotype, or --tensorqtl_use_qc_bed.
    //
    if (run_tensor) {
        def tensorqtl_meta = [id: params.tensorqtl_cis_id ?: 'tensorqtl_cis']
        def ch_tensorqtl_plink = channel.empty()
        def ch_tensorqtl_pheno = channel.empty()
        def ch_tensorqtl_covar = channel.empty()

        if (params.genotype_input) {
            if (params.tensorqtl_use_qc_bed || run_modality) {
                if (run_qc) {
                    ch_tensorqtl_plink = ch_shared_qc_bed
                } else {
                    log.warn "tensorqtl_use_qc_bed/modality=true but QC did not run; falling back to samplesheet genotype bed."
                    ch_tensorqtl_plink = ch_shared_geno_plink
                }
            } else {
                ch_tensorqtl_plink = ch_shared_geno_plink
            }

            def ch_param_tqtl_pheno = params.tensorqtl_cis_phenotype
                ? ch_tensorqtl_plink.map { meta, _bed, _bim, _fam ->
                    [meta, file(params.tensorqtl_cis_phenotype, checkIfExists: true)]
                }
                : channel.empty()
            ch_tensorqtl_pheno = ch_molqtl_pheno.ifEmpty(ch_param_tqtl_pheno)

            def ch_param_tqtl_covar = params.tensorqtl_cis_covariates
                ? ch_tensorqtl_plink.map { meta, _bed, _bim, _fam ->
                    [meta, file(params.tensorqtl_cis_covariates, checkIfExists: true)]
                }
                : channel.empty()
            ch_tensorqtl_covar = ch_molqtl_covar.ifEmpty(ch_param_tqtl_covar)
        } else {
            if (params.tensorqtl_use_qc_bed) {
                ch_tensorqtl_plink = ch_shared_qc_bed.map { _meta, bed, bim, fam ->
                    [tensorqtl_meta, bed, bim, fam]
                }
            } else if (params.tensorqtl_cis_bed && params.tensorqtl_cis_bim && params.tensorqtl_cis_fam) {
                ch_tensorqtl_plink = channel.of([
                    tensorqtl_meta,
                    file(params.tensorqtl_cis_bed, checkIfExists: true),
                    file(params.tensorqtl_cis_bim, checkIfExists: true),
                    file(params.tensorqtl_cis_fam, checkIfExists: true)
                ])
            }

            ch_tensorqtl_pheno = params.tensorqtl_cis_phenotype
                ? channel.of([tensorqtl_meta, file(params.tensorqtl_cis_phenotype, checkIfExists: true)])
                : channel.empty()

            ch_tensorqtl_covar = params.tensorqtl_cis_covariates
                ? channel.of([tensorqtl_meta, file(params.tensorqtl_cis_covariates, checkIfExists: true)])
                : channel.empty()

            if (!params.tensorqtl_use_qc_bed &&
                (!params.tensorqtl_cis_bed || !params.tensorqtl_cis_bim || !params.tensorqtl_cis_fam)) {
                log.warn "run_tensorqtl_cis=true but missing --tensorqtl_cis_bed/bim/fam (or --tensorqtl_use_qc_bed); channels empty."
            }
            if (!params.tensorqtl_cis_phenotype && !params.run_phenotype_prepare) {
                log.warn "run_tensorqtl_cis=true but missing --tensorqtl_cis_phenotype (and phenotype_prepare did not run)."
            }
        }

        ch_tensorqtl_pheno = ch_tensorqtl_pheno.ifEmpty(
            ch_prepared_pheno.map { _meta, bed -> [tensorqtl_meta, bed] }
        )
        ch_tensorqtl_covar = ch_tensorqtl_covar.ifEmpty(
            ch_peer_cov_tensor.map { _meta, cov -> [tensorqtl_meta, cov] }
        )

        MOLQTL_MAP_TENSORQTL(
            ch_tensorqtl_plink,
            ch_tensorqtl_pheno,
            ch_tensorqtl_covar
        )
        ch_versions = ch_versions.mix(MOLQTL_MAP_TENSORQTL.out.versions)
        ch_qtl_cis_for_finemap = ch_qtl_cis_for_finemap.mix(
            MOLQTL_MAP_TENSORQTL.out.cis_qtl.map { meta, cis -> [meta + [engine: 'tensorqtl'], cis] }
        )
    }

    //
    // Optional QTLtools cis-molQTL (default OFF)
    // Prefer --genotype_input molqtl_phenotype[/covariates] (+ optional QC bed).
    //
    if (run_qtltools) {
        def qtltools_meta = [id: params.qtltools_cis_id ?: 'qtltools_cis']
        def ch_qtltools_plink = channel.empty()
        def ch_qtltools_pheno = channel.empty()
        def ch_qtltools_covar = channel.empty()

        if (params.genotype_input) {
            if (params.qtltools_use_qc_bed || run_modality) {
                if (run_qc) {
                    ch_qtltools_plink = ch_shared_qc_bed
                } else {
                    log.warn "qtltools_use_qc_bed/modality=true but QC did not run; falling back to samplesheet genotype bed."
                    ch_qtltools_plink = ch_shared_geno_plink
                }
            } else {
                ch_qtltools_plink = ch_shared_geno_plink
            }

            def ch_param_qt_pheno = params.qtltools_cis_phenotype
                ? ch_qtltools_plink.map { meta, _bed, _bim, _fam ->
                    [meta, file(params.qtltools_cis_phenotype, checkIfExists: true)]
                }
                : channel.empty()
            ch_qtltools_pheno = ch_molqtl_pheno.ifEmpty(ch_param_qt_pheno)

            def ch_param_qt_covar = params.qtltools_cis_covariates
                ? ch_qtltools_plink.map { meta, _bed, _bim, _fam ->
                    [meta, file(params.qtltools_cis_covariates, checkIfExists: true)]
                }
                : channel.empty()
            ch_qtltools_covar = ch_molqtl_covar.ifEmpty(ch_param_qt_covar)
        } else {
            if (params.qtltools_use_qc_bed) {
                ch_qtltools_plink = ch_shared_qc_bed.map { _meta, bed, bim, fam ->
                    [qtltools_meta, bed, bim, fam]
                }
            } else if (params.qtltools_cis_bed && params.qtltools_cis_bim && params.qtltools_cis_fam) {
                ch_qtltools_plink = channel.of([
                    qtltools_meta,
                    file(params.qtltools_cis_bed, checkIfExists: true),
                    file(params.qtltools_cis_bim, checkIfExists: true),
                    file(params.qtltools_cis_fam, checkIfExists: true)
                ])
            }

            ch_qtltools_pheno = params.qtltools_cis_phenotype
                ? channel.of([qtltools_meta, file(params.qtltools_cis_phenotype, checkIfExists: true)])
                : channel.empty()

            ch_qtltools_covar = params.qtltools_cis_covariates
                ? channel.of([qtltools_meta, file(params.qtltools_cis_covariates, checkIfExists: true)])
                : channel.empty()

            if (!params.qtltools_use_qc_bed &&
                (!params.qtltools_cis_bed || !params.qtltools_cis_bim || !params.qtltools_cis_fam)) {
                log.warn "run_qtltools_cis=true but missing --qtltools_cis_bed/bim/fam (or --qtltools_use_qc_bed); channels empty."
            }
            if (!params.qtltools_cis_phenotype && !params.run_phenotype_prepare) {
                log.warn "run_qtltools_cis=true but missing --qtltools_cis_phenotype (and phenotype_prepare did not run)."
            }
        }

        ch_qtltools_pheno = ch_qtltools_pheno.ifEmpty(
            ch_prepared_pheno.map { _meta, bed -> [qtltools_meta, bed] }
        )
        ch_qtltools_covar = ch_qtltools_covar.ifEmpty(
            ch_peer_cov_omiga.map { _meta, cov -> [qtltools_meta, cov] }
        )

        MOLQTL_MAP_QTLTOOLS(
            ch_qtltools_plink,
            ch_qtltools_pheno,
            ch_qtltools_covar
        )
        ch_versions = ch_versions.mix(MOLQTL_MAP_QTLTOOLS.out.versions)
        ch_qtl_cis_for_finemap = ch_qtl_cis_for_finemap.mix(
            MOLQTL_MAP_QTLTOOLS.out.cis_qtl.map { meta, cis -> [meta + [engine: 'qtltools'], cis] }
        )
    }

    //
    // Optional SuSiE fine-mapping (default OFF)
    // Prefer --finemap_susie_sumstats; else reuse cis QTL outputs from OmiGA/tensorQTL when those ran.
    //
    if (params.run_finemap_susie) {
        def finemap_meta = [id: params.finemap_susie_id ?: 'finemap_susie']
        def ch_finemap_sumstats = channel.empty()
        def ch_finemap_ld = channel.empty()

        if (params.finemap_susie_sumstats) {
            ch_finemap_sumstats = channel.of([
                finemap_meta,
                file(params.finemap_susie_sumstats, checkIfExists: true)
            ])
        } else {
            ch_finemap_sumstats = ch_qtl_cis_for_finemap
            if (!run_omiga && !run_tensor && !run_qtltools) {
                log.warn "run_finemap_susie=true but missing --finemap_susie_sumstats (and no cis QTL engine outputs)."
            }
        }

        if (params.finemap_susie_ld) {
            ch_finemap_ld = ch_finemap_sumstats.map { meta, _sumstats ->
                [meta, file(params.finemap_susie_ld, checkIfExists: true)]
            }
        }

        QTL_FINEMAP_SUSIE(
            ch_finemap_sumstats,
            ch_finemap_ld
        )
        ch_versions = ch_versions.mix(QTL_FINEMAP_SUSIE.out.versions)
    }

    //
    // Optional QTL–GWAS coloc.abf (default OFF). hyprcoloc is not wired yet.
    // Prefer --coloc_qtl_sumstats / --coloc_gwas_sumstats; else reuse cis QTL + GWAS standardized.
    //
    if (params.run_coloc) {
        def coloc_meta = [id: params.coloc_id ?: 'coloc']
        def ch_coloc_qtl = params.coloc_qtl_sumstats
            ? channel.of([coloc_meta, file(params.coloc_qtl_sumstats, checkIfExists: true)])
            : ch_qtl_cis_for_finemap
        def ch_coloc_gwas = params.coloc_gwas_sumstats
            ? channel.of([coloc_meta, file(params.coloc_gwas_sumstats, checkIfExists: true)])
            : ch_gwas_std

        if (!params.coloc_qtl_sumstats && !run_omiga && !run_tensor && !run_qtltools) {
            log.warn "run_coloc=true but missing --coloc_qtl_sumstats (and no cis QTL engine outputs)."
        }
        if (!params.coloc_gwas_sumstats && !run_gwas) {
            log.warn "run_coloc=true but missing --coloc_gwas_sumstats (and no GWAS standardized tables)."
        }

        QTL_COLOC(ch_coloc_qtl, ch_coloc_gwas)
        ch_versions = ch_versions.mix(QTL_COLOC.out.versions)
    }

    //
    // Optional HyPrColoc-style multi-trait clustering (default OFF)
    //
    if (params.run_hyprcoloc) {
        def hypr_meta = [id: params.hyprcoloc_id ?: 'hyprcoloc']
        if (!params.hyprcoloc_sumstats) {
            log.warn "run_hyprcoloc=true but missing --hyprcoloc_sumstats."
        } else {
            QTL_HYPRCOLOC(
                channel.of([hypr_meta, file(params.hyprcoloc_sumstats, checkIfExists: true)])
            )
            ch_versions = ch_versions.mix(QTL_HYPRCOLOC.out.versions)
        }
    }

    //
    // Optional SMR / HEIDI (default OFF)
    //
    if (params.run_smr) {
        def smr_meta = [id: params.smr_id ?: 'smr']
        def ch_smr_qtl = params.smr_qtl_sumstats
            ? channel.of([smr_meta, file(params.smr_qtl_sumstats, checkIfExists: true)])
            : ch_qtl_cis_for_finemap
        def ch_smr_gwas = params.smr_gwas_sumstats
            ? channel.of([smr_meta, file(params.smr_gwas_sumstats, checkIfExists: true)])
            : ch_gwas_std
        if (!params.smr_qtl_sumstats && !run_omiga && !run_tensor && !run_qtltools) {
            log.warn "run_smr=true but missing --smr_qtl_sumstats (and no cis QTL engine outputs)."
        }
        if (!params.smr_gwas_sumstats && !run_gwas) {
            log.warn "run_smr=true but missing --smr_gwas_sumstats (and no GWAS standardized tables)."
        }
        QTL_SMR(ch_smr_qtl, ch_smr_gwas)
        ch_versions = ch_versions.mix(QTL_SMR.out.versions)
    }

    //
    // Optional mashr multi-condition shrinkage (default OFF)
    //
    if (params.run_mashr) {
        def mash_meta = [id: params.mashr_id ?: 'mashr']
        if (!params.mashr_sumstats) {
            log.warn "run_mashr=true but missing --mashr_sumstats."
        } else {
            QTL_MASHR(
                channel.of([mash_meta, file(params.mashr_sumstats, checkIfExists: true)])
            )
            ch_versions = ch_versions.mix(QTL_MASHR.out.versions)
        }
    }

    //
    // Optional METAL inverse-variance meta (default OFF)
    //
    if (params.run_metal) {
        def metal_meta = [id: params.metal_id ?: 'metal']
        if (!params.metal_sumstats) {
            log.warn "run_metal=true but missing --metal_sumstats."
        } else {
            QTL_METAL(
                channel.of([metal_meta, file(params.metal_sumstats, checkIfExists: true)])
            )
            ch_versions = ch_versions.mix(QTL_METAL.out.versions)
        }
    }

    //
    // Optional TORUS enrichment prior (default OFF)
    //
    if (params.run_torus) {
        def torus_meta = [id: params.torus_id ?: 'torus']
        if (!params.torus_annot) {
            log.warn "run_torus=true but missing --torus_annot."
        } else {
            QTL_TORUS(
                channel.of([torus_meta, file(params.torus_annot, checkIfExists: true)])
            )
            ch_versions = ch_versions.mix(QTL_TORUS.out.versions)
        }
    }

    //
    // Optional TWAS (default OFF)
    //
    if (params.run_twas) {
        def twas_meta = [id: params.twas_id ?: 'twas']
        def ch_twas_w = params.twas_weights
            ? channel.of([twas_meta, file(params.twas_weights, checkIfExists: true)])
            : channel.empty()
        def ch_twas_g = params.twas_gwas
            ? channel.of([twas_meta, file(params.twas_gwas, checkIfExists: true)])
            : ch_gwas_std
        def ch_twas_ld = params.twas_ld
            ? channel.of([twas_meta, file(params.twas_ld, checkIfExists: true)])
            : channel.empty()
        if (!params.twas_weights) {
            log.warn "run_twas=true but missing --twas_weights."
        }
        if (!params.twas_gwas && !run_gwas) {
            log.warn "run_twas=true but missing --twas_gwas (and no GWAS standardized tables)."
        }
        if (params.twas_weights) {
            QTL_TWAS(ch_twas_w, ch_twas_g, ch_twas_ld)
            ch_versions = ch_versions.mix(QTL_TWAS.out.versions)
        }
    }

    //
    // Optional LDSC heritability (default OFF)
    //
    if (params.run_ldsc) {
        def ldsc_meta = [id: params.ldsc_id ?: 'ldsc']
        def ch_ldsc_ss = params.ldsc_sumstats
            ? channel.of([ldsc_meta, file(params.ldsc_sumstats, checkIfExists: true)])
            : ch_gwas_std
        def ch_ldsc_annot = params.ldsc_annot
            ? channel.of([ldsc_meta, file(params.ldsc_annot, checkIfExists: true)])
            : channel.empty()
        def ch_ldsc_l2 = params.ldsc_ldscores
            ? channel.of([ldsc_meta, file(params.ldsc_ldscores, checkIfExists: true)])
            : channel.empty()
        if (!params.ldsc_sumstats && !run_gwas) {
            log.warn "run_ldsc=true but missing --ldsc_sumstats (and no GWAS standardized tables)."
        }
        if (params.ldsc_sumstats || run_gwas) {
            QTL_LDSC(ch_ldsc_ss, ch_ldsc_annot, ch_ldsc_l2)
            ch_versions = ch_versions.mix(QTL_LDSC.out.versions)
        }
    }

    //
    // Optional FINEMAP/CAVIAR/DAP-G PIPs (default OFF)
    //
    if (params.run_finemap_extra) {
        def extra_meta = [id: params.finemap_extra_id ?: 'finemap_extra']
        def ch_extra = params.finemap_extra_sumstats
            ? channel.of([extra_meta, file(params.finemap_extra_sumstats, checkIfExists: true)])
            : ch_qtl_cis_for_finemap
        if (!params.finemap_extra_sumstats && !run_omiga && !run_tensor && !run_qtltools) {
            log.warn "run_finemap_extra=true but missing --finemap_extra_sumstats (and no cis QTL engine outputs)."
        }
        QTL_FINEMAP_EXTRA(ch_extra)
        ch_versions = ch_versions.mix(QTL_FINEMAP_EXTRA.out.versions)
    }

    //
    // Optional annotation cache download (default OFF)
    //
    if (params.run_cache) {
        def cache_meta = [id: params.cache_id ?: 'cache']
        def tools = (params.cache_tools ?: 'snpeff,ensemblvep')
            .tokenize(',')
            .collect { tool -> tool.trim().toLowerCase() }
        def ch_vep_info = tools.contains('ensemblvep')
            ? channel.of([
                cache_meta,
                params.annotate_vep_genome ?: (params.genome ?: 'GRCh38'),
                params.annotate_vep_species ?: 'homo_sapiens',
                params.annotate_vep_cache_version ?: '110'
            ])
            : channel.empty()
        def ch_snpeff_info = tools.contains('snpeff')
            ? channel.of([cache_meta, params.annotate_snpeff_db ?: 'GRCh38.99'])
            : channel.empty()
        ANNOTATION_CACHE(ch_vep_info, ch_snpeff_info)
        ch_versions = ch_versions.mix(ANNOTATION_CACHE.out.versions)
        ch_downloaded_vep = ANNOTATION_CACHE.out.vep_cache
        ch_downloaded_snpeff = ANNOTATION_CACHE.out.snpeff_cache
    }

    //
    // Optional FASTA bgzip + faidx/dict (default OFF)
    //
    if (params.run_fasta_index) {
        def fa_meta = [id: params.fasta_index_id ?: 'fasta_index']
        if (!params.fasta_index_fasta) {
            log.warn "run_fasta_index=true but missing --fasta_index_fasta."
        } else {
            REFERENCE_FASTA(
                channel.of([fa_meta, file(params.fasta_index_fasta, checkIfExists: true)])
            )
            ch_versions = ch_versions.mix(REFERENCE_FASTA.out.versions)
            ch_fasta_indexed = REFERENCE_FASTA.out.fasta_fai_gzi_dict
        }
    }

    //
    // Optional Somalier relatedness (default OFF)
    //
    if (params.run_relate) {
        def rel_meta = [id: params.relate_id ?: 'relate']
        def have_rel_fa = params.relate_fasta || fasta_index_ok
        def have_rel_fai = params.relate_fasta_fai || fasta_index_ok
        if (!params.relate_vcf || !have_rel_fa || !have_rel_fai || !params.relate_sites) {
            log.warn "run_relate=true but missing --relate_vcf/--relate_sites and FASTA (--relate_fasta/--relate_fasta_fai or run_fasta_index)."
        } else {
            def ch_rel_fa = params.relate_fasta
                ? channel.of([[id: 'relate_ref'], file(params.relate_fasta, checkIfExists: true)])
                : ch_fasta_indexed.map { _meta, fa, _fai, _gzi, _sizes, _dict -> [[id: 'relate_ref'], fa] }
            def ch_rel_fai = params.relate_fasta_fai
                ? channel.of([[id: 'relate_ref'], file(params.relate_fasta_fai, checkIfExists: true)])
                : ch_fasta_indexed.map { _meta, _fa, fai, _gzi, _sizes, _dict -> [[id: 'relate_ref'], fai] }
            VARIANT_RELATE(
                channel.of([
                    rel_meta,
                    file(params.relate_vcf, checkIfExists: true),
                    params.relate_vcf_tbi
                        ? file(params.relate_vcf_tbi, checkIfExists: true)
                        : []
                ]),
                ch_rel_fa,
                ch_rel_fai,
                channel.of([[id: 'relate_sites'], file(params.relate_sites, checkIfExists: true)]),
                channel.of([
                    rel_meta,
                    params.relate_ped ? file(params.relate_ped, checkIfExists: true) : []
                ])
            )
            ch_versions = ch_versions.mix(VARIANT_RELATE.out.versions)
            ch_somalier_pairs = VARIANT_RELATE.out.pairs_tsv
        }
    }

    //
    // Optional: apply Somalier pairs as relatedness removals on QC bed
    //
    if (params.genotype_qc_use_somalier) {
        def ch_pairs_for_qc = params.genotype_qc_somalier_pairs
            ? channel.of([
                [id: params.relate_id ?: 'relate'],
                file(params.genotype_qc_somalier_pairs, checkIfExists: true)
            ])
            : ch_somalier_pairs
        if (!params.genotype_qc_somalier_pairs && !params.run_relate) {
            log.warn "genotype_qc_use_somalier=true but missing --genotype_qc_somalier_pairs (and run_relate did not run)."
        }
        SOMALIER_OUTLIERS(ch_pairs_for_qc)
        ch_versions = ch_versions.mix(SOMALIER_OUTLIERS.out.versions)
        def ch_som_remove = SOMALIER_OUTLIERS.out.outliers
            .filter { _meta, path -> path.size() > 0 }
        def ch_bed_to_filter = ch_shared_qc_bed.join(ch_som_remove)
        PLINK2_REMOVE_SOMALIER(
            ch_bed_to_filter.map { meta, bed, bim, fam, _out -> [meta, bed, bim, fam] },
            ch_bed_to_filter.map { _meta, _bed, _bim, _fam, out -> out }
        )
        ch_versions = ch_versions.mix(PLINK2_REMOVE_SOMALIER.out.versions)
        ch_shared_qc_bed = ch_shared_qc_bed
            .join(SOMALIER_OUTLIERS.out.outliers)
            .filter { _meta, _bed, _bim, _fam, path -> path.size() == 0 }
            .map { meta, bed, bim, fam, _path -> [meta, bed, bim, fam] }
            .mix(
                PLINK2_REMOVE_SOMALIER.out.remove_bed
                    .join(PLINK2_REMOVE_SOMALIER.out.remove_bim)
                    .join(PLINK2_REMOVE_SOMALIER.out.remove_fam)
            )
        ch_shared_geno_plink = ch_shared_qc_bed
    }

    //
    // Optional VCF annotation (default OFF)
    //
    if (run_annotate_flag) {
        def ann_meta = [id: params.annotate_id ?: 'annotate']
        def ann_vcf = params.annotate_vcf ?: vcf_prep_vcf
        def ann_tbi = params.annotate_vcf_tbi ?: vcf_prep_tbi
        if (!ann_vcf) {
            log.warn "run_annotate=true but missing --annotate_vcf."
        } else {
            def ch_ann_vcf = channel.of([
                ann_meta,
                file(ann_vcf, checkIfExists: true),
                ann_tbi ? file(ann_tbi, checkIfExists: true) : []
            ])
            def ch_ann_fa = params.annotate_fasta
                ? channel.of([[id: 'annotate_fasta'], file(params.annotate_fasta, checkIfExists: true)])
                : (fasta_index_ok
                    ? ch_fasta_indexed.map { _meta, fa, _fai, _gzi, _sizes, _dict -> [[id: 'annotate_fasta'], fa] }
                    : channel.of([[id: 'annotate_fasta'], []]))
            def ch_vep_cache = params.annotate_vep_cache
                ? channel.of([[id: 'vep_cache'], file(params.annotate_vep_cache, checkIfExists: true)])
                : (params.run_cache ? ch_downloaded_vep : channel.of([[id: 'vep_cache'], []]))
            def ch_snpeff_cache = params.annotate_snpeff_cache
                ? channel.of([[id: 'snpeff_cache'], file(params.annotate_snpeff_cache, checkIfExists: true)])
                : (params.run_cache ? ch_downloaded_snpeff : channel.of([[id: 'snpeff_cache'], []]))
            VARIANT_ANNOTATE(ch_ann_vcf, ch_ann_fa, ch_vep_cache, ch_snpeff_cache)
            ch_versions = ch_versions.mix(VARIANT_ANNOTATE.out.versions)
        }
    }

    //
    // Optional VCF phasing (default OFF)
    //
    if (run_phase_flag) {
        def phase_meta = [id: params.phase_id ?: 'phase']
        def phase_vcf = params.phase_vcf ?: vcf_prep_vcf
        def phase_tbi = params.phase_vcf_tbi ?: vcf_prep_tbi
        if (!phase_vcf) {
            log.warn "run_phase=true but missing --phase_vcf."
        } else {
            def ch_phase_vcf = channel.of([
                phase_meta,
                file(phase_vcf, checkIfExists: true),
                phase_tbi ? file(phase_tbi, checkIfExists: true) : []
            ])
            def ch_phase_ref = (params.phase_ref_vcf)
                ? channel.of([
                    phase_meta,
                    file(params.phase_ref_vcf, checkIfExists: true),
                    params.phase_ref_vcf_tbi
                        ? file(params.phase_ref_vcf_tbi, checkIfExists: true)
                        : []
                ])
                : channel.empty()
            def ch_phase_map = params.phase_map
                ? channel.of([phase_meta, file(params.phase_map, checkIfExists: true)])
                : channel.empty()
            VARIANT_PHASE(ch_phase_vcf, ch_phase_ref, ch_phase_map)
            ch_versions = ch_versions.mix(VARIANT_PHASE.out.versions)
        }
    }

    //
    // Optional genotype imputation (default OFF)
    //
    if (run_impute_flag) {
        def imp_meta = [id: params.impute_id ?: 'impute']
        def imp_vcf = params.impute_vcf ?: vcf_prep_vcf
        def imp_tbi = params.impute_vcf_tbi ?: vcf_prep_tbi
        if (!imp_vcf || !params.impute_panel) {
            log.warn "run_impute=true but missing --impute_vcf and/or --impute_panel."
        } else {
            def ch_imp_vcf = channel.of([
                imp_meta,
                file(imp_vcf, checkIfExists: true),
                imp_tbi ? file(imp_tbi, checkIfExists: true) : []
            ])
            def ch_imp_panel = channel.of([
                imp_meta,
                file(params.impute_panel, checkIfExists: true),
                params.impute_panel_tbi
                    ? file(params.impute_panel_tbi, checkIfExists: true)
                    : []
            ])
            def ch_imp_map = params.impute_map
                ? channel.of([imp_meta, file(params.impute_map, checkIfExists: true)])
                : channel.empty()
            VARIANT_IMPUTE(ch_imp_vcf, ch_imp_panel, ch_imp_map)
            ch_versions = ch_versions.mix(VARIANT_IMPUTE.out.versions)
        }
    }

    //
    // Optional chained annotate → phase → impute (default OFF)
    //
    if (run_vcf_prep) {
        def prep_meta = [id: params.vcf_prep_id ?: 'vcf_prep']
        if (!vcf_prep_vcf) {
            log.warn "run_vcf_prep=true but missing --vcf_prep_vcf (or --annotate_vcf/--phase_vcf/--impute_vcf)."
        } else if (!params.vcf_prep_skip_impute && !params.impute_panel) {
            log.warn "run_vcf_prep=true but missing --impute_panel (or set --vcf_prep_skip_impute)."
        } else {
            def ch_prep_vcf = channel.of([
                prep_meta,
                file(vcf_prep_vcf, checkIfExists: true),
                vcf_prep_tbi ? file(vcf_prep_tbi, checkIfExists: true) : []
            ])
            def ch_prep_fa = params.annotate_fasta
                ? channel.of([prep_meta, file(params.annotate_fasta, checkIfExists: true)])
                : (fasta_index_ok
                    ? ch_fasta_indexed.map { _meta, fa, _fai, _gzi, _sizes, _dict -> [prep_meta, fa] }
                    : channel.of([prep_meta, []]))
            def ch_prep_vep = params.annotate_vep_cache
                ? channel.of([prep_meta, file(params.annotate_vep_cache, checkIfExists: true)])
                : (params.run_cache ? ch_downloaded_vep.map { _meta, cache -> [prep_meta, cache] } : channel.of([prep_meta, []]))
            def ch_prep_snpeff = params.annotate_snpeff_cache
                ? channel.of([prep_meta, file(params.annotate_snpeff_cache, checkIfExists: true)])
                : (params.run_cache ? ch_downloaded_snpeff.map { _meta, cache -> [prep_meta, cache] } : channel.of([prep_meta, []]))
            def ch_prep_pref = params.phase_ref_vcf
                ? channel.of([
                    prep_meta,
                    file(params.phase_ref_vcf, checkIfExists: true),
                    params.phase_ref_vcf_tbi ? file(params.phase_ref_vcf_tbi, checkIfExists: true) : []
                ])
                : channel.empty()
            def ch_prep_pmap = params.phase_map
                ? channel.of([prep_meta, file(params.phase_map, checkIfExists: true)])
                : channel.empty()
            def ch_prep_panel = params.impute_panel
                ? channel.of([
                    prep_meta,
                    file(params.impute_panel, checkIfExists: true),
                    params.impute_panel_tbi ? file(params.impute_panel_tbi, checkIfExists: true) : []
                ])
                : channel.empty()
            def ch_prep_imap = params.impute_map
                ? channel.of([prep_meta, file(params.impute_map, checkIfExists: true)])
                : channel.empty()
            VARIANT_VCF_PREP(
                ch_prep_vcf,
                ch_prep_fa,
                ch_prep_vep,
                ch_prep_snpeff,
                ch_prep_pref,
                ch_prep_pmap,
                ch_prep_panel,
                ch_prep_imap
            )
            ch_versions = ch_versions.mix(VARIANT_VCF_PREP.out.versions)
        }
    }

    //
    // Optional GLIMPSE2 BAM/GL imputation (default OFF)
    //
    if (params.run_impute_bam) {
        def bam_meta = [id: params.impute_bam_id ?: 'impute_bam']
        if (!params.impute_bam_input || !params.impute_panel) {
            log.warn "run_impute_bam=true but missing --impute_bam_input and/or --impute_panel."
        } else {
            def ch_bam_fa = params.impute_bam_fasta
                ? channel.of([
                    bam_meta,
                    file(params.impute_bam_fasta, checkIfExists: true),
                    params.impute_bam_fasta_fai
                        ? file(params.impute_bam_fasta_fai, checkIfExists: true)
                        : []
                ])
                : (fasta_index_ok
                    ? ch_fasta_indexed.map { _meta, fa, fai, _gzi, _sizes, _dict -> [bam_meta, fa, fai] }
                    : channel.of([bam_meta, [], []]))
            VARIANT_IMPUTE_BAM(
                channel.of([
                    bam_meta,
                    file(params.impute_bam_input, checkIfExists: true),
                    params.impute_bam_index
                        ? file(params.impute_bam_index, checkIfExists: true)
                        : []
                ]),
                channel.of([
                    bam_meta,
                    file(params.impute_panel, checkIfExists: true),
                    params.impute_panel_tbi
                        ? file(params.impute_panel_tbi, checkIfExists: true)
                        : []
                ]),
                ch_bam_fa,
                params.impute_map
                    ? channel.of([bam_meta, file(params.impute_map, checkIfExists: true)])
                    : channel.of([bam_meta, []])
            )
            ch_versions = ch_versions.mix(VARIANT_IMPUTE_BAM.out.versions)
        }
    }

    //
    // Optional SV calling (default OFF)
    //
    if (params.run_sv) {
        def sv_meta = [id: params.sv_id ?: 'sv']
        if (!params.sv_bam || !params.sv_bam_index || !(params.sv_fasta || fasta_index_ok) || !(params.sv_fasta_fai || fasta_index_ok)) {
            log.warn "run_sv=true but missing --sv_bam/--sv_bam_index and FASTA (--sv_fasta/--sv_fasta_fai or run_fasta_index)."
        } else {
            def ch_sv_fa = params.sv_fasta
                ? channel.of([[id: 'sv_ref'], file(params.sv_fasta, checkIfExists: true)])
                : ch_fasta_indexed.map { _meta, fa, _fai, _gzi, _sizes, _dict -> [[id: 'sv_ref'], fa] }
            def ch_sv_fai = params.sv_fasta_fai
                ? channel.of([[id: 'sv_ref'], file(params.sv_fasta_fai, checkIfExists: true)])
                : ch_fasta_indexed.map { _meta, _fa, fai, _gzi, _sizes, _dict -> [[id: 'sv_ref'], fai] }
            VARIANT_SV(
                channel.of([
                    sv_meta,
                    file(params.sv_bam, checkIfExists: true),
                    file(params.sv_bam_index, checkIfExists: true)
                ]),
                ch_sv_fa,
                ch_sv_fai
            )
            ch_versions = ch_versions.mix(VARIANT_SV.out.versions)
        }
    }

    //
    // Optional STR genotyping (default OFF)
    //
    if (params.run_str) {
        def str_meta = [id: params.str_id ?: 'str']
        if (!params.str_bam || !params.str_bam_index || !(params.str_fasta || fasta_index_ok) || !(params.str_fasta_fai || fasta_index_ok)) {
            log.warn "run_str=true but missing --str_bam/--str_bam_index and FASTA (--str_fasta/--str_fasta_fai or run_fasta_index)."
        } else {
            def ch_str_bam = channel.of([
                str_meta,
                file(params.str_bam, checkIfExists: true),
                file(params.str_bam_index, checkIfExists: true)
            ])
            def ch_str_fa = params.str_fasta
                ? channel.of([[id: 'str_ref'], file(params.str_fasta, checkIfExists: true)])
                : ch_fasta_indexed.map { _meta, fa, _fai, _gzi, _sizes, _dict -> [[id: 'str_ref'], fa] }
            def ch_str_fai = params.str_fasta_fai
                ? channel.of([[id: 'str_ref'], file(params.str_fasta_fai, checkIfExists: true)])
                : ch_fasta_indexed.map { _meta, _fa, fai, _gzi, _sizes, _dict -> [[id: 'str_ref'], fai] }
            def ch_str_cat = params.str_catalog
                ? channel.of([str_meta, file(params.str_catalog, checkIfExists: true)])
                : channel.empty()
            def ch_str_reg = params.str_regions
                ? channel.of([str_meta, file(params.str_regions, checkIfExists: true)])
                : channel.empty()
            def ch_str_rep = params.str_repeats
                ? channel.of([str_meta, file(params.str_repeats, checkIfExists: true)])
                : ch_str_reg
            def ch_str_hip = params.str_hipstr_bed
                ? channel.of([str_meta, file(params.str_hipstr_bed, checkIfExists: true)])
                : ch_str_reg
            VARIANT_STR(
                ch_str_bam,
                ch_str_fa,
                ch_str_fai,
                ch_str_cat,
                ch_str_reg,
                ch_str_rep,
                ch_str_hip
            )
            ch_versions = ch_versions.mix(VARIANT_STR.out.versions)
        }
    }

    //
    // Optional cis-QTL postprocess: unified table + BH q-values (default OFF)
    //
    if (params.run_qtl_postprocess) {
        def ch_post_in = params.qtl_postprocess_input
            ? channel.of([
                [id: params.qtl_postprocess_id ?: 'qtl_postprocess'],
                file(params.qtl_postprocess_input, checkIfExists: true)
            ])
            : ch_qtl_cis_for_finemap

        if (!params.qtl_postprocess_input && !run_omiga && !run_tensor && !run_qtltools) {
            log.warn "run_qtl_postprocess=true but missing --qtl_postprocess_input (and no cis QTL engine outputs)."
        }

        QTL_POSTPROCESS_CIS(ch_post_in)
        ch_versions = ch_versions.mix(QTL_POSTPROCESS_CIS.out.versions)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'variant2qtl_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'variant2qtl'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
