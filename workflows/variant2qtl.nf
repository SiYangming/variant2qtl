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
include { QTL_FINEMAP_SUSIE         } from '../subworkflows/local/qtl_finemap_susie/main'

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

    if (params.run_genotype_ingest || params.run_genotype_qc || params.run_gwas_benchmark || params.genotype_input) {
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

            if (params.run_genotype_ingest) {
                if (!params.genotype_ingest_vcf) {
                    log.warn "run_genotype_ingest=true but missing --genotype_ingest_vcf; channels empty."
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

                if ((params.run_genotype_qc || params.run_gwas_benchmark) &&
                    (!params.gwas_benchmark_bed || !params.gwas_benchmark_bim || !params.gwas_benchmark_fam)) {
                    log.warn "genotype_qc/gwas_benchmark enabled but missing --gwas_benchmark_bed/bim/fam; channels empty."
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
        if (params.run_genotype_qc) {
            GENOTYPE_QC(ch_gwas_plink_raw)
            ch_versions = ch_versions.mix(GENOTYPE_QC.out.versions)
            ch_gwas_plink = GENOTYPE_QC.out.bed
            ch_shared_qc_bed = GENOTYPE_QC.out.bed
            // QC changes SNP set — force VCF rebuild from filtered bed
            ch_gwas_vcf = channel.empty()
        }

        ch_shared_geno_plink = ch_gwas_plink
        ch_shared_geno_vcf = ch_gwas_vcf

        if (params.run_gwas_benchmark) {
            if (!params.genotype_input && !params.gwas_benchmark_phenotype) {
                log.warn "run_gwas_benchmark=true but missing --gwas_benchmark_phenotype (or --genotype_input)."
            }
            GWAS_BENCHMARK_PARALLEL(
                ch_gwas_plink,
                ch_gwas_vcf,
                ch_gwas_pheno,
                ch_gwas_covar
            )
            ch_versions = ch_versions.mix(GWAS_BENCHMARK_PARALLEL.out.versions)
        }
    }

    //
    // Optional OmiGA cis-molQTL (default OFF)
    // Prefer --genotype_input molqtl_phenotype[/covariates] (+ optional QC bed).
    // Legacy: --omiga_cis_bed/bim/fam + --omiga_cis_phenotype, or --omiga_cis_use_qc_bed.
    //
    if (params.run_omiga_cis) {
        def omiga_meta = [id: params.omiga_cis_id ?: 'omiga_cis']
        def ch_omiga_plink = channel.empty()
        def ch_omiga_pheno = channel.empty()
        def ch_omiga_covar = channel.empty()
        def ch_omiga_vcf = channel.empty()

        if (params.genotype_input) {
            if (params.omiga_cis_use_qc_bed) {
                if (params.run_genotype_qc) {
                    ch_omiga_plink = ch_shared_qc_bed
                } else {
                    log.warn "omiga_cis_use_qc_bed=true but run_genotype_qc=false; falling back to samplesheet genotype bed."
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
            if (!params.omiga_cis_phenotype) {
                log.warn "run_omiga_cis=true but missing --omiga_cis_phenotype."
            }
        }

        MOLQTL_MAP_OMIGA(
            ch_omiga_plink,
            ch_omiga_pheno,
            ch_omiga_covar,
            ch_omiga_vcf
        )
        ch_versions = ch_versions.mix(MOLQTL_MAP_OMIGA.out.versions)
        ch_qtl_cis_for_finemap = ch_qtl_cis_for_finemap.mix(MOLQTL_MAP_OMIGA.out.cis_qtl)
    }

    //
    // Optional tensorQTL cis-molQTL (default OFF)
    // Prefer --genotype_input molqtl_phenotype[/covariates] (+ optional QC bed).
    // Legacy: --tensorqtl_cis_bed/bim/fam + --tensorqtl_cis_phenotype, or --tensorqtl_use_qc_bed.
    //
    if (params.run_tensorqtl_cis) {
        def tensorqtl_meta = [id: params.tensorqtl_cis_id ?: 'tensorqtl_cis']
        def ch_tensorqtl_plink = channel.empty()
        def ch_tensorqtl_pheno = channel.empty()
        def ch_tensorqtl_covar = channel.empty()

        if (params.genotype_input) {
            if (params.tensorqtl_use_qc_bed) {
                if (params.run_genotype_qc) {
                    ch_tensorqtl_plink = ch_shared_qc_bed
                } else {
                    log.warn "tensorqtl_use_qc_bed=true but run_genotype_qc=false; falling back to samplesheet genotype bed."
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
            if (!params.tensorqtl_cis_phenotype) {
                log.warn "run_tensorqtl_cis=true but missing --tensorqtl_cis_phenotype."
            }
        }

        MOLQTL_MAP_TENSORQTL(
            ch_tensorqtl_plink,
            ch_tensorqtl_pheno,
            ch_tensorqtl_covar
        )
        ch_versions = ch_versions.mix(MOLQTL_MAP_TENSORQTL.out.versions)
        ch_qtl_cis_for_finemap = ch_qtl_cis_for_finemap.mix(MOLQTL_MAP_TENSORQTL.out.cis_qtl)
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
            if (!params.run_omiga_cis && !params.run_tensorqtl_cis) {
                log.warn "run_finemap_susie=true but missing --finemap_susie_sumstats (and no OmiGA/tensorQTL cis outputs)."
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
