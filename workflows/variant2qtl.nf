/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variant2qtl_pipeline'
include { GWAS_BENCHMARK_PARALLEL   } from '../subworkflows/local/gwas_benchmark_parallel/main'
include { GENOTYPE_QC               } from '../subworkflows/local/genotype_qc/main'
include { GENOTYPE_INGEST_HARMONIZE } from '../subworkflows/local/genotype_ingest_harmonize/main'
include { MOLQTL_MAP_OMIGA          } from '../subworkflows/local/molqtl_map_omiga/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT2QTL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { entry -> entry[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Optional genotype ingest → QC → GWAS benchmark (all default OFF).
    // Prefer ingest VCF→bed over raw --gwas_benchmark_bed when ingest is ON.
    // When run_genotype_qc=true, filtered bed feeds the benchmark (VCF cleared so
    // fan-out regenerates from QC bed). Params-bed path still works when ingest OFF.
    //
    if (params.run_genotype_ingest || params.run_genotype_qc || params.run_gwas_benchmark) {
        def gwas_meta = [id: params.gwas_benchmark_id ?: 'gwas_benchmark']

        def ch_gwas_plink_raw = Channel.empty()
        def ch_gwas_vcf = Channel.empty()

        if (params.run_genotype_ingest) {
            if (!params.genotype_ingest_vcf) {
                log.warn "run_genotype_ingest=true but missing --genotype_ingest_vcf; channels empty."
            } else {
                def ch_ingest_vcf = Channel.of([
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
                ? Channel.of([
                    gwas_meta,
                    file(params.gwas_benchmark_bed, checkIfExists: true),
                    file(params.gwas_benchmark_bim, checkIfExists: true),
                    file(params.gwas_benchmark_fam, checkIfExists: true)
                ])
                : Channel.empty()

            ch_gwas_vcf = params.gwas_benchmark_vcf
                ? Channel.of([gwas_meta, file(params.gwas_benchmark_vcf, checkIfExists: true)])
                : Channel.empty()

            if ((params.run_genotype_qc || params.run_gwas_benchmark) &&
                (!params.gwas_benchmark_bed || !params.gwas_benchmark_bim || !params.gwas_benchmark_fam)) {
                log.warn "genotype_qc/gwas_benchmark enabled but missing --gwas_benchmark_bed/bim/fam; channels empty."
            }
        }

        def ch_gwas_pheno = params.gwas_benchmark_phenotype
            ? Channel.of([gwas_meta, file(params.gwas_benchmark_phenotype, checkIfExists: true)])
            : Channel.empty()

        def ch_gwas_covar = params.gwas_benchmark_covariates
            ? Channel.of([gwas_meta, file(params.gwas_benchmark_covariates, checkIfExists: true)])
            : Channel.empty()

        def ch_gwas_plink = ch_gwas_plink_raw
        if (params.run_genotype_qc) {
            GENOTYPE_QC(ch_gwas_plink_raw)
            ch_versions = ch_versions.mix(GENOTYPE_QC.out.versions)
            ch_gwas_plink = GENOTYPE_QC.out.bed
            // QC changes SNP set — force VCF rebuild from filtered bed
            ch_gwas_vcf = Channel.empty()
        }

        if (params.run_gwas_benchmark) {
            if (!params.gwas_benchmark_phenotype) {
                log.warn "run_gwas_benchmark=true but missing --gwas_benchmark_phenotype."
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
    // Inputs: --omiga_cis_bed/bim/fam + --omiga_cis_phenotype (+ optional covariates).
    //
    if (params.run_omiga_cis) {
        def omiga_meta = [id: params.omiga_cis_id ?: 'omiga_cis']

        def ch_omiga_plink = (params.omiga_cis_bed && params.omiga_cis_bim && params.omiga_cis_fam)
            ? Channel.of([
                omiga_meta,
                file(params.omiga_cis_bed, checkIfExists: true),
                file(params.omiga_cis_bim, checkIfExists: true),
                file(params.omiga_cis_fam, checkIfExists: true)
            ])
            : Channel.empty()

        def ch_omiga_pheno = params.omiga_cis_phenotype
            ? Channel.of([omiga_meta, file(params.omiga_cis_phenotype, checkIfExists: true)])
            : Channel.empty()

        def ch_omiga_covar = params.omiga_cis_covariates
            ? Channel.of([omiga_meta, file(params.omiga_cis_covariates, checkIfExists: true)])
            : Channel.empty()

        if (!params.omiga_cis_bed || !params.omiga_cis_bim || !params.omiga_cis_fam) {
            log.warn "run_omiga_cis=true but missing --omiga_cis_bed/bim/fam; channels empty."
        }
        if (!params.omiga_cis_phenotype) {
            log.warn "run_omiga_cis=true but missing --omiga_cis_phenotype."
        }

        MOLQTL_MAP_OMIGA(
            ch_omiga_plink,
            ch_omiga_pheno,
            ch_omiga_covar
        )
        ch_versions = ch_versions.mix(MOLQTL_MAP_OMIGA.out.versions)
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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'variant2qtl_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
