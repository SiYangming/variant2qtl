//
// Parallel GWAS engine benchmark: GEMMA, EMMAX, TASSEL, rMVP, OmiGA.
// Requires includeConfig for genotype_to_gwas_formats PLINK_RECODE aliases:
//   includeConfig 'subworkflows/local/genotype_to_gwas_formats/nextflow.config'
//
// Select engines with params.gwas_benchmark_engines (comma-separated).
// Default: gemma,emmax,tassel,rmvp,omiga
// On arm64 hosts, omit emmax if the amd64 emmax binary hits Illegal instruction.
//

include { PHENOCOVAR_ADAPT                           } from '../../../modules/local/utils/phenocovar_adapt/main'
include { ASSOC_STANDARDIZE                          } from '../../../modules/local/utils/assoc_standardize/main'
include { GENOTYPE_TO_GWAS_FORMATS                   } from '../genotype_to_gwas_formats/main'
include { GEMMA_RELATEDNESS                          } from '../../../modules/local/gemma/relatedness/main'
include { GEMMA_LMM                                  } from '../../../modules/local/gemma/lmm/main'
include { EMMAX_KIN                                  } from '../../../modules/local/emmax/kin/main'
include { EMMAX_ASSOC                                } from '../../../modules/local/emmax/assoc/main'
include { TASSEL_MLM                                 } from '../../../modules/local/tassel/mlm/main'
include { RMVP_GWAS                                  } from '../../../modules/local/rmvp/gwas/main'
include { OMIGA_GWAS                                 } from '../../../modules/local/omiga/gwas/main'

workflow GWAS_BENCHMARK_PARALLEL {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_vcf         // channel: [ meta, vcf ] (empty → convert from bed)
    ch_phenotype   // channel: [ meta, phe ] plink-style FID IID trait header
    ch_covariates  // channel: [ meta, covar ] plink-style FID IID cov... (may be empty)

    main:
    // Local / nf-core plink modules publish versions via topic("versions").
    // Do not mix those tuples into Path-based ch_versions (breaks softwareVersionsToYAML).
    ch_versions = channel.empty()

    def engines = (params.gwas_benchmark_engines ?: 'gemma,emmax,tassel,rmvp,omiga')
        .tokenize(',')
        .collect { token -> token.trim().toLowerCase() }
        .findAll { token -> token }

    // Align optional covariates to phenotype meta (empty path list if absent)
    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _phe, cov -> [meta, cov ?: []] }

    PHENOCOVAR_ADAPT(ch_phenotype, ch_cov_aligned)

    GENOTYPE_TO_GWAS_FORMATS(ch_plink, ch_vcf)

    ch_gemma_assoc    = channel.empty()
    ch_emmax_assoc    = channel.empty()
    ch_tassel_results = channel.empty()
    ch_rmvp_results   = channel.empty()
    ch_omiga_gwas     = channel.empty()
    ch_raw_assoc      = channel.empty()

    // --- GEMMA: relatedness → LMM ---
    if (engines.contains('gemma')) {
        GEMMA_RELATEDNESS(GENOTYPE_TO_GWAS_FORMATS.out.bed)

        ch_gemma_cov = PHENOCOVAR_ADAPT.out.gemma_pheno
            .join(PHENOCOVAR_ADAPT.out.gemma_covar, remainder: true)
            .map { meta, _phe, cov -> [meta, cov ?: []] }

        GEMMA_LMM(
            GENOTYPE_TO_GWAS_FORMATS.out.bed,
            PHENOCOVAR_ADAPT.out.gemma_pheno,
            GEMMA_RELATEDNESS.out.relatedness,
            ch_gemma_cov
        )
        ch_gemma_assoc = GEMMA_LMM.out.assoc
        ch_raw_assoc = ch_raw_assoc.mix(ch_gemma_assoc.map { meta, f -> [meta, f, 'gemma'] })
    }

    // --- EMMAX: kinship → assoc ---
    if (engines.contains('emmax')) {
        EMMAX_KIN(GENOTYPE_TO_GWAS_FORMATS.out.tped)

        ch_emmax_cov = PHENOCOVAR_ADAPT.out.emmax_pheno
            .join(PHENOCOVAR_ADAPT.out.emmax_covar, remainder: true)
            .map { meta, _phe, cov -> [meta, cov ?: []] }

        EMMAX_ASSOC(
            GENOTYPE_TO_GWAS_FORMATS.out.tped,
            PHENOCOVAR_ADAPT.out.emmax_pheno,
            EMMAX_KIN.out.kinship,
            ch_emmax_cov
        )
        ch_emmax_assoc = EMMAX_ASSOC.out.assoc
        ch_raw_assoc = ch_raw_assoc.mix(ch_emmax_assoc.map { meta, f -> [meta, f, 'emmax'] })
    }

    // --- TASSEL / rMVP / OmiGA ---
    if (engines.contains('tassel')) {
        ch_tassel_kin = GENOTYPE_TO_GWAS_FORMATS.out.vcf
            .map { meta, _vcf -> [meta, []] }

        TASSEL_MLM(
            GENOTYPE_TO_GWAS_FORMATS.out.vcf,
            PHENOCOVAR_ADAPT.out.tassel_pheno,
            ch_tassel_kin
        )
        ch_tassel_results = TASSEL_MLM.out.results
        ch_raw_assoc = ch_raw_assoc.mix(ch_tassel_results.map { meta, f -> [meta, f, 'tassel'] })
    }

    if (engines.contains('rmvp')) {
        RMVP_GWAS(
            GENOTYPE_TO_GWAS_FORMATS.out.vcf,
            PHENOCOVAR_ADAPT.out.rmvp_pheno
        )
        ch_rmvp_results = RMVP_GWAS.out.results
        ch_raw_assoc = ch_raw_assoc.mix(ch_rmvp_results.map { meta, f -> [meta, f, 'rmvp'] })
    }

    if (engines.contains('omiga')) {
        ch_omiga_cov = PHENOCOVAR_ADAPT.out.omiga_pheno
            .join(PHENOCOVAR_ADAPT.out.omiga_covar, remainder: true)
            .map { meta, _phe, cov -> [meta, cov ?: []] }

        OMIGA_GWAS(
            GENOTYPE_TO_GWAS_FORMATS.out.bed,
            PHENOCOVAR_ADAPT.out.omiga_pheno,
            ch_omiga_cov
        )
        ch_omiga_gwas = OMIGA_GWAS.out.gwas
        ch_raw_assoc = ch_raw_assoc.mix(ch_omiga_gwas.map { meta, f -> [meta, f, 'omiga'] })
    }

    if (params.gwas_benchmark_standardize) {
        ASSOC_STANDARDIZE(ch_raw_assoc)
        ch_standardized = ASSOC_STANDARDIZE.out.standardized
    } else {
        ch_standardized = channel.empty()
    }

    emit:
    gemma_assoc    = ch_gemma_assoc
    emmax_assoc    = ch_emmax_assoc
    tassel_results = ch_tassel_results
    rmvp_results   = ch_rmvp_results
    omiga_gwas     = ch_omiga_gwas
    standardized   = ch_standardized
    raw_assoc      = ch_raw_assoc
    versions       = ch_versions
}
