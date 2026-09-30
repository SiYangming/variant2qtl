//
// OmiGA cis-molQTL mapping.
// Optionally fans genotypes through genotype_to_analysis_format (BGEN) before OMIGA_CIS.
// Enable from the main workflow with params.run_omiga_cis (default false).
// Versions via topic("versions") on OMIGA_CIS — do not mix into Path ch_versions.
//

include { OMIGA_CIS                     } from '../../../modules/local/omiga/cis/main'
include { GENOTYPE_TO_ANALYSIS_FORMAT   } from '../genotype_to_analysis_format/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// OMIGA_CIS.out.versions, GENOTYPE_TO_ANALYSIS_FORMAT.out.versions


workflow MOLQTL_MAP_OMIGA {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_phenotype   // channel: [ meta, phenotype BED/OPF ]
    ch_covariates  // channel: [ meta, covar ] (may be empty path list)
    ch_vcf         // channel: [ meta, vcf ] (optional; empty OK)

    main:
    ch_versions = channel.empty()

    GENOTYPE_TO_ANALYSIS_FORMAT(ch_plink, ch_vcf)
    ch_versions = ch_versions.mix(GENOTYPE_TO_ANALYSIS_FORMAT.out.versions)

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    OMIGA_CIS(
        GENOTYPE_TO_ANALYSIS_FORMAT.out.bed,
        ch_phenotype,
        ch_cov_aligned
    )

    emit:
    cis_qtl  = OMIGA_CIS.out.cis_qtl
    outdir   = OMIGA_CIS.out.outdir
    bed      = GENOTYPE_TO_ANALYSIS_FORMAT.out.bed
    bgen     = GENOTYPE_TO_ANALYSIS_FORMAT.out.bgen
    versions = ch_versions
}
