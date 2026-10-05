//
// MatrixEQTL cis-molQTL mapping.
// Enable from the main workflow with params.run_matrixeqtl_cis (default false).
// Versions via topic("versions") on MATRIXEQTL_CIS — do not mix into Path ch_versions.
//

include { MATRIXEQTL_CIS } from '../../../modules/local/matrixeqtl/cis/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// MATRIXEQTL_CIS.out.versions


workflow MOLQTL_MAP_MATRIXEQTL {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_phenotype   // channel: [ meta, phenotype BED ]
    ch_covariates  // channel: [ meta, covar ] (may be empty path list)

    main:
    ch_versions = channel.empty()

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    MATRIXEQTL_CIS(
        ch_plink,
        ch_phenotype,
        ch_cov_aligned
    )

    emit:
    cis_qtl  = MATRIXEQTL_CIS.out.cis_qtl
    outdir   = MATRIXEQTL_CIS.out.outdir
    versions = ch_versions
}
