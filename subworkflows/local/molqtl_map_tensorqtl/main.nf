//
// tensorQTL cis-molQTL mapping (CPU-friendly default; GPU via module ext.args).
// Enable from the main workflow with params.run_tensorqtl_cis (default false).
// Versions via topic("versions") on TENSORQTL_CIS — do not mix into Path ch_versions.
//

include { TENSORQTL_CIS } from '../../../modules/local/tensorqtl/cis/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// TENSORQTL_CIS.out.versions


workflow MOLQTL_MAP_TENSORQTL {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_phenotype   // channel: [ meta, phenotype BED ]
    ch_covariates  // channel: [ meta, covar ] (may be empty path list)

    main:
    ch_versions = channel.empty()

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    TENSORQTL_CIS(
        ch_plink,
        ch_phenotype,
        ch_cov_aligned
    )

    emit:
    cis_qtl  = TENSORQTL_CIS.out.cis_qtl
    outdir   = TENSORQTL_CIS.out.outdir
    versions = ch_versions
}
