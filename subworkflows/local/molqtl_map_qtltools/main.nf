//
// QTLtools cis-molQTL mapping.
// Enable from the main workflow with params.run_qtltools_cis (default false).
// Versions via topic("versions") on QTLTOOLS_CIS — do not mix into Path ch_versions.
//

include { QTLTOOLS_CIS } from '../../../modules/local/qtltools/cis/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// QTLTOOLS_CIS.out.versions


workflow MOLQTL_MAP_QTLTOOLS {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_phenotype   // channel: [ meta, phenotype BED ]
    ch_covariates  // channel: [ meta, covar ] (may be empty path list)

    main:
    ch_versions = channel.empty()

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    QTLTOOLS_CIS(
        ch_plink,
        ch_phenotype,
        ch_cov_aligned
    )

    emit:
    cis_qtl  = QTLTOOLS_CIS.out.cis_qtl
    outdir   = QTLTOOLS_CIS.out.outdir
    versions = ch_versions
}
