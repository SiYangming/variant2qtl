//
// PEER hidden covariates from a molecular phenotype BED.
// Enable with params.run_peer (default false).
// Versions via topic("versions") on PEER_FACTORS — do not mix into Path ch_versions.
//

include { PEER_FACTORS } from '../../../modules/local/peer/factors/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// PEER_FACTORS.out.versions


workflow COVARIATE_PEER {
    take:
    ch_phenotype   // channel: [ meta, FastQTL BED ]
    ch_covariates  // channel: [ meta, covariates ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    PEER_FACTORS(
        ch_phenotype,
        ch_cov_aligned
    )

    emit:
    cov_tensorqtl = PEER_FACTORS.out.cov_tensorqtl
    cov_omiga     = PEER_FACTORS.out.cov_omiga
    outdir        = PEER_FACTORS.out.outdir
    versions      = ch_versions
}
