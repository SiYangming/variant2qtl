//
// FINEMAP / CAVIAR / DAP-G style PIPs from sumstats (identity LD).
// Enable with params.run_finemap_extra (default false).
// Versions via topic("versions") on FINEMAP_EXTRA — do not mix into Path ch_versions.
//

include { FINEMAP_EXTRA } from '../../../modules/local/finemap/extra/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// FINEMAP_EXTRA.out.versions


workflow QTL_FINEMAP_EXTRA {
    take:
    ch_sumstats  // channel: [ meta, sumstats ]

    main:
    ch_versions = channel.empty()

    FINEMAP_EXTRA(ch_sumstats)

    emit:
    pip            = FINEMAP_EXTRA.out.pip
    credible_sets  = FINEMAP_EXTRA.out.credible_sets
    outdir         = FINEMAP_EXTRA.out.outdir
    versions       = ch_versions
}
