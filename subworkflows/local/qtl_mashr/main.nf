//
// Multi-condition mashr shrinkage.
// Enable with params.run_mashr (default false).
// Versions via topic("versions") on MASHR_FIT — do not mix into Path ch_versions.
//

include { MASHR_FIT } from '../../../modules/local/mashr/fit/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// MASHR_FIT.out.versions


workflow QTL_MASHR {
    take:
    ch_effects  // channel: [ meta, effects ]

    main:
    ch_versions = channel.empty()

    MASHR_FIT(ch_effects)

    emit:
    mash     = MASHR_FIT.out.mash
    lfsr     = MASHR_FIT.out.lfsr
    outdir   = MASHR_FIT.out.outdir
    versions = ch_versions
}
