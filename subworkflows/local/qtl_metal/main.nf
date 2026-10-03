//
// Inverse-variance meta-analysis of multi-cohort sumstats.
// Enable with params.run_metal (default false).
// Versions via topic("versions") on METAL_IVW — do not mix into Path ch_versions.
//

include { METAL_IVW } from '../../../modules/local/metal/ivw/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// METAL_IVW.out.versions


workflow QTL_METAL {
    take:
    ch_cohorts  // channel: [ meta, cohorts ]

    main:
    ch_versions = channel.empty()

    METAL_IVW(ch_cohorts)

    emit:
    metal    = METAL_IVW.out.metal
    outdir   = METAL_IVW.out.outdir
    versions = ch_versions
}
