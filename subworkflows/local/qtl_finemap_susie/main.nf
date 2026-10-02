//
// SuSiE fine-mapping from association / QTL summary statistics.
// Enable from the main workflow with params.run_finemap_susie (default false).
// Versions via topic("versions") on SUSIE_FINEMAP — do not mix into Path ch_versions.
//

include { SUSIE_FINEMAP } from '../../../modules/local/susie/finemap/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// SUSIE_FINEMAP.out.versions


workflow QTL_FINEMAP_SUSIE {
    take:
    ch_sumstats  // channel: [ meta, sumstats ]
    ch_ld        // channel: [ meta, ld ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_ld_aligned = ch_sumstats
        .join(ch_ld, remainder: true)
        .map { meta, _sumstats, ld -> [meta, ld ?: []] }

    SUSIE_FINEMAP(
        ch_sumstats,
        ch_ld_aligned
    )

    emit:
    pip            = SUSIE_FINEMAP.out.pip
    credible_sets  = SUSIE_FINEMAP.out.credible_sets
    susie_log      = SUSIE_FINEMAP.out.log
    outdir         = SUSIE_FINEMAP.out.outdir
    versions       = ch_versions
}
