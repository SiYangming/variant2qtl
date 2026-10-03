//
// Multi-trait HyPrColoc-style clustering.
// Enable with params.run_hyprcoloc (default false).
// Versions via topic("versions") on HYPRCOLOC — do not mix into Path ch_versions.
//

include { HYPRCOLOC } from '../../../modules/local/coloc/hyprcoloc/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// HYPRCOLOC.out.versions


workflow QTL_HYPRCOLOC {
    take:
    ch_traits  // channel: [ meta, long-format trait sumstats ]

    main:
    ch_versions = channel.empty()

    HYPRCOLOC(ch_traits)

    emit:
    clusters = HYPRCOLOC.out.clusters
    pairs    = HYPRCOLOC.out.pairs
    outdir   = HYPRCOLOC.out.outdir
    versions = ch_versions
}
