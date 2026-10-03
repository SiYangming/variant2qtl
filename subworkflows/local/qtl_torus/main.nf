//
// Annotation enrichment prior from QTL z-scores.
// Enable with params.run_torus (default false).
// Versions via topic("versions") on TORUS_ENRICH — do not mix into Path ch_versions.
//

include { TORUS_ENRICH } from '../../../modules/local/torus/enrich/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// TORUS_ENRICH.out.versions


workflow QTL_TORUS {
    take:
    ch_annot  // channel: [ meta, annot ]

    main:
    ch_versions = channel.empty()

    TORUS_ENRICH(ch_annot)

    emit:
    prior    = TORUS_ENRICH.out.prior
    enrich   = TORUS_ENRICH.out.enrich
    outdir   = TORUS_ENRICH.out.outdir
    versions = ch_versions
}
