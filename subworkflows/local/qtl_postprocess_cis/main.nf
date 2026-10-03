//
// Harmonise cis-QTL tables and compute BH q-values.
// Enable from the main workflow with params.run_qtl_postprocess (default false).
// Versions via topic("versions") on QTL_POSTPROCESS — do not mix into Path ch_versions.
//

include { QTL_POSTPROCESS } from '../../../modules/local/utils/qtl_postprocess/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// QTL_POSTPROCESS.out.versions


workflow QTL_POSTPROCESS_CIS {
    take:
    ch_cis  // channel: [ meta, cis_qtl ]

    main:
    ch_versions = channel.empty()

    QTL_POSTPROCESS(ch_cis)

    emit:
    standardized = QTL_POSTPROCESS.out.standardized
    versions     = ch_versions
}
