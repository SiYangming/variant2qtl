//
// LD-score regression style heritability.
// Enable with params.run_ldsc (default false).
// Versions via topic("versions") on LDSC_H2 — do not mix into Path ch_versions.
//

include { LDSC_H2 } from '../../../modules/local/ldsc/h2/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// LDSC_H2.out.versions


workflow QTL_LDSC {
    take:
    ch_sumstats  // channel: [ meta, sumstats ]
    ch_annot     // channel: [ meta, annot ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_annot_aligned = ch_sumstats
        .join(ch_annot, remainder: true)
        .map { meta, _sumstats, annot -> [meta, annot ?: []] }

    LDSC_H2(
        ch_sumstats,
        ch_annot_aligned
    )

    emit:
    h2           = LDSC_H2.out.h2
    partitioned  = LDSC_H2.out.partitioned
    outdir       = LDSC_H2.out.outdir
    versions     = ch_versions
}
