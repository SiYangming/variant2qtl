//
// FUSION/PrediXcan-style TWAS from gene weights and GWAS sumstats.
// Enable with params.run_twas (default false).
// Versions via topic("versions") on TWAS_FUSION — do not mix into Path ch_versions.
//

include { TWAS_FUSION } from '../../../modules/local/twas/fusion/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// TWAS_FUSION.out.versions


workflow QTL_TWAS {
    take:
    ch_weights  // channel: [ meta, weights ]
    ch_gwas     // channel: [ meta, gwas ]

    main:
    ch_versions = channel.empty()

    ch_pairs = ch_weights
        .combine(ch_gwas)
        .map { meta, weights, _gwas_meta, gwas ->
            [meta, weights, gwas]
        }

    TWAS_FUSION(
        ch_pairs.map { meta, weights, _gwas -> [meta, weights] },
        ch_pairs.map { meta, _weights, gwas -> [meta, gwas] }
    )

    emit:
    twas     = TWAS_FUSION.out.twas
    outdir   = TWAS_FUSION.out.outdir
    versions = ch_versions
}
