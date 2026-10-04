//
// FUSION/PrediXcan/MetaXcan-style TWAS from gene weights and GWAS sumstats.
// Accepts long TSV weights or PredictDB (.db) SQLite via TWAS_PREDICTDB.
// Enable with params.run_twas (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { TWAS_PREDICTDB } from '../../../modules/local/twas/predictdb/main'
include { TWAS_FUSION    } from '../../../modules/local/twas/fusion/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// TWAS_PREDICTDB.out.versions
// TWAS_FUSION.out.versions


workflow QTL_TWAS {
    take:
    ch_weights  // channel: [ meta, weights (.tsv or PredictDB .db) ]
    ch_gwas     // channel: [ meta, gwas ]
    ch_ld       // channel: [ meta, ld ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_w_branched = ch_weights.branch { _meta, weights ->
        db: weights.name.toLowerCase().endsWith('.db')
        tsv: true
    }

    TWAS_PREDICTDB(ch_w_branched.db)
    ch_weights_tsv = ch_w_branched.tsv.mix(TWAS_PREDICTDB.out.weights)

    ch_pairs = ch_weights_tsv
        .combine(ch_gwas)
        .map { meta, weights, _gwas_meta, gwas ->
            [meta, weights, gwas]
        }

    ch_ld_aligned = ch_pairs
        .join(ch_ld, remainder: true)
        .map { meta, _weights, _gwas, ld -> [meta, ld ?: []] }

    TWAS_FUSION(
        ch_pairs.map { meta, weights, _gwas -> [meta, weights] },
        ch_pairs.map { meta, _weights, gwas -> [meta, gwas] },
        ch_ld_aligned
    )

    emit:
    twas     = TWAS_FUSION.out.twas
    outdir   = TWAS_FUSION.out.outdir
    versions = ch_versions
}
