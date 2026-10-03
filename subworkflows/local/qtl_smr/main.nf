//
// SMR and HEIDI from QTL + GWAS summary statistics.
// Enable with params.run_smr (default false).
// Versions via topic("versions") on SMR_HEIDI — do not mix into Path ch_versions.
//

include { SMR_HEIDI } from '../../../modules/local/smr/heidi/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// SMR_HEIDI.out.versions


workflow QTL_SMR {
    take:
    ch_qtl   // channel: [ meta, qtl_sumstats ]
    ch_gwas  // channel: [ meta, gwas_sumstats ]

    main:
    ch_versions = channel.empty()

    ch_pairs = ch_qtl
        .combine(ch_gwas)
        .map { meta, qtl, _gwas_meta, gwas ->
            [meta, qtl, gwas]
        }

    SMR_HEIDI(
        ch_pairs.map { meta, qtl, _gwas -> [meta, qtl] },
        ch_pairs.map { meta, _qtl, gwas -> [meta, gwas] }
    )

    emit:
    smr      = SMR_HEIDI.out.smr
    outdir   = SMR_HEIDI.out.outdir
    versions = ch_versions
}
