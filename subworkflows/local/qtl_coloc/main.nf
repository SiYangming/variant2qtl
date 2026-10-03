//
// QTL–GWAS colocalisation via coloc.abf.
// Enable from the main workflow with params.run_coloc (default false).
// Versions via topic("versions") on COLOC_ABF — do not mix into Path ch_versions.
//

include { COLOC_ABF } from '../../../modules/local/coloc/abf/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// COLOC_ABF.out.versions


workflow QTL_COLOC {
    take:
    ch_qtl   // channel: [ meta, qtl_sumstats ]
    ch_gwas  // channel: [ meta, gwas_sumstats ]

    main:
    ch_versions = channel.empty()

    // Broadcast GWAS sumstats onto each QTL file (meta ids often differ).
    // combine() of two [meta, path] channels flattens to four arguments.
    ch_pairs = ch_qtl
        .combine(ch_gwas)
        .map { meta, qtl, _gwas_meta, gwas ->
            [meta, qtl, gwas]
        }

    COLOC_ABF(
        ch_pairs.map { meta, qtl, _gwas -> [meta, qtl] },
        ch_pairs.map { meta, _qtl, gwas -> [meta, gwas] }
    )

    emit:
    summary  = COLOC_ABF.out.summary
    abf      = COLOC_ABF.out.abf
    coloc_log = COLOC_ABF.out.log
    outdir   = COLOC_ABF.out.outdir
    versions = ch_versions
}
