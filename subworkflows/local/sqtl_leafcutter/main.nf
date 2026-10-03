//
// LeafCutter intron counts → FastQTL BED + phenotype_group.
// Enable with params.run_sqtl_leafcutter (default false).
// Versions via topic("versions") on LEAFCUTTER_PREPARE — do not mix into Path ch_versions.
//

include { LEAFCUTTER_PREPARE } from '../../../modules/local/leafcutter/prepare/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// LEAFCUTTER_PREPARE.out.versions


workflow SQTL_LEAFCUTTER {
    take:
    ch_counts  // channel: [ meta, perind counts ]
    ch_genes   // channel: [ meta, gene BED ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_genes_aligned = ch_counts
        .join(ch_genes, remainder: true)
        .map { meta, _counts, genes -> [meta, genes ?: []] }

    LEAFCUTTER_PREPARE(
        ch_counts,
        ch_genes_aligned
    )

    emit:
    bed              = LEAFCUTTER_PREPARE.out.bed
    phenotype_group  = LEAFCUTTER_PREPARE.out.phenotype_group
    outdir           = LEAFCUTTER_PREPARE.out.outdir
    versions         = ch_versions
}
