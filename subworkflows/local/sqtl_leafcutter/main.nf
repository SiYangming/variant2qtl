//
// LeafCutter intron counts → FastQTL BED + phenotype_group.
// Optional BAM → regtools junctions extract → leafcutter clusterregtools.
// Enable with params.run_sqtl_leafcutter and/or params.run_leafcutter_cluster (default false).
// Mix Path versions from REGTOOLS only — do not mix topic versions into ch_versions.
//

include { REGTOOLS_JUNCTIONSEXTRACT    } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTERREGTOOLS   } from '../../../modules/nf-core/leafcutter/clusterregtools/main'
include { LEAFCUTTER_PREPARE           } from '../../../modules/local/leafcutter/prepare/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// LEAFCUTTER_PREPARE.out.versions, LEAFCUTTER_CLUSTERREGTOOLS.out.versions
// REGTOOLS_JUNCTIONSEXTRACT.out.versions


workflow SQTL_LEAFCUTTER {
    take:
    ch_counts  // channel: [ meta, perind counts ] (empty when clustering from BAM)
    ch_genes   // channel: [ meta, gene BED ] (may be empty)
    ch_bam     // channel: [ meta, bam, bai ] (empty when using counts)

    main:
    ch_versions = channel.empty()
    ch_counts_use = ch_counts

    if (params.sqtl_bam || params.run_leafcutter_cluster) {
        REGTOOLS_JUNCTIONSEXTRACT(
            ch_bam,
            params.sqtl_junc_strand ?: 'XS'
        )
        ch_versions = ch_versions.mix(REGTOOLS_JUNCTIONSEXTRACT.out.versions)
        LEAFCUTTER_CLUSTERREGTOOLS(
            REGTOOLS_JUNCTIONSEXTRACT.out.junc.map { meta, junc -> [meta, junc] }
        )
        ch_counts_use = LEAFCUTTER_CLUSTERREGTOOLS.out.counts
    }

    ch_genes_aligned = ch_counts_use
        .join(ch_genes, remainder: true)
        .map { meta, _counts, genes -> [meta, genes ?: []] }

    LEAFCUTTER_PREPARE(
        ch_counts_use,
        ch_genes_aligned
    )

    emit:
    bed              = LEAFCUTTER_PREPARE.out.bed
    phenotype_group  = LEAFCUTTER_PREPARE.out.phenotype_group
    outdir           = LEAFCUTTER_PREPARE.out.outdir
    versions         = ch_versions
}
