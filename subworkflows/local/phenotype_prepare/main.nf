//
// Sample intersect, missingness filter, INV/quantile, FastQTL BED.
// Enable with params.run_phenotype_prepare (default false).
// Versions via topic("versions") on PHENOTYPE_PREPARE — do not mix into Path ch_versions.
//

include { PHENOTYPE_PREPARE } from '../../../modules/local/utils/phenotype_prepare/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// PHENOTYPE_PREPARE.out.versions


workflow PHENOTYPE_PREPARE_SWF {
    take:
    ch_matrix    // channel: [ meta, matrix ]
    ch_gene_map  // channel: [ meta, gene_bed ] (may be empty)
    ch_samples   // channel: [ meta, samples/fam ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_gene_aligned = ch_matrix
        .join(ch_gene_map, remainder: true)
        .map { meta, _matrix, gene_map -> [meta, gene_map ?: []] }

    ch_samples_aligned = ch_matrix
        .join(ch_samples, remainder: true)
        .map { meta, _matrix, samples -> [meta, samples ?: []] }

    PHENOTYPE_PREPARE(
        ch_matrix,
        ch_gene_aligned,
        ch_samples_aligned
    )

    emit:
    bed      = PHENOTYPE_PREPARE.out.bed
    outdir   = PHENOTYPE_PREPARE.out.outdir
    versions = ch_versions
}
