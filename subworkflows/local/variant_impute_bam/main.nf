//
// Optional GLIMPSE2 imputation from BAM/CRAM or genotype-likelihood VCF.
// Enable with params.run_impute_bam (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { GLIMPSE2_PHASE } from '../../../modules/nf-core/glimpse2/phase/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// GLIMPSE2_PHASE.out.versions


workflow VARIANT_IMPUTE_BAM {
    take:
    ch_input  // channel: [ meta, bam_or_vcf, index ]
    ch_panel  // channel: [ meta, panel_vcf, tbi ]
    ch_fasta  // channel: [ meta, fasta, fai ] (dummy empty paths ok)
    ch_map    // channel: [ meta, map ] (may be empty)

    main:
    ch_versions = channel.empty()
    def region = (params.impute_bam_region ?: params.impute_region ?: '1').toString()
    def suffix = params.impute_bam_suffix ?: 'vcf.gz'

    ch_map_aligned = ch_input
        .join(ch_map, remainder: true)
        .map { meta, _input, _idx, gmap -> [meta, gmap ?: []] }

    ch_joined = ch_input
        .join(ch_panel)
        .join(ch_map_aligned)
        .map { meta, input, idx, panel, panel_tbi, gmap ->
            [meta, input, idx ?: [], [], [], region, region, panel, panel_tbi ?: [], gmap]
        }

    GLIMPSE2_PHASE(
        ch_joined,
        ch_fasta,
        suffix
    )

    emit:
    imputed  = GLIMPSE2_PHASE.out.phased_variants
    versions = ch_versions
}
