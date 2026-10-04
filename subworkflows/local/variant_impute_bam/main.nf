//
// Optional GLIMPSE2 imputation (chunk → phase → ligate).
// Enable with params.run_impute_bam (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { BAM_VCF_IMPUTE_GLIMPSE2 } from '../../../subworkflows/nf-core/bam_vcf_impute_glimpse2/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// BAM_VCF_IMPUTE_GLIMPSE2.out.versions
// GLIMPSE2_CHUNK.out.versions
// GLIMPSE2_SPLITREFERENCE.out.versions
// GLIMPSE2_PHASE.out.versions
// GLIMPSE2_LIGATE.out.versions
// BCFTOOLS_INDEX.out.versions


workflow VARIANT_IMPUTE_BAM {
    take:
    ch_input  // channel: [ meta, bam_or_vcf, index ]
    ch_panel  // channel: [ meta, panel_vcf, tbi ]
    ch_fasta  // channel: [ meta, fasta, fai ] (dummy empty paths ok)
    ch_map    // channel: [ meta, map ] (may be empty)

    main:
    ch_versions = channel.empty()
    def region_raw = (params.impute_bam_region ?: params.impute_region ?: '1').toString()
    def region = region_raw.contains(':') ? region_raw : "${region_raw}:1-999999999"
    def do_chunk = params.impute_bam_chunk ? true : false
    def chunk_model = params.impute_bam_chunk_model ?: 'sequential'
    def do_split = params.impute_bam_split_ref ? true : false

    ch_glimpse_input = ch_input.map { meta, input, idx ->
        [meta, input, idx ?: [], [], []]
    }

    ch_ref = ch_panel.map { meta, panel, tbi ->
        [meta, panel, tbi ?: [], region]
    }

    ch_map_aligned = ch_panel
        .join(ch_map, remainder: true)
        .map { meta, _panel, _tbi, gmap -> [meta, gmap ?: []] }

    ch_chunks = do_chunk
        ? channel.empty()
        : ch_panel.map { meta, _panel, _tbi -> [meta, region, region] }

    BAM_VCF_IMPUTE_GLIMPSE2(
        ch_glimpse_input,
        ch_ref,
        ch_chunks,
        ch_map_aligned,
        ch_fasta,
        do_chunk,
        chunk_model,
        do_split
    )

    emit:
    imputed  = BAM_VCF_IMPUTE_GLIMPSE2.out.vcf_index.map { meta, vcf, _idx -> [meta, vcf] }
    versions = ch_versions
}
