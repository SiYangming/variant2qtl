//
// Optional VCF phasing via nf-core SHAPEIT5 phase_common.
// Enable with params.run_phase (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { SHAPEIT5_PHASECOMMON } from '../../../modules/nf-core/shapeit5/phasecommon/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// SHAPEIT5_PHASECOMMON.out.versions


workflow VARIANT_PHASE {
    take:
    ch_vcf // channel: [ meta, vcf, tbi ]
    ch_ref // channel: [ meta, vcf, tbi ] (optional; may be empty)
    ch_map // channel: [ meta, map ] (optional; may be empty)

    main:
    ch_versions = channel.empty()

    def region = (params.phase_region ?: '1').toString()

    ch_ref_aligned = ch_vcf
        .join(ch_ref, remainder: true)
        .map { meta, _vcf, _tbi, ref_vcf, ref_tbi ->
            [meta, ref_vcf ?: [], ref_tbi ?: []]
        }

    ch_map_aligned = ch_vcf
        .join(ch_map, remainder: true)
        .map { meta, _vcf, _tbi, gmap ->
            [meta, gmap ?: []]
        }

    ch_phase_in = ch_vcf
        .join(ch_ref_aligned)
        .join(ch_map_aligned)
        .map { meta, vcf, tbi, ref_vcf, ref_tbi, gmap ->
            [meta, vcf, tbi, [], region, ref_vcf, ref_tbi, [], [], gmap]
        }

    SHAPEIT5_PHASECOMMON(ch_phase_in)

    emit:
    phased   = SHAPEIT5_PHASECOMMON.out.phased_variant
    versions = ch_versions
}
