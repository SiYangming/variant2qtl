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

    ch_vcf_keyed = ch_vcf.map { meta, vcf, tbi ->
        [meta.group_id ?: meta.id, meta, vcf, tbi]
    }

    ch_ref_aligned = ch_vcf_keyed
        .join(
            ch_ref.map { meta, ref_vcf, ref_tbi -> [meta.id, ref_vcf, ref_tbi] },
            remainder: true
        )
        .map { row ->
            def meta = row[1]
            def ref_vcf = (row.size() > 4 && row[4]) ? row[4] : []
            def ref_tbi = (row.size() > 5 && row[5]) ? row[5] : []
            [meta, ref_vcf, ref_tbi]
        }

    ch_map_aligned = ch_vcf_keyed
        .join(
            ch_map.map { meta, gmap -> [meta.id, gmap] },
            remainder: true
        )
        .map { row ->
            def meta = row[1]
            def gmap = (row.size() > 4 && row[4]) ? row[4] : []
            [meta, gmap]
        }

    ch_phase_in = ch_vcf
        .join(ch_ref_aligned)
        .join(ch_map_aligned)
        .map { meta, vcf, tbi, ref_vcf, ref_tbi, gmap ->
            def region = (meta.region ?: params.phase_region ?: '1').toString()
            [meta, vcf, tbi, [], region, ref_vcf, ref_tbi, [], [], gmap]
        }

    SHAPEIT5_PHASECOMMON(ch_phase_in)

    emit:
    phased   = SHAPEIT5_PHASECOMMON.out.phased_variant
    versions = ch_versions
}
