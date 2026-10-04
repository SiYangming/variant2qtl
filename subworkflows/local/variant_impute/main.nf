//
// Optional genotype imputation via nf-core Beagle5 / Minimac4 / GLIMPSE phase.
// Enable with params.run_impute (default false). Engine via params.impute_engine.
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { BEAGLE5_BEAGLE  } from '../../../modules/nf-core/beagle5/beagle/main'
include { MINIMAC4_IMPUTE } from '../../../modules/nf-core/minimac4/impute/main'
include { GLIMPSE_PHASE   } from '../../../modules/nf-core/glimpse/phase/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// BEAGLE5_BEAGLE.out.versions
// MINIMAC4_IMPUTE.out.versions
// GLIMPSE_PHASE.out.versions


workflow VARIANT_IMPUTE {
    take:
    ch_vcf    // channel: [ meta, vcf, tbi ]
    ch_panel  // channel: [ meta, panel_vcf_or_msav, tbi_or_empty ]
    ch_map    // channel: [ meta, map ] (optional; may be empty)

    main:
    ch_versions = channel.empty()
    ch_imputed = channel.empty()

    def engine = (params.impute_engine ?: 'beagle5').toString().trim().toLowerCase()
    def region = (params.impute_region ?: '1').toString()

    ch_map_aligned = ch_vcf
        .join(ch_map, remainder: true)
        .map { meta, _vcf, _tbi, gmap -> [meta, gmap ?: []] }

    ch_joined = ch_vcf
        .join(ch_panel)
        .join(ch_map_aligned)
        .map { meta, vcf, tbi, panel, panel_tbi, gmap ->
            [meta, vcf, tbi, panel, panel_tbi ?: [], gmap]
        }

    if (engine == 'beagle5' || engine == 'beagle') {
        BEAGLE5_BEAGLE(
            ch_joined.map { meta, vcf, tbi, panel, panel_tbi, gmap ->
                [meta, vcf, tbi, panel, panel_tbi, gmap, [], [], region]
            }
        )
        ch_imputed = BEAGLE5_BEAGLE.out.vcf
    }
    else if (engine == 'minimac4' || engine == 'minimac') {
        MINIMAC4_IMPUTE(
            ch_joined.map { meta, vcf, tbi, panel, _panel_tbi, gmap ->
                [meta, vcf, tbi, panel, [], [], gmap, region]
            }
        )
        ch_imputed = MINIMAC4_IMPUTE.out.vcf
    }
    else if (engine == 'glimpse' || engine == 'glimpse1') {
        GLIMPSE_PHASE(
            ch_joined.map { meta, vcf, tbi, panel, panel_tbi, gmap ->
                [meta, vcf, tbi, [], region, region, panel, panel_tbi, gmap]
            }
        )
        ch_imputed = GLIMPSE_PHASE.out.phased_variants
    }
    else {
        error("Unknown impute_engine '${engine}'. Use beagle5, minimac4, or glimpse.")
    }

    emit:
    imputed  = ch_imputed
    versions = ch_versions
}
