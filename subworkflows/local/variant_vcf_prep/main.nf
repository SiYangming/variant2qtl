//
// Chain annotate → phase → impute on one VCF (skippable steps).
// Enable with params.run_vcf_prep (default false).
// Downstream steps consume the previous step's VCF, not a parallel copy of the input.
//

include { VARIANT_ANNOTATE } from '../variant_annotate/main'
include { VARIANT_PHASE    } from '../variant_phase/main'
include { VARIANT_IMPUTE   } from '../variant_impute/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// VARIANT_ANNOTATE.out.versions
// VARIANT_PHASE.out.versions
// VARIANT_IMPUTE.out.versions


workflow VARIANT_VCF_PREP {
    take:
    ch_vcf          // channel: [ meta, vcf, tbi ]
    ch_fasta        // channel: [ meta, fasta ] (dummy empty path ok)
    ch_vep_cache    // channel: [ meta, cache ]
    ch_snpeff_cache // channel: [ meta, cache ]
    ch_phase_ref    // channel: [ meta, vcf, tbi ] (may be empty)
    ch_phase_map    // channel: [ meta, map ] (may be empty)
    ch_impute_panel // channel: [ meta, panel, tbi ]
    ch_impute_map   // channel: [ meta, map ] (may be empty)

    main:
    ch_versions = channel.empty()
    def skip_ann = params.vcf_prep_skip_annotate
    def skip_phase = params.vcf_prep_skip_phase
    def skip_impute = params.vcf_prep_skip_impute

    ch_work = ch_vcf

    if (!skip_ann) {
        VARIANT_ANNOTATE(
            ch_work,
            ch_fasta,
            ch_vep_cache,
            ch_snpeff_cache
        )
        ch_versions = ch_versions.mix(VARIANT_ANNOTATE.out.versions)
        ch_work = VARIANT_ANNOTATE.out.vcf.map { meta, vcf -> [meta, vcf, []] }
    }

    if (!skip_phase) {
        VARIANT_PHASE(
            ch_work,
            ch_phase_ref,
            ch_phase_map
        )
        ch_versions = ch_versions.mix(VARIANT_PHASE.out.versions)
        ch_work = VARIANT_PHASE.out.phased.map { meta, vcf -> [meta, vcf, []] }
    }

    ch_imputed = channel.empty()
    if (!skip_impute) {
        VARIANT_IMPUTE(
            ch_work,
            ch_impute_panel,
            ch_impute_map
        )
        ch_versions = ch_versions.mix(VARIANT_IMPUTE.out.versions)
        ch_imputed = VARIANT_IMPUTE.out.imputed
    }

    emit:
    vcf      = ch_work
    imputed  = ch_imputed
    versions = ch_versions
}
