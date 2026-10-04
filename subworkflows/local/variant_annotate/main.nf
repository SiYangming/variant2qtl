//
// Optional VCF annotation via nf-core Ensembl VEP and/or snpEff modules.
// Enable with params.run_annotate (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { SNPEFF_SNPEFF  } from '../../../modules/nf-core/snpeff/snpeff/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// ENSEMBLVEP_VEP.out.versions
// SNPEFF_SNPEFF.out.versions


workflow VARIANT_ANNOTATE {
    take:
    ch_vcf          // channel: [ meta, vcf, tbi ]
    ch_fasta        // channel: [ meta, fasta ] (use [[id], []] when absent)
    ch_vep_cache    // channel: [ meta, cache ] (use [[id], []] when absent)
    ch_snpeff_cache // channel: [ meta, cache ] (use [[id], []] when absent)

    main:
    ch_versions = channel.empty()

    def tools = (params.annotate_tools ?: 'snpeff')
        .tokenize(',')
        .collect { tool -> tool.trim().toLowerCase() }

    ch_current = ch_vcf.map { meta, vcf, _tbi -> [meta, vcf] }

    if (tools.contains('ensemblvep')) {
        def genome = params.annotate_vep_genome ?: (params.genome ?: 'GRCh38')
        def species = params.annotate_vep_species ?: 'homo_sapiens'
        def cache_version = params.annotate_vep_cache_version ?: '110'

        ENSEMBLVEP_VEP(
            ch_current.map { meta, vcf -> [meta, vcf, []] },
            genome,
            species,
            cache_version,
            ch_vep_cache,
            ch_fasta,
            []
        )
        ch_current = ENSEMBLVEP_VEP.out.vcf
    }

    if (tools.contains('snpeff')) {
        def db = params.annotate_snpeff_db ?: 'GRCh38.99'
        SNPEFF_SNPEFF(
            ch_current,
            db,
            ch_snpeff_cache
        )
        ch_current = SNPEFF_SNPEFF.out.vcf
    }

    emit:
    vcf      = ch_current
    versions = ch_versions
}
