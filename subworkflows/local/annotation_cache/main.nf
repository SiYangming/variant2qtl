//
// Optional Ensembl VEP / snpEff cache download.
// Enable with params.run_cache (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { SNPEFF_DOWNLOAD     } from '../../../modules/nf-core/snpeff/download/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// ENSEMBLVEP_DOWNLOAD.out.versions
// SNPEFF_DOWNLOAD.out.versions


workflow ANNOTATION_CACHE {
    take:
    ch_vep_info   // channel: [ meta, assembly, species, cache_version ] (may be empty)
    ch_snpeff_info // channel: [ meta, snpeff_db ] (may be empty)

    main:
    ch_versions = channel.empty()
    ch_vep_cache = channel.empty()
    ch_snpeff_cache = channel.empty()

    def tools = (params.cache_tools ?: 'snpeff,ensemblvep')
        .tokenize(',')
        .collect { tool -> tool.trim().toLowerCase() }

    if (tools.contains('ensemblvep')) {
        ENSEMBLVEP_DOWNLOAD(
            ch_vep_info,
            params.cache_vep_preflight ?: false
        )
        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
    }

    if (tools.contains('snpeff')) {
        SNPEFF_DOWNLOAD(ch_snpeff_info)
        ch_snpeff_cache = SNPEFF_DOWNLOAD.out.cache
    }

    emit:
    vep_cache    = ch_vep_cache
    snpeff_cache = ch_snpeff_cache
    versions     = ch_versions
}
