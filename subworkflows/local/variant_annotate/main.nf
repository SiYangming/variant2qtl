//
// Optional VCF annotation: vcfanno (resource INFO) → VEP/snpEff (scatter-gather) → filter.
// Enable with params.run_annotate (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { VCFANNO                            } from '../../../modules/nf-core/vcfanno/main'
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF     } from '../../../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP     } from '../../../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// VCFANNO.out.versions
// VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.versions
// ENSEMBLVEP_VEP.out.versions
// SNPEFF_SNPEFF.out.versions
// BCFTOOLS_PLUGINSCATTER.out.versions
// BCFTOOLS_CONCAT.out.versions
// BCFTOOLS_SORT.out.versions
// HTSLIB_BGZIPTABIX.out.versions
// VCF_FILTER_BCFTOOLS_ENSEMBLVEP.out.versions
// BCFTOOLS_VIEW.out.versions
// ENSEMBLVEP_FILTERVEP.out.versions


workflow VARIANT_ANNOTATE {
    take:
    ch_vcf          // channel: [ meta, vcf, tbi ]
    ch_fasta        // channel: [ meta, fasta ] (use [[id], []] when absent)
    ch_vep_cache    // channel: [ meta, cache ] (use [[id], []] when absent)
    ch_snpeff_cache // channel: [ meta, cache ] (use [[id], []] when absent)

    main:
    ch_versions = channel.empty()

    def tools = (params.annotate_tools ?: '')
        .tokenize(',')
        .collect { tool -> tool.trim().toLowerCase() }
        .findAll { tool -> tool in ['ensemblvep', 'snpeff'] }

    def run_vcfanno = params.annotate_vcfanno && params.annotate_vcfanno_toml
    ch_current_vcf_tbi = ch_vcf

    if (run_vcfanno) {
        def resource_paths = (params.annotate_vcfanno_resources ?: '')
            .tokenize(',')
            .collect { path_str -> path_str.trim() }
            .findAll { path_str -> path_str }
            .collect { path_str -> file(path_str, checkIfExists: true) }

        def toml = file(params.annotate_vcfanno_toml, checkIfExists: true)
        def lua = params.annotate_vcfanno_lua
            ? file(params.annotate_vcfanno_lua, checkIfExists: true)
            : []

        VCFANNO(
            ch_current_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf, tbi, []] },
            toml,
            lua,
            resource_paths
        )
        ch_current_vcf_tbi = VCFANNO.out.vcf.join(VCFANNO.out.tbi, failOnDuplicate: true, failOnMismatch: true)
    }

    ch_current = ch_current_vcf_tbi.map { meta, vcf, _tbi -> [meta, vcf] }

    if (tools) {
        def genome = params.annotate_vep_genome ?: (params.genome ?: 'GRCh38')
        def species = params.annotate_vep_species ?: 'homo_sapiens'
        def cache_version = params.annotate_vep_cache_version ?: '110'
        def snpeff_db = params.annotate_snpeff_db ?: 'GRCh38.99'
        def sites_per_chunk = params.annotate_sites_per_chunk

        def ch_effect_in = ch_current_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf, tbi, []] }
        // Value-like channels so scatter shards can reuse fasta/caches
        def ch_fasta_c = ch_fasta.collect()
        def ch_vep_cache_c = ch_vep_cache.collect()
        def ch_snpeff_cache_c = ch_snpeff_cache.collect()

        VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(
            ch_effect_in,
            ch_fasta_c,
            genome,
            species,
            cache_version,
            ch_vep_cache_c,
            [],
            snpeff_db,
            ch_snpeff_cache_c,
            tools,
            sites_per_chunk
        )
        ch_current = VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.vcf_tbi.map { meta, vcf, _tbi -> [meta, vcf] }
    }

    if (params.annotate_filter) {
        def do_bcf = params.annotate_filter_bcftools ? true : false
        def do_vep = (params.annotate_filter_vep && params.annotate_filter_vep_file) ? true : false
        if (do_bcf || do_vep) {
            def ch_feat = params.annotate_filter_vep_file
                ? channel.of([[id: 'vep_filter'], file(params.annotate_filter_vep_file, checkIfExists: true)])
                : channel.of([[id: 'vep_filter'], []])
            VCF_FILTER_BCFTOOLS_ENSEMBLVEP(
                ch_current,
                ch_feat,
                do_bcf,
                do_vep
            )
            ch_current = VCF_FILTER_BCFTOOLS_ENSEMBLVEP.out.vcf
        }
    }

    emit:
    vcf      = ch_current
    versions = ch_versions
}
