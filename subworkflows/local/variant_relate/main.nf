//
// Optional Somalier extract + relate (sample relatedness QC).
// Enable with params.run_relate (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// SOMALIER_EXTRACT.out.versions
// SOMALIER_RELATE.out.versions


workflow VARIANT_RELATE {
    take:
    ch_vcf   // channel: [ meta, vcf, tbi ]
    ch_fasta // channel: [ meta, fasta ]
    ch_fai   // channel: [ meta, fai ]
    ch_sites // channel: [ meta, sites_vcf ]
    ch_ped   // channel: [ meta, ped ] (may be empty list path)

    main:
    ch_versions = channel.empty()

    SOMALIER_EXTRACT(
        ch_vcf,
        ch_fasta,
        ch_fai,
        ch_sites
    )

    ch_relate_in = SOMALIER_EXTRACT.out.extract
        .join(ch_ped, remainder: true)
        .map { meta, extract, ped ->
            [meta, extract instanceof List ? extract : [extract], ped ?: []]
        }

    SOMALIER_RELATE(
        ch_relate_in,
        []
    )

    emit:
    extract     = SOMALIER_EXTRACT.out.extract
    html        = SOMALIER_RELATE.out.html
    pairs_tsv   = SOMALIER_RELATE.out.pairs_tsv
    samples_tsv = SOMALIER_RELATE.out.samples_tsv
    versions    = ch_versions
}
