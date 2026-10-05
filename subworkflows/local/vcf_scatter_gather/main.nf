//
// Scatter a VCF over BED intervals, or gather scattered VCFs by meta.group_id.
// Versions via Path from BED_TO_REGION — do not mix topic versions into ch_versions.
//

include { BED_SCATTER_BEDTOOLS } from '../../../subworkflows/nf-core/bed_scatter_bedtools/main'
include { VCF_GATHER_BCFTOOLS  } from '../../../subworkflows/nf-core/vcf_gather_bcftools/main'
include { BED_TO_REGION        } from '../../../modules/local/utils/bed_to_region/main'
include { BCFTOOLS_INDEX       } from '../../../modules/nf-core/bcftools/index/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// BED_TO_REGION.out.versions, BCFTOOLS_INDEX.out.versions
// BED_SCATTER_BEDTOOLS.out.versions, VCF_GATHER_BCFTOOLS.out.versions


workflow VCF_SCATTER_REGIONS {
    take:
    ch_vcf // channel: [ meta, vcf, tbi ]
    ch_bed // channel: [ meta, bed, scatter_count ]

    main:
    ch_versions = channel.empty()

    BED_SCATTER_BEDTOOLS(ch_bed)
    BED_TO_REGION(BED_SCATTER_BEDTOOLS.out.scattered_beds)
    ch_versions = ch_versions.mix(BED_TO_REGION.out.versions)

    ch_scattered = ch_vcf
        .combine(BED_TO_REGION.out.region)
        .map { meta, vcf, tbi, _bmeta, region_file, scatter_count ->
            def region = region_file.text.trim()
            def chunk = region_file.name.replaceFirst(/\.region\.txt$/, '')
            [
                meta + [group_id: meta.id, id: "${meta.id}_${chunk}", region: region, scatter_count: scatter_count],
                vcf,
                tbi
            ]
        }

    emit:
    vcf      = ch_scattered
    versions = ch_versions
}

workflow VCF_GATHER_REGIONS {
    take:
    ch_vcf // channel: [ meta, vcf ]  (meta.group_id + meta.scatter_count required)

    main:
    BCFTOOLS_INDEX(ch_vcf)
    VCF_GATHER_BCFTOOLS(
        ch_vcf
            .join(BCFTOOLS_INDEX.out.index)
            .map { meta, vcf, index -> [meta, vcf, index, meta.scatter_count] },
        ['group_id'],
        false
    )

    ch_out = VCF_GATHER_BCFTOOLS.out.vcf_index.map { meta, vcf, tbi ->
        def v = vcf instanceof List ? vcf[0] : vcf
        def t = tbi instanceof List ? (tbi ? tbi[0] : []) : (tbi ?: [])
        [meta + [id: (meta.group_id ?: meta.id)], v, t]
    }

    emit:
    vcf      = ch_out
    versions = channel.empty()
}
