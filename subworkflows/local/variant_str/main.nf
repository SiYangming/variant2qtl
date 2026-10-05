//
// Optional STR genotyping from BAM (ExpansionHunter / GangSTR / HipSTR / TRGT).
// Enable with params.run_str (default false).
// Versions via topic("versions") on nf-core callers — do not mix into Path ch_versions.
//

include { EXPANSIONHUNTER } from '../../../modules/nf-core/expansionhunter/main'
include { GANGSTR         } from '../../../modules/nf-core/gangstr/main'
include { HIPSTR          } from '../../../modules/nf-core/hipstr/main'
include { TRGT_GENOTYPE   } from '../../../modules/nf-core/trgt/genotype/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_STR } from '../../../modules/nf-core/bcftools/index/main'
include { TRTOOLS_MERGESTR } from '../../../modules/nf-core/trtools/mergestr/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// EXPANSIONHUNTER.out.versions
// GANGSTR.out.versions
// HIPSTR.out.versions
// TRGT_GENOTYPE.out.versions
// BCFTOOLS_INDEX_STR.out.versions
// TRTOOLS_MERGESTR.out.versions


workflow VARIANT_STR {
    take:
    ch_bam      // channel: [ meta, bam, bai ]
    ch_fasta    // channel: [ meta, fasta ]
    ch_fai      // channel: [ meta, fai ]
    ch_catalog  // channel: [ meta, catalog ]
    ch_regions  // channel: [ meta, regions ]
    ch_repeats  // channel: [ meta, repeats ]
    ch_hipstr   // channel: [ meta, hipstr_bed ]

    main:
    ch_versions = channel.empty()
    ch_vcf = channel.empty()

    def engines = (params.str_engines ?: 'expansionhunter')
        .tokenize(',')
        .collect { engine -> engine.trim().toLowerCase() }

    if (engines.contains('expansionhunter')) {
        EXPANSIONHUNTER(
            ch_bam,
            ch_fasta,
            ch_fai,
            ch_catalog
        )
        ch_vcf = ch_vcf.mix(EXPANSIONHUNTER.out.vcf)
    }

    if (engines.contains('gangstr')) {
        GANGSTR(
            ch_bam
                .combine(ch_regions)
                .map { meta, bam, bai, _reg_meta, regions ->
                    [meta, bam, bai, regions]
                },
            ch_fasta.map { _meta, fasta -> fasta },
            ch_fai.map { _meta, fai -> fai }
        )
        ch_vcf = ch_vcf.mix(GANGSTR.out.vcf)
    }

    if (engines.contains('hipstr')) {
        HIPSTR(
            ch_bam.map { meta, bam, bai -> [meta, bam, bai, []] },
            ch_fasta
                .combine(ch_fai)
                .map { meta, fasta, _fai_meta, fai -> [meta, fasta, fai] },
            ch_hipstr,
            channel.of([[], [], []]),
            channel.of([[], [], []]),
            [],
            [],
            [],
            false,
            false,
            false
        )
        ch_vcf = ch_vcf.mix(HIPSTR.out.vcf)
    }

    if (engines.contains('trgt')) {
        TRGT_GENOTYPE(
            ch_bam.map { meta, bam, bai -> [meta, bam, bai, []] },
            ch_fasta,
            ch_fai,
            ch_repeats
        )
        ch_vcf = ch_vcf.mix(TRGT_GENOTYPE.out.vcf)
    }

    if (params.str_merge != false) {
        ch_grouped = ch_vcf
            .map { meta, vcf -> [meta.id, meta, vcf] }
            .groupTuple()
            .map { _id, metas, vcfs ->
                def files = vcfs instanceof Collection ? vcfs.flatten() : [vcfs]
                [metas[0], files]
            }
            .branch { _meta, files ->
                one: files.size() <= 1
                more: true
            }
        ch_one = ch_grouped.one.map { meta, files -> [meta, files[0]] }
        ch_flat = ch_grouped.more.transpose().map { meta, vcf ->
            def stem = vcf.name.replaceFirst(/\.vcf(\.gz)?$/, '')
            [meta + [id: "${meta.id}_${stem}", group_id: meta.id], vcf]
        }
        BCFTOOLS_INDEX_STR(ch_flat)
        TRTOOLS_MERGESTR(
            ch_flat
                .join(BCFTOOLS_INDEX_STR.out.index)
                .map { meta, vcf, tbi -> [meta.group_id ?: meta.id, meta, vcf, tbi] }
                .groupTuple()
                .map { _id, metas, vcfs, tbis -> [metas[0], vcfs, tbis] }
        )
        ch_vcf = ch_one.mix(
            TRTOOLS_MERGESTR.out.vcf.map { meta, vcf ->
                [meta + [id: meta.group_id ?: meta.id], vcf]
            }
        )
    }

    emit:
    vcf      = ch_vcf
    versions = ch_versions
}
