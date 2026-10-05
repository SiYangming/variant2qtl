//
// Optional germline SV calling from BAM (smoove / manta / delly).
// Enable with params.run_sv (default false).
// Versions via topic("versions") on nf-core callers — do not mix into Path ch_versions.
//

include { SMOOVE_CALL    } from '../../../modules/nf-core/smoove/call/main'
include { MANTA_GERMLINE } from '../../../modules/nf-core/manta/germline/main'
include { DELLY_CALL     } from '../../../modules/nf-core/delly/call/main'
include { GUNZIP as GUNZIP_SV } from '../../../modules/nf-core/gunzip/main'
include { SURVIVOR_MERGE } from '../../../modules/nf-core/survivor/merge/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// SMOOVE_CALL.out.versions
// MANTA_GERMLINE.out.versions
// DELLY_CALL.out.versions
// GUNZIP_SV.out.versions
// SURVIVOR_MERGE.out.versions


workflow VARIANT_SV {
    take:
    ch_bam    // channel: [ meta, bam, bai ]
    ch_fasta  // channel: [ meta, fasta ]
    ch_fai    // channel: [ meta, fai ]

    main:
    ch_versions = channel.empty()
    ch_vcf = channel.empty()

    def engines = (params.sv_engines ?: 'smoove')
        .tokenize(',')
        .collect { engine -> engine.trim().toLowerCase() }

    if (engines.contains('smoove')) {
        SMOOVE_CALL(
            ch_bam.map { meta, bam, bai -> [meta, bam, bai, []] },
            ch_fasta,
            ch_fai
        )
        ch_vcf = ch_vcf.mix(SMOOVE_CALL.out.vcf)
    }

    if (engines.contains('manta')) {
        MANTA_GERMLINE(
            ch_bam.map { meta, bam, bai -> [meta, bam, bai, [], []] },
            ch_fasta,
            ch_fai,
            []
        )
        ch_vcf = ch_vcf.mix(MANTA_GERMLINE.out.diploid_sv_vcf)
    }

    if (engines.contains('delly')) {
        DELLY_CALL(
            ch_bam.map { meta, bam, bai -> [meta, bam, bai, [], [], []] },
            ch_fasta,
            ch_fai,
            'vcf'
        )
        ch_vcf = ch_vcf.mix(DELLY_CALL.out.bcf)
    }

    // Merge multi-engine calls (SURVIVOR). Gzipped VCFs are gunzipped first.
    if (params.sv_merge != false) {
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
        ch_by_ext = ch_flat.branch { _meta, vcf ->
            gz: vcf.toString().endsWith('.gz')
            rest: true
        }
        GUNZIP_SV(ch_by_ext.gz)
        ch_uncomp = ch_by_ext.rest.mix(GUNZIP_SV.out.gunzip)
        SURVIVOR_MERGE(
            ch_uncomp
                .map { meta, vcf -> [meta.group_id ?: meta.id, meta, vcf] }
                .groupTuple()
                .map { _id, metas, vcfs -> [metas[0], vcfs] },
            params.sv_merge_max_dist ?: 1000,
            params.sv_merge_min_callers ?: 1,
            1,
            1,
            0,
            params.sv_merge_min_size ?: 0
        )
        ch_vcf = ch_one.mix(
            SURVIVOR_MERGE.out.vcf.map { meta, vcf ->
                [meta + [id: meta.group_id ?: meta.id], vcf]
            }
        )
    }

    emit:
    vcf      = ch_vcf
    versions = ch_versions
}
