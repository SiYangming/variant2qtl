//
// Optional germline SV calling from BAM (smoove / manta / delly).
// Enable with params.run_sv (default false).
// Versions via topic("versions") on nf-core callers — do not mix into Path ch_versions.
//

include { SMOOVE_CALL    } from '../../../modules/nf-core/smoove/call/main'
include { MANTA_GERMLINE } from '../../../modules/nf-core/manta/germline/main'
include { DELLY_CALL     } from '../../../modules/nf-core/delly/call/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// SMOOVE_CALL.out.versions
// MANTA_GERMLINE.out.versions
// DELLY_CALL.out.versions


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

    emit:
    vcf      = ch_vcf
    versions = ch_versions
}
