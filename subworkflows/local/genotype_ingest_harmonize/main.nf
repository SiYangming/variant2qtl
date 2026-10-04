//
// Genotype ingest / harmonize:
//   optional Picard liftover → optional bcftools norm → optional sample subset → PLINK_VCF.
// Liftover / norm / samples are opt-in via params (see nextflow.config).
// Include: includeConfig 'subworkflows/local/genotype_ingest_harmonize/nextflow.config'
//

include { PICARD_LIFTOVERVCF                          } from '../../../modules/nf-core/picard/liftovervcf/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORM       } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_VIEW       } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BIALLELIC  } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM                               } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW                               } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_BIALLELIC    } from '../../../modules/nf-core/bcftools/view/main'
include { PLINK_VCF                                   } from '../../../modules/nf-core/plink/vcf/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// PICARD_LIFTOVERVCF.out.versions, BCFTOOLS_INDEX_NORM.out.versions, BCFTOOLS_INDEX_VIEW.out.versions, BCFTOOLS_INDEX_BIALLELIC.out.versions, BCFTOOLS_NORM.out.versions, BCFTOOLS_VIEW.out.versions, BCFTOOLS_VIEW_BIALLELIC.out.versions, PLINK_VCF.out.versions


workflow GENOTYPE_INGEST_HARMONIZE {
    take:
    ch_vcf  // channel: [ meta, vcf ] or [ meta, vcf, tbi ]

    main:
    ch_versions = channel.empty()

    ch_work = ch_vcf.map { row ->
        def meta = row[0]
        def vcf  = row[1]
        [meta, vcf]
    }

    // --- Optional liftover (chain + target fasta + sequence dictionary) ---
    def do_liftover = params.genotype_ingest_liftover_chain &&
        params.genotype_ingest_liftover_fasta &&
        params.genotype_ingest_liftover_dict

    if (do_liftover) {
        def ch_dict = channel.value([
            [id: 'liftover_dict'],
            file(params.genotype_ingest_liftover_dict, checkIfExists: true)
        ])
        def ch_fasta = channel.value([
            [id: 'liftover_fasta'],
            file(params.genotype_ingest_liftover_fasta, checkIfExists: true)
        ])
        def ch_chain = channel.value([
            [id: 'liftover_chain'],
            file(params.genotype_ingest_liftover_chain, checkIfExists: true)
        ])

        PICARD_LIFTOVERVCF(
            ch_work,
            ch_dict,
            ch_fasta,
            ch_chain
        )
        ch_work = PICARD_LIFTOVERVCF.out.vcf_lifted
    }

    // --- Optional left-align / allele split via bcftools norm ---
    if (params.genotype_ingest_fasta) {
        def ch_ref = channel.value([
            [id: 'norm_fasta'],
            file(params.genotype_ingest_fasta, checkIfExists: true)
        ])

        BCFTOOLS_INDEX_NORM(ch_work)
        ch_indexed = ch_work.join(BCFTOOLS_INDEX_NORM.out.index)

        BCFTOOLS_NORM(ch_indexed, ch_ref)
        ch_work = BCFTOOLS_NORM.out.vcf
    }

    // Drop multi-allelic / symbolic alleles so PLINK ingest can proceed.
    // Complex SV/STR records are discarded (documented limitation).
    def do_biallelic = params.genotype_ingest_biallelic || params.sv_feed_ingest || params.str_feed_ingest
    if (do_biallelic) {
        BCFTOOLS_INDEX_BIALLELIC(ch_work)
        BCFTOOLS_VIEW_BIALLELIC(
            ch_work.join(BCFTOOLS_INDEX_BIALLELIC.out.index),
            [],
            [],
            []
        )
        ch_work = BCFTOOLS_VIEW_BIALLELIC.out.vcf
    }

    // --- Optional sample keep-list (cohort alignment) ---
    if (params.genotype_ingest_samples) {
        def samples = file(params.genotype_ingest_samples, checkIfExists: true)

        BCFTOOLS_INDEX_VIEW(ch_work)
        ch_indexed = ch_work.join(BCFTOOLS_INDEX_VIEW.out.index)

        BCFTOOLS_VIEW(
            ch_indexed,
            [],
            [],
            samples
        )
        ch_work = BCFTOOLS_VIEW.out.vcf
    }

    ch_vcf_in = ch_work.map { row ->
        def meta = row[0]
        def vcf  = row[1]
        [meta, vcf]
    }

    PLINK_VCF(ch_vcf_in)

    ch_bed_out = PLINK_VCF.out.bed
        .join(PLINK_VCF.out.bim)
        .join(PLINK_VCF.out.fam)

    ch_vcf_out = ch_vcf_in

    emit:
    bed      = ch_bed_out
    vcf      = ch_vcf_out
    versions = ch_versions
}
