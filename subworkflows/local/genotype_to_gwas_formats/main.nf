//
// Minimal genotype format fan-out for the GWAS benchmark branch.
// Emits bed/bim/fam (passthrough), vcf.gz (passthrough or PLINK_RECODE),
// and tped/tfam via PLINK_RECODE --recode transpose.
//
// Include nextflow.config from this directory (or equivalent process withName rules):
//   includeConfig 'subworkflows/local/genotype_to_gwas_formats/nextflow.config'
//

include { PLINK_RECODE as PLINK_RECODE_TRANSPOSE } from '../../../modules/nf-core/plink/recode/main'
include { PLINK_RECODE as PLINK_RECODE_VCF       } from '../../../modules/nf-core/plink/recode/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// PLINK_RECODE_TRANSPOSE.out.versions, PLINK_RECODE_VCF.out.versions


workflow GENOTYPE_TO_GWAS_FORMATS {
    take:
    ch_plink  // channel: [ meta, bed, bim, fam ]
    ch_vcf    // channel: [ meta, vcf ] (empty → convert from bed)

    main:
    // versions via topic("versions") on PLINK_RECODE — keep Path YAML channel empty
    ch_versions = channel.empty()

    // Passthrough PLINK binary
    ch_bed = ch_plink

    // tped/tfam via --recode transpose (ext.args in nextflow.config)
    PLINK_RECODE_TRANSPOSE(ch_plink)
    ch_tped = PLINK_RECODE_TRANSPOSE.out.tped
        .join(PLINK_RECODE_TRANSPOSE.out.tfam)

    // VCF: passthrough if provided for meta.id, else convert bed → vcf.gz
    ch_joined = ch_plink
        .join(ch_vcf, remainder: true)
        .map { meta, bed, bim, fam, vcf ->
            [meta, bed, bim, fam, vcf]
        }
        .branch { _meta, _bed, _bim, _fam, vcf ->
            has_vcf: vcf != null
            need_vcf: true
        }

    ch_vcf_passthrough = ch_joined.has_vcf
        .map { meta, _bed, _bim, _fam, vcf -> [meta, vcf] }

    ch_bed_for_vcf = ch_joined.need_vcf
        .map { meta, bed, bim, fam, _vcf -> [meta, bed, bim, fam] }

    PLINK_RECODE_VCF(ch_bed_for_vcf)

    ch_vcf_converted = PLINK_RECODE_VCF.out.vcfgz
        .mix(PLINK_RECODE_VCF.out.vcf)

    ch_vcf_out = ch_vcf_passthrough.mix(ch_vcf_converted)

    emit:
    bed      = ch_bed       // channel: [ meta, bed, bim, fam ]
    vcf      = ch_vcf_out   // channel: [ meta, vcf.gz|vcf ]
    tped     = ch_tped      // channel: [ meta, tped, tfam ]
    versions = ch_versions  // channel: versions topic / tuples
}
