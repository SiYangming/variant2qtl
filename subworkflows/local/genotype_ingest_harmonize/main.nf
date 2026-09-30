//
// Minimal genotype ingest: VCF → PLINK bed/bim/fam via PLINK_VCF.
// No liftover / allele-flip / SV-STR yet (full genotype_ingest_harmonize is P0 later).
// Include: includeConfig 'subworkflows/local/genotype_ingest_harmonize/nextflow.config'
//

include { PLINK_VCF } from '../../../modules/nf-core/plink/vcf/main'

workflow GENOTYPE_INGEST_HARMONIZE {
    take:
    ch_vcf  // channel: [ meta, vcf ] or [ meta, vcf, tbi ]

    main:
    // PLINK_VCF uses topic("versions") — keep Path YAML channel empty
    ch_versions = channel.empty()

    // Drop optional tbi; PLINK_VCF only needs [ meta, vcf ]
    ch_vcf_in = ch_vcf.map { row ->
        def meta = row[0]
        def vcf  = row[1]
        [meta, vcf]
    }

    PLINK_VCF(ch_vcf_in)

    ch_bed_out = PLINK_VCF.out.bed
        .join(PLINK_VCF.out.bim)
        .join(PLINK_VCF.out.fam)

    // Passthrough input VCF (same meta) for downstream engines that prefer VCF
    ch_vcf_out = ch_vcf_in

    emit:
    bed      = ch_bed_out  // channel: [ meta, bed, bim, fam ]
    vcf      = ch_vcf_out  // channel: [ meta, vcf ]
    versions = ch_versions // channel: Path versions.yml (empty; plink via topic)
}
