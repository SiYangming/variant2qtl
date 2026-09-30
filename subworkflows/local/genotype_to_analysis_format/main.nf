//
// Convert genotypes to analysis-friendly formats for OmiGA / tensorQTL / etc.
// Always emits PLINK bed passthrough; optional VCF (from bed via PLINK_VCF reverse
// is not here — pass VCF in) and optional BGEN via PLINK2_VCF2BGEN when VCF present.
// Include: includeConfig 'subworkflows/local/genotype_to_analysis_format/nextflow.config'
//

include { PLINK2_VCF2BGEN } from '../../../modules/nf-core/plink2/vcf2bgen/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// PLINK2_VCF2BGEN.out.versions

workflow GENOTYPE_TO_ANALYSIS_FORMAT {
    take:
    ch_plink  // channel: [ meta, bed, bim, fam ]
    ch_vcf    // channel: [ meta, vcf ] (may be empty)

    main:
    ch_versions = channel.empty()

    ch_bed = ch_plink
    ch_bgen = channel.empty()
    ch_sample = channel.empty()

    if (params.genotype_to_bgen) {
        // PLINK2_VCF2BGEN needs dosage field + flags; use defaults when VCF available
        ch_vcf_in = ch_vcf.map { meta, vcf ->
            [meta, vcf, params.genotype_bgen_dosage ?: 'DS', true, 'double-id']
        }

        PLINK2_VCF2BGEN(ch_vcf_in)
        ch_versions = ch_versions.mix(PLINK2_VCF2BGEN.out.versions)
        ch_bgen = PLINK2_VCF2BGEN.out.bgen_file
        ch_sample = PLINK2_VCF2BGEN.out.sample_file
    }

    emit:
    bed      = ch_bed      // channel: [ meta, bed, bim, fam ]
    vcf      = ch_vcf      // channel: [ meta, vcf ]
    bgen     = ch_bgen     // channel: [ meta, bgen ] (optional)
    sample   = ch_sample   // channel: [ meta, sample ] (optional)
    versions = ch_versions
}
