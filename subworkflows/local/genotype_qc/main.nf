//
// Minimal genotype QC: MAF / HWE / missingness via PLINK2_FILTER → filtered bed/bim/fam.
// Include: includeConfig 'subworkflows/local/genotype_qc/nextflow.config'
//

include { PLINK2_FILTER } from '../../../modules/nf-core/plink2/filter/main'

workflow GENOTYPE_QC {
    take:
    ch_plink  // channel: [ meta, bed, bim, fam ]

    main:
    ch_versions = Channel.empty()

    // PLINK2_FILTER expects [ meta, genotype, variant, sample ] with matching basenames
    ch_plink_in = ch_plink.map { meta, bed, bim, fam -> [meta, bed, bim, fam] }

    PLINK2_FILTER(ch_plink_in)

    ch_bed_out = PLINK2_FILTER.out.bed
        .join(PLINK2_FILTER.out.bim)
        .join(PLINK2_FILTER.out.fam)

    // Path-based versions.yml from plink2 filter
    ch_versions = ch_versions.mix(PLINK2_FILTER.out.versions)

    emit:
    bed      = ch_bed_out   // channel: [ meta, bed, bim, fam ]
    versions = ch_versions  // channel: path(versions.yml)
}
