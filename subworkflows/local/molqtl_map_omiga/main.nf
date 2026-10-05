//
// OmiGA molQTL mapping (cis, trans, independent-cis).
// Optionally fans genotypes through genotype_to_analysis_format (BGEN) first.
// Enable from the main workflow with params.run_omiga_cis / run_omiga_trans /
// run_omiga_independent_cis (all default false). Do not mix topic versions
// into Path ch_versions.
//

include { OMIGA_CIS                          } from '../../../modules/local/omiga/cis/main'
include { OMIGA_CIS as OMIGA_TRANS           } from '../../../modules/local/omiga/cis/main'
include { OMIGA_CIS as OMIGA_INDEPENDENT     } from '../../../modules/local/omiga/cis/main'
include { GENOTYPE_TO_ANALYSIS_FORMAT        } from '../genotype_to_analysis_format/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// OMIGA_CIS.out.versions, OMIGA_TRANS.out.versions, OMIGA_INDEPENDENT.out.versions
// GENOTYPE_TO_ANALYSIS_FORMAT.out.versions


workflow MOLQTL_MAP_OMIGA {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_phenotype   // channel: [ meta, phenotype BED/OPF ]
    ch_covariates  // channel: [ meta, covar ] (may be empty path list)
    ch_vcf         // channel: [ meta, vcf ] (optional; empty OK)

    main:
    ch_versions = channel.empty()
    ch_cis_qtl = channel.empty()
    ch_outdir = channel.empty()

    GENOTYPE_TO_ANALYSIS_FORMAT(ch_plink, ch_vcf)
    ch_versions = ch_versions.mix(GENOTYPE_TO_ANALYSIS_FORMAT.out.versions)

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    def do_trans = params.run_omiga_trans
    def do_indep = params.run_omiga_independent_cis
    def do_cis = params.run_omiga_cis || (!do_trans && !do_indep)

    if (do_cis) {
        OMIGA_CIS(
            GENOTYPE_TO_ANALYSIS_FORMAT.out.bed,
            ch_phenotype,
            ch_cov_aligned
        )
        ch_cis_qtl = ch_cis_qtl.mix(OMIGA_CIS.out.cis_qtl)
        ch_outdir = ch_outdir.mix(OMIGA_CIS.out.outdir)
    }

    if (do_trans) {
        OMIGA_TRANS(
            GENOTYPE_TO_ANALYSIS_FORMAT.out.bed,
            ch_phenotype,
            ch_cov_aligned
        )
        ch_cis_qtl = ch_cis_qtl.mix(OMIGA_TRANS.out.cis_qtl)
        ch_outdir = ch_outdir.mix(OMIGA_TRANS.out.outdir)
    }

    if (do_indep) {
        OMIGA_INDEPENDENT(
            GENOTYPE_TO_ANALYSIS_FORMAT.out.bed,
            ch_phenotype,
            ch_cov_aligned
        )
        ch_cis_qtl = ch_cis_qtl.mix(OMIGA_INDEPENDENT.out.cis_qtl)
        ch_outdir = ch_outdir.mix(OMIGA_INDEPENDENT.out.outdir)
    }

    emit:
    cis_qtl  = ch_cis_qtl
    outdir   = ch_outdir
    bed      = GENOTYPE_TO_ANALYSIS_FORMAT.out.bed
    bgen     = GENOTYPE_TO_ANALYSIS_FORMAT.out.bgen
    versions = ch_versions
}
