//
// Minimal OmiGA cis-molQTL mapping stub.
// Enable from the main workflow with params.run_omiga_cis (default false).
// Versions are emitted via topic("versions") on OMIGA_CIS — do not mix into Path ch_versions.
//

include { OMIGA_CIS } from '../../../modules/local/omiga/cis/main'

workflow MOLQTL_MAP_OMIGA {
    take:
    ch_plink       // channel: [ meta, bed, bim, fam ]
    ch_phenotype   // channel: [ meta, phenotype BED/OPF ]
    ch_covariates  // channel: [ meta, covar ] (may be empty path list)

    main:
    // Local omiga modules publish versions via topic("versions").
    // Do not mix those tuples into Path-based ch_versions (breaks softwareVersionsToYAML).
    ch_versions = channel.empty()

    ch_cov_aligned = ch_phenotype
        .join(ch_covariates, remainder: true)
        .map { meta, _pheno, cov -> [meta, cov ?: []] }

    OMIGA_CIS(
        ch_plink,
        ch_phenotype,
        ch_cov_aligned
    )

    emit:
    cis_qtl  = OMIGA_CIS.out.cis_qtl
    outdir   = OMIGA_CIS.out.outdir
    versions = ch_versions
}
