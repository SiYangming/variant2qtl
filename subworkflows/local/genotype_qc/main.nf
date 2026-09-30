//
// Genotype QC: MAF/HWE/geno (PLINK2_FILTER), optional het outliers, relatedness, PCA.
// Enable extras via params.genotype_qc_run_{het,relatedness,pca} (default false).
// Include: includeConfig 'subworkflows/local/genotype_qc/nextflow.config'
//

include { PLINK2_FILTER                         } from '../../../modules/nf-core/plink2/filter/main'
include { PLINK2_HET                            } from '../../../modules/nf-core/plink2/het/main'
include { PLINK2_REMOVE as PLINK2_REMOVE_HET    } from '../../../modules/nf-core/plink2/remove/main'
include { PLINK2_REMOVE as PLINK2_REMOVE_REL    } from '../../../modules/nf-core/plink2/remove/main'
include { PLINK_GENOME                          } from '../../../modules/nf-core/plink/genome/main'
include { HET_OUTLIERS                          } from '../../../modules/local/utils/het_outliers/main'
include { RELATEDNESS_OUTLIERS                  } from '../../../modules/local/utils/relatedness_outliers/main'
include { PLINK2_PCA_BFILE                      } from '../../../modules/local/utils/plink2_pca_bfile/main'

workflow GENOTYPE_QC {
    take:
    ch_plink  // channel: [ meta, bed, bim, fam ]

    main:
    ch_versions = channel.empty()

    ch_plink_in = ch_plink.map { meta, bed, bim, fam -> [meta, bed, bim, fam] }

    PLINK2_FILTER(ch_plink_in)

    ch_bed = PLINK2_FILTER.out.bed
        .join(PLINK2_FILTER.out.bim)
        .join(PLINK2_FILTER.out.fam)

    ch_versions = ch_versions.mix(PLINK2_FILTER.out.versions)

    ch_het_report    = channel.empty()
    ch_genome_report = channel.empty()
    ch_eigenvec      = channel.empty()
    ch_eigenval      = channel.empty()

    // --- Optional heterozygosity outlier removal ---
    if (params.genotype_qc_run_het) {
        PLINK2_HET(ch_bed)
        ch_versions = ch_versions.mix(PLINK2_HET.out.versions)
        ch_het_report = PLINK2_HET.out.het

        HET_OUTLIERS(PLINK2_HET.out.het)
        ch_versions = ch_versions.mix(HET_OUTLIERS.out.versions)

        ch_remove = HET_OUTLIERS.out.outliers
            .filter { _meta, path -> path.size() > 0 }

        ch_keep = ch_bed
            .join(HET_OUTLIERS.out.outliers)
            .filter { _meta, _bed, _bim, _fam, path -> path.size() == 0 }
            .map { meta, bed, bim, fam, _path -> [meta, bed, bim, fam] }

        ch_to_remove = ch_bed.join(ch_remove)

        PLINK2_REMOVE_HET(
            ch_to_remove.map { meta, bed, bim, fam, _outliers -> [meta, bed, bim, fam] },
            ch_to_remove.map { _meta, _bed, _bim, _fam, outliers -> outliers }
        )
        ch_versions = ch_versions.mix(PLINK2_REMOVE_HET.out.versions)

        ch_removed_bed = PLINK2_REMOVE_HET.out.remove_bed
            .join(PLINK2_REMOVE_HET.out.remove_bim)
            .join(PLINK2_REMOVE_HET.out.remove_fam)

        ch_bed = ch_keep.mix(ch_removed_bed)
    }

    // --- Optional relatedness (IBD) report + remove second of each pair ---
    if (params.genotype_qc_run_relatedness) {
        PLINK_GENOME(ch_bed)
        ch_genome_report = PLINK_GENOME.out.genome

        RELATEDNESS_OUTLIERS(PLINK_GENOME.out.genome)
        ch_versions = ch_versions.mix(RELATEDNESS_OUTLIERS.out.versions)

        ch_rel_remove = RELATEDNESS_OUTLIERS.out.outliers
            .filter { _meta, path -> path.size() > 0 }

        ch_rel_keep = ch_bed
            .join(RELATEDNESS_OUTLIERS.out.outliers)
            .filter { _meta, _bed, _bim, _fam, path -> path.size() == 0 }
            .map { meta, bed, bim, fam, _path -> [meta, bed, bim, fam] }

        ch_rel_to_remove = ch_bed.join(ch_rel_remove)

        PLINK2_REMOVE_REL(
            ch_rel_to_remove.map { meta, bed, bim, fam, _outliers -> [meta, bed, bim, fam] },
            ch_rel_to_remove.map { _meta, _bed, _bim, _fam, outliers -> outliers }
        )
        ch_versions = ch_versions.mix(PLINK2_REMOVE_REL.out.versions)

        ch_rel_removed = PLINK2_REMOVE_REL.out.remove_bed
            .join(PLINK2_REMOVE_REL.out.remove_bim)
            .join(PLINK2_REMOVE_REL.out.remove_fam)

        ch_bed = ch_rel_keep.mix(ch_rel_removed)
    }

    // --- Optional PCA on filtered bed ---
    if (params.genotype_qc_run_pca) {
        PLINK2_PCA_BFILE(ch_bed)
        ch_versions = ch_versions.mix(PLINK2_PCA_BFILE.out.versions)
        ch_eigenvec = PLINK2_PCA_BFILE.out.eigenvec
        ch_eigenval = PLINK2_PCA_BFILE.out.eigenval
    }

    emit:
    bed      = ch_bed
    het      = ch_het_report
    genome   = ch_genome_report
    eigenvec = ch_eigenvec
    eigenval = ch_eigenval
    versions = ch_versions
}
