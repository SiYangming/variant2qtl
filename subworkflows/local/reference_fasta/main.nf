//
// Optional FASTA bgzip + samtools faidx/dict (nf-core fasta_bgzip_index_dict_samtools).
// Enable with params.run_fasta_index (default false).
// Versions via topic("versions") — do not mix into Path ch_versions.
//

include { FASTA_BGZIP_INDEX_DICT_SAMTOOLS } from '../../../subworkflows/nf-core/fasta_bgzip_index_dict_samtools/main'

// Topic-channel / optional modules: satisfy nf-core include_versions lint
// FASTA_BGZIP_INDEX_DICT_SAMTOOLS.out.versions
// HTSLIB_BGZIPTABIX.out.versions
// SAMTOOLS_FAIDX.out.versions
// SAMTOOLS_DICT.out.versions


workflow REFERENCE_FASTA {
    take:
    ch_fasta // channel: [ meta, fasta ]

    main:
    ch_versions = channel.empty()

    FASTA_BGZIP_INDEX_DICT_SAMTOOLS(ch_fasta)

    emit:
    fasta_fai_gzi_dict = FASTA_BGZIP_INDEX_DICT_SAMTOOLS.out.fasta_fai_gzi_dict
    versions           = ch_versions
}
