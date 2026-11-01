include { PREPARE_MAG_PROTEIN_DATABASE } from '../../../modules/local/prepare_mag_protein_database/main.nf'
include { DIAMOND } from '../../../modules/core/diamond/main.nf'
include { FILTER_DIAMOND_BEST_HITS } from '../../../modules/local/filter_diamond_best_hits/main.nf'
include { CREATE_TRANSCRIPT_MAG_GTF } from '../../../modules/local/create_transcript_mag_gtf/main.nf'

workflow MAG_ANNOTATION {

    take:
    pep
    gff3
    mag_protein_dir_ch

    main:

    PREPARE_MAG_PROTEIN_DATABASE(
        mag_protein_dir_ch,
        "mag_annotation"
    )

    DIAMOND(
        pep,
        PREPARE_MAG_PROTEIN_DATABASE.out.db,
        "coassembly",
        "mag_annotation"
    )

    FILTER_DIAMOND_BEST_HITS(
        DIAMOND.out.hits,
        "coassembly",
        "mag_annotation"
    )

    CREATE_TRANSCRIPT_MAG_GTF(
        FILTER_DIAMOND_BEST_HITS.out.best_hits,
        gff3,
        "coassembly",
        "mag_annotation"
    )

    emit:
    transcript_mag_map = CREATE_TRANSCRIPT_MAG_GTF.out.map
    transcript_mag_gtf = CREATE_TRANSCRIPT_MAG_GTF.out.gtf
    diamond_hits       = DIAMOND.out.hits
    diamond_best_hits  = FILTER_DIAMOND_BEST_HITS.out.best_hits

    prepare_versions = PREPARE_MAG_PROTEIN_DATABASE.out.versions
    diamond_log      = DIAMOND.out.log
    diamond_versions = DIAMOND.out.versions
    filter_versions  = FILTER_DIAMOND_BEST_HITS.out.versions
    gtf_log          = CREATE_TRANSCRIPT_MAG_GTF.out.log
    gtf_versions     = CREATE_TRANSCRIPT_MAG_GTF.out.versions
}
