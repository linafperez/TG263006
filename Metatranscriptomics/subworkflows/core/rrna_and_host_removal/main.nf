include { SORTMERNA } from '../../../modules/core/sortmerna/main.nf'
include { STAR } from '../../../modules/core/star/main.nf'

workflow RRNA_AND_HOST_REMOVAL {

    take:
    trimmed_reads_ch
    rrna_db_dir_ch 
    rrna_idx_dir_ch
    star_index_ch

    main:

    SORTMERNA(
        trimmed_reads_ch,
        rrna_db_dir_ch,
        rrna_idx_dir_ch,
        "sortmerna_non_rrna"
    )

    non_rrna_reads_ch = SORTMERNA.out.reads
        .map { meta, read1, read2 ->
            tuple(
                meta,
                [
                    read1,
                    read2
                ]
            )
        }

    STAR(
        non_rrna_reads_ch,
        star_index_ch,
        "star_non_human"
    )

    emit:
    non_human_reads = STAR.out.reads

    sortmerna_log = SORTMERNA.out.log
    star_log      = STAR.out.log

    sortmerna_versions = SORTMERNA.out.versions
    star_versions      = STAR.out.versions
}
