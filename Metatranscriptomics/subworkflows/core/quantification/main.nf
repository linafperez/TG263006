include { SALMON_INDEX; SALMON_QUANT } from '../../../modules/core/salmon/main.nf'

workflow QUANTIFICATION {

    take:
    transcripts
    non_human_reads_ch

    main:

    SALMON_INDEX(
        transcripts,
        "coassembly",
        "coassembly"
    )

    /*
     * The same transcriptome index is reused for every sample.
     */
    salmon_index_ch = SALMON_INDEX.out.index

    salmon_reads_ch = non_human_reads_ch
        .map { meta, read1, read2 ->
            tuple(
                meta,
                [
                    read1,
                    read2
                ]
            )
        }

    SALMON_QUANT(
        salmon_reads_ch,
        salmon_index_ch,
        "salmon_quant"
    )

    salmon_quant_dirs_ch = SALMON_QUANT.out.quant_dir
        .map { meta, quant_dir ->
            quant_dir
        }
        .collect()

    emit:
    quant_dirs = salmon_quant_dirs_ch
    quant      = SALMON_QUANT.out.quant

    salmon_index_log      = SALMON_INDEX.out.log
    salmon_quant_log      = SALMON_QUANT.out.log
    salmon_index_versions = SALMON_INDEX.out.versions
    salmon_quant_versions = SALMON_QUANT.out.versions
}
