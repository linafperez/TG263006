nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QIIME2 DENOISING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Inputs:
        demultiplexed paired-end sequences
        primer sequences in FASTA format

    Steps:
        CUTADAPT_TRIM
        CUTADAPT_SUMMARIZE
        DADA2
        DENOISING_SUMMARY

    Outputs:
        cutadapt_sequences
        cutadapt_summary
        feature_table
        representative_sequences
        denoising_stats
        feature_table_summary
        representative_sequences_summary
        denoising_stats_summary

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    CUTADAPT_TRIM
} from '../../../modules/core/qiime_denoising/cutadapt'


include {
    CUTADAPT_SUMMARIZE
} from '../../../modules/core/qiime_denoising/cutadapt_summary'


include {
    DADA2
} from '../../../modules/core/qiime_denoising/dada2'


include {
    DENOISING_SUMMARY
} from '../../../modules/core/qiime_denoising/summarize'


workflow QIIME_DENOISING {

    take:
    demux
    primers


    main:

    CUTADAPT_TRIM(
        demux,
        primers
    )


    CUTADAPT_SUMMARIZE(
        CUTADAPT_TRIM.out.cutadapt_sequences
    )


    DADA2(
        CUTADAPT_TRIM.out.cutadapt_sequences
    )


    DENOISING_SUMMARY(
        DADA2.out.table,
        DADA2.out.representative_sequences,
        DADA2.out.stats
    )


    emit:

    cutadapt_sequences =
        CUTADAPT_TRIM.out.cutadapt_sequences

    cutadapt_log =
        CUTADAPT_TRIM.out.log

    cutadapt_summary =
        CUTADAPT_SUMMARIZE.out.summary

    table =
        DADA2.out.table

    rep_seqs =
        DADA2.out.representative_sequences

    stats =
        DADA2.out.stats

    feature_table_summary =
        DENOISING_SUMMARY.out.feature_table_summary

    rep_seqs_summary =
        DENOISING_SUMMARY.out.representative_sequences_summary

    denoising_stats_summary =
        DENOISING_SUMMARY.out.denoising_stats_summary
}
