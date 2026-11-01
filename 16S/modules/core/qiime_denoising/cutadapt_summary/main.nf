nextflow.enable.dsl = 2

process CUTADAPT_SUMMARIZE {

    tag "QIIME2 Cutadapt quality summary"

    label "qiime"

    publishDir "${params.outdir}/QC/cutadapt",
        mode: "copy"

    input:
    path cutadapt_sequences

    output:
    path "cutadapt_summary.qzv", emit: summary

    script:
    """
    qiime demux summarize \
        --i-data ${cutadapt_sequences} \
        --o-visualization cutadapt_summary.qzv
    """
}


workflow QIIME_CUTADAPT_SUMMARY {

    take:
    cutadapt_sequences

    main:

    CUTADAPT_SUMMARIZE(
        cutadapt_sequences
    )

    emit:

    summary =
        CUTADAPT_SUMMARIZE.out.summary
}
