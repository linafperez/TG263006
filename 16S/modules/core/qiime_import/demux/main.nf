nextflow.enable.dsl = 2

process DEMUX_SUMMARIZE {

    tag "QIIME2 Demux Summary"

    label "qiime"

    publishDir "${params.outdir}/QC/demux",
        mode: "copy"

    input:
    path demux

    output:
    path "demux.qzv", emit: summary

    script:

    """
    qiime demux summarize \
        --i-data ${demux} \
        --o-visualization demux.qzv
    """
}

workflow QIIME_DEMUX_SUMMARY {

    take:
    demux

    main:
    DEMUX_SUMMARIZE(demux)

    emit:
    summary = DEMUX_SUMMARIZE.out.summary
}
