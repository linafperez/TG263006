nextflow.enable.dsl = 2

process DENOISING_SUMMARY {

    tag "QIIME2 denoising summary"

    label "qiime"

    publishDir "${params.outdir}/QC/dada2",
        mode: "copy"

    input:
    path feature_table
    path rep_seqs
    path denoising_stats

    output:
    path "feature_table.qzv", emit: feature_table_summary
    path "rep_seqs.qzv", emit: representative_sequences_summary
    path "denoising_stats.qzv", emit: denoising_stats_summary

    script:

    """
    set -euo pipefail

    qiime feature-table summarize \\
        --i-table ${feature_table} \\
        --o-visualization feature_table.qzv

    qiime feature-table tabulate-seqs \\
        --i-data ${rep_seqs} \\
        --o-visualization rep_seqs.qzv

    qiime metadata tabulate \\
        --m-input-file ${denoising_stats} \\
        --o-visualization denoising_stats.qzv
    """
}


workflow QIIME_DENOISING_SUMMARY {

    take:
    feature_table
    rep_seqs
    denoising_stats

    main:

    DENOISING_SUMMARY(
        feature_table,
        rep_seqs,
        denoising_stats
    )

    emit:

    feature_table_summary = DENOISING_SUMMARY.out.feature_table_summary
    rep_seqs_summary = DENOISING_SUMMARY.out.representative_sequences_summary
    denoising_stats_summary = DENOISING_SUMMARY.out.denoising_stats_summary
}
