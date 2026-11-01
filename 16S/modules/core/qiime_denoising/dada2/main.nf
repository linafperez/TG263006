nextflow.enable.dsl = 2

process DADA2 {

    tag "QIIME2 DADA2"

    label "qiime"

    publishDir "${params.outdir}/Qiime/dada2",
        mode: "copy"

    input:
    path sequences

    output:
    path "feature_table.qza", emit: table
    path "rep_seqs.qza", emit: representative_sequences
    path "denoising_stats.qza", emit: stats

    script:

    """
    qiime dada2 denoise-paired \\
        --i-demultiplexed-seqs ${sequences} \\
        --p-trunc-len-f ${params.trunc_len_f} \\
        --p-trunc-len-r ${params.trunc_len_r} \\
        --p-trim-left-f ${params.trim_left_f} \\
        --p-trim-left-r ${params.trim_left_r} \\
        --p-max-ee-f ${params.max_ee_f} \\
        --p-max-ee-r ${params.max_ee_r} \\
        --p-trunc-q ${params.trunc_q} \\
        --p-min-overlap ${params.min_overlap} \\
        --p-pooling-method ${params.pooling_method} \\
        --p-chimera-method ${params.chimera_method} \\
        --p-n-threads ${task.cpus} \\
        --p-n-reads-learn ${params.n_reads_learn} \\
        --p-hashed-feature-ids ${params.hashed_feature_ids} \\
        --p-retain-all-samples ${params.retain_all_samples} \\
        --o-table feature_table.qza \\
        --o-representative-sequences rep_seqs.qza \\
        --o-denoising-stats denoising_stats.qza
    """
}

workflow QIIME_DADA2 {

    take:
    sequences

    main:
    DADA2(sequences)

    emit:
    table = DADA2.out.table
    rep_seqs = DADA2.out.representative_sequences
    stats = DADA2.out.stats
}
