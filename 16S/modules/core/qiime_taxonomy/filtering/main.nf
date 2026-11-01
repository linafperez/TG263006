nextflow.enable.dsl = 2


process TAXONOMIC_FILTERING {

    tag "QIIME2 taxonomic filtering"

    label "qiime"

    publishDir "${params.outdir}/Qiime/taxonomic_filtering",
        mode: "copy"

    input:
    path feature_table
    path representative_sequences
    path taxonomy

    output:
    path "filtered_table.qza", emit: filtered_table
    path "filtered_rep_seqs.qza", emit: filtered_sequences

    script:

    """
    qiime taxa filter-table \\
        --i-table ${feature_table} \\
        --i-taxonomy ${taxonomy} \\
        --p-exclude '${params.taxonomy_exclude}' \\
        --o-filtered-table filtered_table.qza

    qiime taxa filter-seqs \\
        --i-sequences ${representative_sequences} \\
        --i-taxonomy ${taxonomy} \\
        --p-exclude '${params.taxonomy_exclude}' \\
        --o-filtered-sequences filtered_rep_seqs.qza
    """
}


workflow QIIME_TAXONOMIC_FILTERING {

    take:
    feature_table
    representative_sequences
    taxonomy

    main:

    TAXONOMIC_FILTERING(
        feature_table,
        representative_sequences,
        taxonomy
    )

    emit:

    filtered_table = TAXONOMIC_FILTERING.out.filtered_table
    filtered_sequences = TAXONOMIC_FILTERING.out.filtered_sequences
}
