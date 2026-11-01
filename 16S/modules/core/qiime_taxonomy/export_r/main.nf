nextflow.enable.dsl = 2


process QIIME_EXPORT_R {

    tag "QIIME2 export for R"

    label "qiime"

    publishDir "${params.outdir}/Qiime/export",
        mode: "copy"

    input:
    path filtered_table
    path taxonomy

    output:
    path "counts_raw.tsv", emit: counts
    path "taxonomy.tsv", emit: taxonomy_table

    script:

    """
    set -euo pipefail

    mkdir -p exported_table
    mkdir -p exported_taxonomy

    qiime tools export \\
        --input-path ${filtered_table} \\
        --output-path exported_table

    qiime tools export \\
        --input-path ${taxonomy} \\
        --output-path exported_taxonomy

    biom convert \\
        -i exported_table/feature-table.biom \\
        -o feature_table.tsv \\
        --to-tsv

    tail -n +2 feature_table.tsv > counts_raw.tsv

    cp exported_taxonomy/taxonomy.tsv taxonomy.tsv
    """
}
