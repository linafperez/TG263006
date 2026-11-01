process FILTER_DIAMOND_BEST_HITS {

    tag "${assembly_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"

    publishDir { "${params.outdir}/mag_annotation/filter_diamond_best_hits/${assembly_id}" }, mode: 'copy'

    input:
    path diamond_hits
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}_diamond_best_hits.tsv", emit: best_hits
    path "${assembly_id}_filter_diamond_best_hits_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_diamond_best_hits.py \
        --input ${diamond_hits} \
        --output ${assembly_id}_diamond_best_hits.tsv

    printf '"${task.process}":\n    python: %s\n' \
        "\$(python3 --version 2>&1 | sed 's/Python //')" \
        > ${assembly_id}_filter_diamond_best_hits_versions.yml
    """

    stub:
    """
    printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n' \
        > ${assembly_id}_diamond_best_hits.tsv

    printf '"${task.process}":\n    python: 3.11\n' \
        > ${assembly_id}_filter_diamond_best_hits_versions.yml
    """
}
