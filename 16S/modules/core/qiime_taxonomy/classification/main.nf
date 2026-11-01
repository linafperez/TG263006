nextflow.enable.dsl = 2


process TAXONOMIC_CLASSIFICATION {

    tag "QIIME2 taxonomic classification"

    label "qiime"

    publishDir "${params.outdir}/Qiime/taxonomy",
        mode: "copy"

    input:
    path representative_sequences
    path classifier

    output:
    path "taxonomy.qza", emit: taxonomy
    path "taxonomy.qzv", emit: taxonomy_visualization

    script:

    """
    qiime feature-classifier classify-sklearn \\
        --i-classifier ${classifier} \\
        --i-reads ${representative_sequences} \\
        --p-n-jobs ${task.cpus} \\
        --p-read-orientation ${params.taxonomy_read_orientation} \\
        --p-confidence ${params.taxonomy_confidence} \\
        --o-classification taxonomy.qza

    qiime metadata tabulate \\
        --m-input-file taxonomy.qza \\
        --o-visualization taxonomy.qzv
    """
}


workflow QIIME_TAXONOMY {

    take:
    representative_sequences
    classifier

    main:

    TAXONOMIC_CLASSIFICATION(
        representative_sequences,
        classifier
    )

    emit:

    taxonomy = TAXONOMIC_CLASSIFICATION.out.taxonomy
    taxonomy_visualization = TAXONOMIC_CLASSIFICATION.out.taxonomy_visualization
}
