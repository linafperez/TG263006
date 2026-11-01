nextflow.enable.dsl = 2

process TRAIN_CLASSIFIER {

    tag "QIIME2 Naive Bayes classifier"

    label "qiime"

    publishDir "${params.outdir}/Qiime/classifier",
        mode: "copy"

    input:
    path reference_reads
    path reference_taxonomy

    output:
    path "classifier.qza", emit: classifier

    script:

    """
    set -euo pipefail

    qiime feature-classifier fit-classifier-naive-bayes \\
        --i-reference-reads ${reference_reads} \\
        --i-reference-taxonomy ${reference_taxonomy} \\
        --o-classifier classifier.qza
    """
}


workflow CLASSIFIER_TRAINING {

    take:
    reference_reads
    reference_taxonomy

    main:
    TRAIN_CLASSIFIER(
        reference_reads,
        reference_taxonomy
    )

    emit:
    classifier = TRAIN_CLASSIFIER.out.classifier
}
