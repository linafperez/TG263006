nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QIIME2 TAXONOMY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Inputs:
        SILVA dereplicated reference sequences
        SILVA dereplicated taxonomy
        ASV feature table
        representative ASV sequences

    Steps:
        SILVA_EXTRACTION
        TRAIN_CLASSIFIER
        TAXONOMIC_CLASSIFICATION
        TAXONOMIC_FILTERING
        QIIME_EXPORT_R

    Outputs:
        reference_reads
        classifier
        taxonomy
        taxonomy_visualization
        filtered_table
        filtered_sequences
        counts_raw
        taxonomy_table

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    SILVA_EXTRACTION
} from '../../../modules/core/qiime_taxonomy/extracting'


include {
    TRAIN_CLASSIFIER
} from '../../../modules/core/qiime_taxonomy/training'


include {
    TAXONOMIC_CLASSIFICATION
} from '../../../modules/core/qiime_taxonomy/classification'


include {
    TAXONOMIC_FILTERING
} from '../../../modules/core/qiime_taxonomy/filtering'


include {
    QIIME_EXPORT_R
} from '../../../modules/core/qiime_taxonomy/export_r'


workflow QIIME_TAXONOMY {

    take:
    reference_sequences
    reference_taxonomy
    feature_table
    representative_sequences
    primers


    main:

    /*
     * Extract the amplicon region from SILVA.
     */

    SILVA_EXTRACTION(
        reference_sequences,
        primers
    )


    /*
     * Train the Naive Bayes classifier.
     */

    TRAIN_CLASSIFIER(
        SILVA_EXTRACTION.out.reference_reads,
        reference_taxonomy
    )


    /*
     * Assign taxonomy to representative ASVs.
     */

    TAXONOMIC_CLASSIFICATION(
        representative_sequences,
        TRAIN_CLASSIFIER.out.classifier
    )


    /*
     * Remove unwanted taxonomic groups from
     * both feature table and representative sequences.
     */

    TAXONOMIC_FILTERING(
        feature_table,
        representative_sequences,
        TAXONOMIC_CLASSIFICATION.out.taxonomy
    )


    /*
     * Export the final taxonomically filtered table
     * and taxonomy to standard formats for R.
     */

    QIIME_EXPORT_R(
        TAXONOMIC_FILTERING.out.filtered_table,
        TAXONOMIC_CLASSIFICATION.out.taxonomy
    )


    emit:

    reference_reads =
        SILVA_EXTRACTION.out.reference_reads

    classifier =
        TRAIN_CLASSIFIER.out.classifier

    taxonomy =
        TAXONOMIC_CLASSIFICATION.out.taxonomy

    taxonomy_visualization =
        TAXONOMIC_CLASSIFICATION.out.taxonomy_visualization

    filtered_table =
        TAXONOMIC_FILTERING.out.filtered_table

    filtered_sequences =
        TAXONOMIC_FILTERING.out.filtered_sequences

    counts_raw =
        QIIME_EXPORT_R.out.counts

    taxonomy_table =
        QIIME_EXPORT_R.out.taxonomy_table
}
