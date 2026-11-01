nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ABUNDANCE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Inputs:
        counts_raw.tsv
        taxonomy.tsv
        validated metadata channel

    Steps:
        ABUNDANCE_PREPROCESSING_R
        TAXONOMIC_COMPOSITION
        TAXONOMIC_COMPOSITION_GROUPED [optional]

    Outputs:
        raw counts
        filtered counts
        relative abundance
        prevalence report
        abundance preprocessing summary
        R serialized object
        taxonomic composition
        mean taxon abundance
        taxonomic composition summary

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    ABUNDANCE_PREPROCESSING_R
} from '../../../modules/local/abundance_preprocessing'


include {
    TAXONOMIC_COMPOSITION;
    TAXONOMIC_COMPOSITION_GROUPED
} from '../../../modules/local/taxonomic_composition'


workflow ABUNDANCE {

    take:
    counts
    taxonomy
    metadata


    main:

    /*
     * Prevalence filtering and relative-abundance
     * transformation.
     */

    ABUNDANCE_PREPROCESSING_R(
        counts
    )


    /*
     * Taxonomic composition by sample.
     */

    TAXONOMIC_COMPOSITION(
        ABUNDANCE_PREPROCESSING_R.out.relative_abundance,
        taxonomy
    )


    /*
     * Optional grouped taxonomic composition.
     *
     * Runs only when both metadata and a grouping
     * variable have been provided.
     */

    if (
        params.metadata &&
        params.composition_group_column
    ) {

        TAXONOMIC_COMPOSITION_GROUPED(
            TAXONOMIC_COMPOSITION.out.composition,
            metadata
        )
    }


    emit:

    /*
     * Preserve raw taxonomically filtered ASV counts.
     *
     * These are used by ALPHA_DIVERSITY.
     */

    counts_raw =
        counts

    filtered_counts =
        ABUNDANCE_PREPROCESSING_R.out.filtered_counts

    relative_abundance =
        ABUNDANCE_PREPROCESSING_R.out.relative_abundance

    prevalence_report =
        ABUNDANCE_PREPROCESSING_R.out.prevalence_report

    abundance_summary =
        ABUNDANCE_PREPROCESSING_R.out.summary

    abundance_rds =
        ABUNDANCE_PREPROCESSING_R.out.rds

    composition =
        TAXONOMIC_COMPOSITION.out.composition

    taxon_mean_abundance =
        TAXONOMIC_COMPOSITION.out.mean_abundance

    composition_summary =
        TAXONOMIC_COMPOSITION.out.summary
}
