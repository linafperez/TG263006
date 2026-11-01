nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALPHA DIVERSITY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Input:
        taxonomically filtered raw ASV counts

    Steps:
        ALPHA_DIVERSITY
        ALPHA_DIVERSITY_PLOT [optional]

    Metrics:
        Observed richness
        Chao1
        Shannon
        Simpson

    Important:
        Alpha diversity is calculated BEFORE prevalence filtering.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    ALPHA_DIVERSITY as ALPHA_DIVERSITY_PROCESS;
    ALPHA_DIVERSITY_PLOT
} from '../../../modules/local/alpha_diversity'


workflow ALPHA_DIVERSITY {

    take:
    counts
    metadata


    main:

    /*
     * Calculate alpha diversity from raw
     * taxonomically filtered ASV counts.
     */

    ALPHA_DIVERSITY_PROCESS(
        counts
    )


    /*
     * Optional visualization by metadata group.
     */

    if (
        params.metadata &&
        params.alpha_group_column
    ) {

        ALPHA_DIVERSITY_PLOT(
            ALPHA_DIVERSITY_PROCESS.out.metrics_long,
            metadata
        )
    }


    emit:

    metrics =
        ALPHA_DIVERSITY_PROCESS.out.metrics

    metrics_long =
        ALPHA_DIVERSITY_PROCESS.out.metrics_long

    summary =
        ALPHA_DIVERSITY_PROCESS.out.summary
}
