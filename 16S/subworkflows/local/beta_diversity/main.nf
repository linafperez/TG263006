nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BETA DIVERSITY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Input:
        prevalence-filtered ASV counts

    Steps:
        BETA_DIVERSITY
            multiplicative zero replacement
            CLR transformation
            Aitchison distance
            PCoA

        BETA_DIVERSITY_PLOT [optional]

        PERMANOVA / PERMDISP [optional]

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    BETA_DIVERSITY as BETA_DIVERSITY_PROCESS;
    BETA_DIVERSITY_PLOT;
    PERMANOVA
} from '../../../modules/local/beta_diversity'


workflow BETA_DIVERSITY {

    take:
    filtered_counts
    metadata


    main:

    /*
     * Core compositional beta-diversity analysis.
     */

    BETA_DIVERSITY_PROCESS(
        filtered_counts
    )


    /*
     * Optional PCoA visualization using a
     * metadata grouping variable.
     */

    if (
        params.metadata &&
        params.beta_group_column
    ) {

        BETA_DIVERSITY_PLOT(
            BETA_DIVERSITY_PROCESS.out.pcoa_coordinates,
            BETA_DIVERSITY_PROCESS.out.pcoa_variance,
            metadata
        )
    }


    /*
     * Optional PERMANOVA + PERMDISP.
     */

    if (
        params.metadata &&
        params.permanova_column
    ) {

        if (params.permanova_permutations == null) {
            error(
                "permanova_permutations must be provided " +
                "when permanova_column is defined."
            )
        }

        if (
            params.permanova_permutations < 99
        ) {
            error(
                "permanova_permutations must be >= 99."
            )
        }


        PERMANOVA(
            BETA_DIVERSITY_PROCESS.out.distance,
            metadata
        )
    }


    emit:

    zero_replaced =
        BETA_DIVERSITY_PROCESS.out.zero_replaced

    clr =
        BETA_DIVERSITY_PROCESS.out.clr

    distance =
        BETA_DIVERSITY_PROCESS.out.distance

    pcoa_coordinates =
        BETA_DIVERSITY_PROCESS.out.pcoa_coordinates

    pcoa_variance =
        BETA_DIVERSITY_PROCESS.out.pcoa_variance

    summary =
        BETA_DIVERSITY_PROCESS.out.summary
}
