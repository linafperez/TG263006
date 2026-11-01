#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE UTILITIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {
    PIPELINE_INITIALISATION
} from './subworkflows/local/utils_16s_rrna_pipeline'

include {
    PIPELINE_COMPLETION
} from './subworkflows/local/utils_16s_rrna_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN SCIENTIFIC WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {
    RRNA_16S_ANALYSIS
} from './workflows/rrna_16s_analysis'


workflow {

    main:

    /*
     * Pipeline initialisation:
     * parameter validation and execution information.
     */

    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )


    /*
     * Main 16S rRNA analysis.
     */

    RRNA_16S_ANALYSIS()


    /*
     * Pipeline completion.
     */

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        RRNA_16S_ANALYSIS.out.multiqc_report
    )
}
