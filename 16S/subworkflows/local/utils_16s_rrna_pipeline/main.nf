//
// Pipeline utilities specific to the 16S rRNA analysis pipeline.
// No external dependencies required.
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version
    validate_params
    monochrome_logs
    nextflow_cli_args
    outdir
    _input
    help
    help_full
    show_hidden

    main:

    def command = """\
    nextflow run main.nf \\
        -profile singularity \\
        --input samplesheet.csv \\
        --outdir results \\
        --primers data/primers.fasta
    """.stripIndent()

    log.info """
    ========================================================================
     ${workflow.manifest.name ?: '16S rRNA Analysis'} v${workflow.manifest.version ?: '1.0.0'}
    ========================================================================
    Example command:

    ${command}
    """.stripIndent()

    if (validate_params) {
        validateInputParameters()
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    monochrome_logs
    multiqc_report

    main:

    workflow.onComplete {
        pipelineSummary()
    }

    workflow.onError {
        log.error """
        [16S rRNA Analysis] Pipeline failed.
        Please check .nextflow.log and the execution reports for details.
        """.stripIndent()
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateInputParameters() {

    /*
     * Core input
     */

    if (!params.input) {
        error "Missing input: provide a samplesheet with '--input'."
    }

    if (!params.outdir) {
        error "Missing output directory: provide '--outdir'."
    }

    if (!params.primers) {
        error "Missing primers: provide a FASTA file with '--primers'."
    }


    /*
     * SILVA local database
     *
     * Either provide both local artifacts or neither.
     */

    if (
        (params.silva_sequences && !params.silva_taxonomy) ||
        (!params.silva_sequences && params.silva_taxonomy)
    ) {
        error """
        Local SILVA mode requires both:
          --silva_sequences
          --silva_taxonomy

        Otherwise omit both parameters to use automatic SILVA download.
        """.stripIndent()
    }


    /*
     * Taxonomic filtering
     */

    if (!params.taxonomy_exclude) {
        error "Missing '--taxonomy_exclude'."
    }


    /*
     * Abundance preprocessing
     */

    if (params.prevalence_threshold == null) {
        error "Missing '--prevalence_threshold'."
    }

    if (
        params.prevalence_threshold < 0 ||
        params.prevalence_threshold > 1
    ) {
        error "'prevalence_threshold' must be between 0 and 1."
    }


    /*
     * Taxonomic composition
     */

    if (!params.composition_taxonomic_rank) {
        error "Missing '--composition_taxonomic_rank'."
    }

    def valid_composition_ranks = [
        'Domain',
        'Phylum',
        'Class',
        'Order',
        'Family',
        'Genus',
        'Species'
    ]

    if (
        !valid_composition_ranks.contains(
            params.composition_taxonomic_rank
        )
    ) {
        error(
            "Invalid composition_taxonomic_rank: " +
            "${params.composition_taxonomic_rank}. " +
            "Allowed values: " +
            "${valid_composition_ranks.join(', ')}"
        )
    }

    if (params.composition_rare_threshold == null) {
        error "Missing '--composition_rare_threshold'."
    }

    if (
        params.composition_rare_threshold < 0 ||
        params.composition_rare_threshold >= 1
    ) {
        error "'composition_rare_threshold' must be >= 0 and < 1."
    }


    /*
     * Beta diversity
     */

    if (!params.beta_zero_method) {
        error "Missing '--beta_zero_method'."
    }

    def supportedBetaZeroMethods = [
        'CZM',
        'GBM',
        'SQ',
        'BL'
    ]

    if (!(params.beta_zero_method in supportedBetaZeroMethods)) {
        error(
            "Unsupported beta_zero_method: " +
            "${params.beta_zero_method}. Supported methods: " +
            supportedBetaZeroMethods.join(', ')
        )
    }

    if (params.beta_zero_fraction == null) {
        error "Missing '--beta_zero_fraction'."
    }

    if (
        params.beta_zero_fraction <= 0 ||
        params.beta_zero_fraction >= 1
    ) {
        error "'beta_zero_fraction' must be > 0 and < 1."
    }

    if (params.beta_zero_threshold == null) {
        error "Missing '--beta_zero_threshold'."
    }

    if (params.beta_zero_threshold <= 0) {
        error "'beta_zero_threshold' must be > 0."
    }

    if (
        params.beta_min_positive_samples == null ||
        params.beta_min_positive_samples < 1
    ) {
        error(
            "beta_min_positive_samples must be >= 1."
        )
    }

    if (
        params.beta_zero_method == 'GBM' &&
        params.beta_min_positive_samples < 2
    ) {
        error(
            "GBM requires beta_min_positive_samples >= 2."
        )
    }


    /*
     * Optional PERMANOVA
     */

    if (
        params.metadata &&
        params.permanova_column
    ) {

        if (params.permanova_permutations == null) {
            error(
                "Provide '--permanova_permutations' " +
                "when '--permanova_column' is defined."
            )
        }

        if (params.permanova_permutations < 99) {
            error "'permanova_permutations' must be >= 99."
        }
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def pipelineSummary() {

    if (workflow.success) {

        log.info """
        ========================================================================
         [${workflow.manifest.name ?: '16S rRNA Analysis'}]
         Pipeline completed successfully
        ========================================================================
        Results:
          ${params.outdir}
        """.stripIndent()

    } else {

        log.info """
        ========================================================================
         [${workflow.manifest.name ?: '16S rRNA Analysis'}]
         Pipeline completed with errors
        ========================================================================
        """.stripIndent()
    }
}
