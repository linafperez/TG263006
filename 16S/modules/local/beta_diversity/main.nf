nextflow.enable.dsl = 2


process BETA_DIVERSITY {

    tag "R beta diversity"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Beta_diversity",
        mode: "copy"

    input:
    path filtered_counts

    output:
    path "zero_replaced_composition.tsv", emit: zero_replaced
    path "clr_matrix.tsv", emit: clr
    path "aitchison_distance.tsv", emit: distance
    path "pcoa_coordinates.tsv", emit: pcoa_coordinates
    path "pcoa_variance.tsv", emit: pcoa_variance
    path "beta_diversity_summary.txt", emit: summary
    path "aitchison_pcoa.pdf", emit: plot_pdf
    path "aitchison_pcoa.png", emit: plot_png
    path "beta_feature_filter_report.tsv", emit: feature_filter_report
    path "beta_excluded_samples.tsv", emit: excluded_samples


    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/beta_diversity.R \\
        ${filtered_counts} \\
        '${params.beta_zero_method}' \\
        ${params.beta_zero_fraction} \\
        ${params.beta_zero_threshold} \\
        ${params.beta_min_positive_samples}
    """
}


process BETA_DIVERSITY_PLOT {

    tag "R beta diversity visualization"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Beta_diversity",
        mode: "copy"

    input:
    path pcoa_coordinates
    path pcoa_variance
    path metadata

    output:
    path "aitchison_pcoa_by_group.pdf", emit: plot_pdf
    path "aitchison_pcoa_by_group.png", emit: plot_png

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/beta_diversity_plot.R \\
        ${pcoa_coordinates} \\
        ${pcoa_variance} \\
        ${metadata} \\
        '${params.beta_group_column}'
    """
}


process PERMANOVA {

    tag "R PERMANOVA and PERMDISP"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Beta_diversity/permanova",
        mode: "copy"

    input:
    path distance
    path metadata

    output:
    path "permanova_results.tsv", emit: permanova
    path "permdisp_results.tsv", emit: permdisp
    path "dispersion_by_sample.tsv", emit: dispersion
    path "permanova_summary.txt", emit: summary

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/permanova.R \\
        ${distance} \\
        ${metadata} \\
        '${params.permanova_column}' \\
        ${params.permanova_permutations}
    """
}
