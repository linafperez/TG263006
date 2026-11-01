nextflow.enable.dsl = 2


process ALPHA_DIVERSITY {

    tag "R alpha diversity"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Alpha_diversity",
        mode: "copy"

    input:
    path counts

    output:
    path "alpha_diversity.tsv", emit: metrics
    path "alpha_diversity_long.tsv", emit: metrics_long
    path "alpha_diversity_summary.txt", emit: summary

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/alpha_diversity.R \\
        ${counts}
    """
}


process ALPHA_DIVERSITY_PLOT {

    tag "R alpha diversity visualization"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Alpha_diversity",
        mode: "copy"

    input:
    path alpha_metrics
    path metadata

    output:
    path "alpha_diversity_by_group.pdf", emit: plot_pdf
    path "alpha_diversity_by_group.png", emit: plot_png

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/alpha_diversity_plot.R \\
        ${alpha_metrics} \\
        ${metadata} \\
        '${params.alpha_group_column}'
    """
}
