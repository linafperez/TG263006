nextflow.enable.dsl = 2


process TAXONOMIC_COMPOSITION {

    tag "R taxonomic composition"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Taxonomy_composition",
        mode: "copy"

    input:
    path relative_abundance
    path taxonomy

    output:
    path "taxonomic_composition.tsv", emit: composition
    path "taxon_mean_abundance.tsv", emit: mean_abundance
    path "taxonomic_composition_summary.txt", emit: summary
    path "taxonomic_composition_by_sample.pdf", emit: sample_plot_pdf
    path "taxonomic_composition_by_sample.png", emit: sample_plot_png

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/taxonomic_composition.R \\
        ${relative_abundance} \\
        ${taxonomy} \\
        '${params.composition_taxonomic_rank}' \\
        ${params.composition_rare_threshold}
    """
}


process TAXONOMIC_COMPOSITION_GROUPED {

    tag "R grouped taxonomic composition"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Taxonomy_composition",
        mode: "copy"

    input:
    path composition
    path metadata

    output:
    path "taxonomic_composition_by_group.tsv", emit: grouped_table
    path "taxonomic_composition_by_group.pdf", emit: grouped_plot_pdf
    path "taxonomic_composition_by_group.png", emit: grouped_plot_png

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null

    Rscript ${projectDir}/bin/taxonomic_composition_grouped.R \\
        ${composition} \\
        ${metadata} \\
        '${params.composition_group_column}'
    """
}
