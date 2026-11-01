nextflow.enable.dsl = 2


process ABUNDANCE_PREPROCESSING_R {

    tag "R abundance preprocessing"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/Abundance",
        mode: "copy"

    input:
    path counts

    output:
    path "filtered_counts.tsv", emit: filtered_counts
    path "relative_abundance.tsv", emit: relative_abundance
    path "prevalence_report.tsv", emit: prevalence_report
    path "abundance_preprocessing_summary.txt", emit: summary
    path "abundance_data.rds", emit: rds

    script:

    """
    export R_LIBS_USER="\${CONDA_PREFIX}/lib/R/library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null
    
    Rscript ${projectDir}/bin/abundance_preprocessing.R \\
        ${counts} \\
        ${params.prevalence_threshold}
    """
}
