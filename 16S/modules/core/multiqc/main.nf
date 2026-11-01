nextflow.enable.dsl = 2

process MULTIQC {

    tag "MultiQC"

    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data' :
        'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b' }"

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "*.html", emit: report
    path "*_data", emit: data
    path "*_plots", optional: true, emit: plots

    script:
    """
    multiqc \
        --force \
        --filename multiqc_report.html \
        .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_data/.stub
    mkdir -p multiqc_plots
    touch multiqc_report.html
    """
}
