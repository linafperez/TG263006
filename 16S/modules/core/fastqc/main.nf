process FASTQC {

    tag "$sample"

    label "process_medium"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/QC/fastqc",
        mode: "copy"

    input:
    tuple val(sample),
          path(read1),
          path(read2)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip", emit: zip

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        ${read1} \
        ${read2}
    """
}
