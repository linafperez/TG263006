process FASTP {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'staphb/fastp:0.23.4' }"
    publishDir { "${params.outdir}/quality_control/fastp/${meta.id}" }, mode: 'copy'
    input:
    tuple val(meta), path(reads)
    val outdir_name

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_fastp.html"), emit: html
    tuple val(meta), path("${meta.id}_fastp.json"), emit: json
    tuple val(meta), path("${meta.id}_versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read1 = reads[0]
    def read2 = reads[1]

    """
    fastp \
        ${args} \
        --thread ${task.cpus} \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${meta.id}_R1.trimmed.fastq.gz \
        --out2 ${meta.id}_R2.trimmed.fastq.gz \
        --html ${meta.id}_fastp.html \
        --json ${meta.id}_fastp.json

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_R1.trimmed.fastq.gz
    touch ${meta.id}_R2.trimmed.fastq.gz
    touch ${meta.id}_fastp.html
    touch ${meta.id}_fastp.json

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}