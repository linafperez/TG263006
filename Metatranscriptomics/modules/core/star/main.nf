process STAR {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_7' :
        'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7' }"
    publishDir { "${params.outdir}/rrna_and_host_removal/star/${meta.id}" }, mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path star_index
    val outdir_name

    output:
    tuple val(meta), path("${meta.id}_non_human_R1.fastq.gz"), path("${meta.id}_non_human_R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_Log.final.out"), emit: log
    tuple val(meta), path("${meta.id}_versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read1 = reads[0]
    def read2 = reads[1]

    """
    STAR \
        ${args} \
        --runMode alignReads \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --readFilesIn ${read1} ${read2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${meta.id}_ \
        --outReadsUnmapped Fastx \
        --outSAMtype None

    gzip -f ${meta.id}_Unmapped.out.mate1
    gzip -f ${meta.id}_Unmapped.out.mate2

    mv ${meta.id}_Unmapped.out.mate1.gz ${meta.id}_non_human_R1.fastq.gz
    mv ${meta.id}_Unmapped.out.mate2.gz ${meta.id}_non_human_R2.fastq.gz

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_non_human_R1.fastq.gz
    touch ${meta.id}_non_human_R2.fastq.gz
    touch ${meta.id}_Log.final.out

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        star: 2.7.11b
    END_VERSIONS
    """
}
