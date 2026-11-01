process FASTQC {

    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'staphb/fastqc:0.12.1' }"
    publishDir { "${params.outdir}/quality_control/fastqc/${outdir_name.replaceFirst("^fastqc_", "")}/${meta.id}" }, mode: 'copy'
    input:
    tuple val(meta), path(reads, arity: '2')
    val outdir_name

    output:
    tuple val(meta), path("${meta.id}_${outdir_name}_R[12]_fastqc.html", arity: '2'), emit: html
    tuple val(meta), path("${meta.id}_${outdir_name}_R[12]_fastqc.zip", arity: '2'), emit: zip
    tuple val(meta), path("${meta.id}_${outdir_name}_versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${outdir_name}"

    def suffix1 = reads[0].name.endsWith('.gz') ? '.fastq.gz' : '.fastq'
    def suffix2 = reads[1].name.endsWith('.gz') ? '.fastq.gz' : '.fastq'

    def name1 = "${prefix}_R1${suffix1}"
    def name2 = "${prefix}_R2${suffix2}"

    """
    mkdir -p fastqc_inputs

    ln -s "\$PWD/${reads[0]}" "fastqc_inputs/${name1}"
    ln -s "\$PWD/${reads[1]}" "fastqc_inputs/${name2}"

    fastqc \
        ${args} \
        --threads ${task.cpus} \
        --outdir . \
        "fastqc_inputs/${name1}" \
        "fastqc_inputs/${name2}"

    cat <<-END_VERSIONS > ${prefix}_versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}_${outdir_name}_R1_fastqc.html"
    touch "${meta.id}_${outdir_name}_R2_fastqc.html"
    touch "${meta.id}_${outdir_name}_R1_fastqc.zip"
    touch "${meta.id}_${outdir_name}_R2_fastqc.zip"

    cat <<-END_VERSIONS > ${meta.id}_${outdir_name}_versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}
