process SALMON_INDEX {

    tag "${index_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2' :
        'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2' }"

    publishDir { "${params.outdir}/quantification/salmon/index" },
        mode: 'copy',
        saveAs: { filename ->
            if (filename == "${index_id}_salmon_index") {
                return "${index_id}"
            }
            return "${index_id}/${filename}"
        }

    input:
    path transcripts
    val index_id
    val outdir_name

    output:
    path "${index_id}_salmon_index", emit: index
    path "${index_id}_salmon_index.log", emit: log
    path "${index_id}_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    salmon index \
        ${args} \
        -t ${transcripts} \
        -i ${index_id}_salmon_index \
        -p ${task.cpus} \
        > ${index_id}_salmon_index.log 2>&1

    cat <<-END_VERSIONS > ${index_id}_versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${index_id}_salmon_index
    touch ${index_id}_salmon_index/versionInfo.json
    touch ${index_id}_salmon_index.log

    cat <<-END_VERSIONS > ${index_id}_versions.yml
    "${task.process}":
        salmon: 1.10.3
    END_VERSIONS
    """
}


process SALMON_QUANT {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2' :
        'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2' }"

    publishDir { "${params.outdir}/quantification/salmon/quant" },
        mode: 'copy',
        saveAs: { filename ->
            if (filename == 'quant.sf' || filename.endsWith('/quant.sf')) {
                return null
            }

            if (filename == "${meta.id}" || filename.startsWith("${meta.id}/")) {
                return filename
            }

            return "${meta.id}/${filename}"
        }

    input:
    tuple val(meta), path(reads)
    path salmon_index
    val outdir_name

    output:
    tuple val(meta), path("${meta.id}/quant.sf"), emit: quant
    tuple val(meta), path("${meta.id}"), emit: quant_dir
    tuple val(meta), path("${meta.id}_salmon_quant.log"), emit: log
    tuple val(meta), path("${meta.id}_versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read1 = reads[0]
    def read2 = reads[1]

    """
    salmon quant \
        ${args} \
        -i ${salmon_index} \
        -1 ${read1} \
        -2 ${read2} \
        -p ${task.cpus} \
        -o ${meta.id} \
        > ${meta.id}_salmon_quant.log 2>&1

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}

    cat > ${meta.id}/quant.sf << EOF
Name	Length	EffectiveLength	TPM	NumReads
coassembly_fake_transcript_1	1000	800	10.5	50
EOF

    touch ${meta.id}_salmon_quant.log

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        salmon: 1.10.3
    END_VERSIONS
    """
}
