process QUAST {

    tag "${assembly_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/rnaquast:2.3.0--h9ee0642_0' :
        'nanozoo/rnaquast:1.5.2--0ecea11' }"
    publishDir { "${params.outdir}/assembly/rnaquast/${assembly_id}" }, mode: 'copy'
    input:
    path transcripts
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}_rnaquast_report", emit: report_dir
    path "${assembly_id}_rnaquast_report/short_report.txt", emit: txt
    path "${assembly_id}_rnaquast_report/short_report.pdf", optional: true, emit: pdf
    path "${assembly_id}_rnaquast_report/short_report.tex", optional: true, emit: tex
    path "${assembly_id}_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p ${assembly_id}_rnaquast_report

    export MPLCONFIGDIR=\$PWD/${assembly_id}_rnaquast_report/matplotlib_cache
    mkdir -p \$MPLCONFIGDIR

    rnaQUAST.py \
        -c ${transcripts} \
        ${args} \
        -o ${assembly_id}_rnaquast_report \
        -t ${task.cpus}

    echo '"${task.process}":' > ${assembly_id}_versions.yml
    echo '    rnaquast: 2.3.0' >> ${assembly_id}_versions.yml
    """

    stub:
    """
    mkdir -p ${assembly_id}_rnaquast_report

    cat > ${assembly_id}_rnaquast_report/short_report.txt << EOF
SHORT SUMMARY REPORT

METRICS/TRANSCRIPTS                                    ${assembly_id}
Transcripts                                            1
EOF

    cat > ${assembly_id}_rnaquast_report/short_report.tsv << EOF
METRICS/TRANSCRIPTS	${assembly_id}_transcripts
Transcripts	1
Transcripts > 500 bp	0
Transcripts > 1000 bp	0
EOF

    cat > ${assembly_id}_rnaquast_report/short_report.tex << EOF
rnaQUAST fake report for ${assembly_id}
EOF

    echo '"${task.process}":' > ${assembly_id}_versions.yml
    echo '    rnaquast: 2.3.0' >> ${assembly_id}_versions.yml
    """
}
