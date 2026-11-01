process TRANSDECODER {

    tag "${assembly_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.7.1--pl5321hdfd78af_2' :
        'quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_2' }"

    publishDir { "${params.outdir}/transcript_prediction/transdecoder/${outdir_name}" }, mode: 'copy'

    input:
    path transcripts
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}.transdecoder.pep", emit: pep
    path "${assembly_id}.transdecoder.cds", emit: cds
    path "${assembly_id}.transdecoder.gff3", emit: gff3
    path "${assembly_id}_transdecoder_dir", emit: transdecoder_dir
    path "${assembly_id}_transdecoder_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def longorfs_args = task.ext.longorfs_args ?: ''
    def predict_args = task.ext.predict_args ?: ''

    """
    cp ${transcripts} ${assembly_id}.fasta

    TransDecoder.LongOrfs \
        -t ${assembly_id}.fasta \
        ${longorfs_args} \
        -O ${assembly_id}_transdecoder_dir

    TransDecoder.Predict \
        -t ${assembly_id}.fasta \
        ${predict_args} \
        -O ${assembly_id}_transdecoder_dir

    cp ${assembly_id}_transdecoder_dir/${assembly_id}.fasta.transdecoder.pep ${assembly_id}.transdecoder.pep
    cp ${assembly_id}_transdecoder_dir/${assembly_id}.fasta.transdecoder.cds ${assembly_id}.transdecoder.cds
    cp ${assembly_id}_transdecoder_dir/${assembly_id}.fasta.transdecoder.gff3 ${assembly_id}.transdecoder.gff3

    echo '"${task.process}":' > ${assembly_id}_transdecoder_versions.yml
    echo "    transdecoder: \$(TransDecoder.LongOrfs --version 2>&1 | head -n 1 || echo unknown)" >> ${assembly_id}_transdecoder_versions.yml
    """

    stub:
    """
    cat > ${assembly_id}.transdecoder.pep << EOF
>${assembly_id}_fake_orf_1 type:complete
MRTLASLLALAAVAA
EOF

    cat > ${assembly_id}.transdecoder.cds << EOF
>${assembly_id}_fake_orf_1 type:complete
ATGCGTACGCTAGCCTGCTAGCCTTAGCAGCTTAA
EOF

    cat > ${assembly_id}.transdecoder.gff3 << EOF
##gff-version 3
${assembly_id}_fake_transcript_1\tTransDecoder\tCDS\t1\t39\t.\t+\t0\tID=${assembly_id}_fake_orf_1
EOF

    mkdir -p ${assembly_id}_transdecoder_dir

    echo '"${task.process}":' > ${assembly_id}_transdecoder_versions.yml
    echo '    transdecoder: 5.7.1' >> ${assembly_id}_transdecoder_versions.yml
    """
}
