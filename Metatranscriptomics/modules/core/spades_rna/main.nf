process SPADES_RNA {

    tag "${assembly_id}"
    label 'process_very_high'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/spades:4.2.0--h8d6e82b_2' :
        'staphb/spades:4.2.0' }"
    publishDir { "${params.outdir}/assembly/spades_rna/${assembly_id}" }, mode: 'copy'
    input:
    path reads
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}_transcripts.fasta", emit: transcripts
    path "${assembly_id}_spades.log", emit: log
    path "${assembly_id}_dataset.yaml", emit: dataset
    path "${assembly_id}_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory ? task.memory.toGiga() as int : 32

    """
    echo "[" > ${assembly_id}_dataset.yaml

    first=1

    for read1 in *_non_human_R1.fastq.gz; do
        [ -e "\$read1" ] || continue

        read2="\${read1/_non_human_R1.fastq.gz/_non_human_R2.fastq.gz}"

        if [ ! -f "\$read2" ]; then
            echo "ERROR: corresponding R2 file not found for \$read1" >&2
            echo "Expected: \$read2" >&2
            exit 1
        fi

        if [ "\$first" -eq 0 ]; then
            echo "," >> ${assembly_id}_dataset.yaml
        fi

        cat >> ${assembly_id}_dataset.yaml << EOF
  {
    orientation: "fr",
    type: "paired-end",
    left reads: [
      "\$read1"
    ],
    right reads: [
      "\$read2"
    ]
  }
EOF

        first=0
    done

    echo "]" >> ${assembly_id}_dataset.yaml

    echo "Dataset file:"
    cat ${assembly_id}_dataset.yaml

    ulimit -n 65535

    spades.py \
        ${args} \
        --dataset ${assembly_id}_dataset.yaml \
        --threads ${task.cpus} \
        --memory ${memory_gb} \
        -o spades_output \
        > ${assembly_id}_spades.log 2>&1

    cp spades_output/transcripts.fasta ${assembly_id}_transcripts.fasta

    echo '"${task.process}":' > ${assembly_id}_versions.yml
    echo "    spades: \$(spades.py --version 2>&1 | sed 's/SPAdes genome assembler v//')" >> ${assembly_id}_versions.yml
    """

    stub:
    """
    cat > ${assembly_id}_transcripts.fasta << EOF
>${assembly_id}_fake_transcript_1
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAACGATCGATCGATCGTAA
EOF

    touch ${assembly_id}_spades.log

    cat > ${assembly_id}_dataset.yaml << EOF
[
  {
    orientation: "fr",
    type: "paired-end",
    left reads: [
      "fake_R1.fastq.gz"
    ],
    right reads: [
      "fake_R2.fastq.gz"
    ]
  }
]
EOF

    echo '"${task.process}":' > ${assembly_id}_versions.yml
    echo '    spades: 4.2.0' >> ${assembly_id}_versions.yml
    """
}
