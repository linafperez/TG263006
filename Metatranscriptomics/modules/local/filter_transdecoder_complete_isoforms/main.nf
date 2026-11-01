process FILTER_TRANSDECODER_COMPLETE_ISOFORMS {

    tag "${assembly_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"

    publishDir { "${params.outdir}/transcript_prediction/filter_transdecoder_complete_isoforms/${outdir_name}" }, mode: 'copy'

    input:
    path pep
    path cds
    path gff3
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}.complete.transdecoder.pep", emit: pep
    path "${assembly_id}.complete.transdecoder.cds", emit: cds
    path "${assembly_id}.complete.transdecoder.gff3", emit: gff3
    path "${assembly_id}_filter_transdecoder_complete_isoforms.log", emit: log
    path "${assembly_id}_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_transdecoder_complete_isoforms.py \
        --pep ${pep} \
        --cds ${cds} \
        --gff3 ${gff3} \
        --output-pep ${assembly_id}.complete.transdecoder.pep \
        --output-cds ${assembly_id}.complete.transdecoder.cds \
        --output-gff3 ${assembly_id}.complete.transdecoder.gff3 \
        > ${assembly_id}_filter_transdecoder_complete_isoforms.log 2>&1

    echo '"${task.process}":' > ${assembly_id}_versions.yml
    echo "    python: \$(python --version | sed 's/Python //')" >> ${assembly_id}_versions.yml
    """

    stub:
    """
    cat > ${assembly_id}.complete.transdecoder.pep << EOF
>${assembly_id}_fake_orf_1 type:complete
MRTLASLLALAAVAA
EOF

    cat > ${assembly_id}.complete.transdecoder.cds << EOF
>${assembly_id}_fake_orf_1 type:complete
ATGCGTACGCTAGCCTGCTAGCCTTAGCAGCTTAA
EOF

    cat > ${assembly_id}.complete.transdecoder.gff3 << EOF
##gff-version 3
${assembly_id}_fake_transcript_1	TransDecoder	mRNA	1	39	.	+	.	ID=${assembly_id}_fake_orf_1
${assembly_id}_fake_transcript_1	TransDecoder	CDS	1	39	.	+	0	ID=${assembly_id}_fake_orf_1.cds;Parent=${assembly_id}_fake_orf_1
EOF

    touch ${assembly_id}_filter_transdecoder_complete_isoforms.log

    echo '"${task.process}":' > ${assembly_id}_versions.yml
    echo '    python: 3.11' >> ${assembly_id}_versions.yml
    """
}
