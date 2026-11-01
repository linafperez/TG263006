process DIAMOND {

    tag "${assembly_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.10--h43eeafb_0' :
        'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0' }"

    publishDir { "${params.outdir}/mag_annotation/diamond/${assembly_id}" }, mode: 'copy'

    input:
    path proteins
    path diamond_db
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}_diamond.daa", emit: daa
    path "${assembly_id}_diamond_hits.tsv", emit: hits
    path "diamond.log", emit: log
    path "${assembly_id}_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    diamond blastp \
        ${args} \
        --db ${diamond_db} \
        --query ${proteins} \
        --out ${assembly_id}_diamond.daa \
        --outfmt 100 \
        --threads ${task.cpus} \
        --log

    diamond view \
        --daa ${assembly_id}_diamond.daa \
        --out ${assembly_id}_diamond_hits.tsv \
        --outfmt 6 \
            qseqid \
            sseqid \
            pident \
            length \
            mismatch \
            gapopen \
            qstart \
            qend \
            sstart \
            send \
            evalue \
            bitscore

    cat <<-END_VERSIONS > ${assembly_id}_versions.yml
    "${task.process}":
        diamond: \$(diamond version | head -n 1 | sed 's/diamond version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${assembly_id}_diamond.daa

    cat > ${assembly_id}_diamond_hits.tsv << EOF2
${assembly_id}_fake_orf_1       fake_MAG_protein_1      99.0    100     0       0       1       100     1       100    1e-50    200
EOF2

    touch diamond.log

    cat <<-END_VERSIONS > ${assembly_id}_versions.yml
    "${task.process}":
        diamond: 2.1.10
    END_VERSIONS
    """
}
