process CREATE_TRANSCRIPT_MAG_GTF {

    tag "${assembly_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"

    publishDir { "${params.outdir}/mag_annotation/create_transcript_mag_gtf/${assembly_id}" }, mode: 'copy'

    input:
    path diamond_best_hits
    path transdecoder_gff3
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}_transcript_mag.gtf", emit: gtf
    path "${assembly_id}_transcript_mag_map.tsv", emit: map
    path "${assembly_id}_create_transcript_mag_gtf.log", emit: log
    path "${assembly_id}_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_transcript_mag_gtf.py \
        --diamond-best-hits ${diamond_best_hits} \
        --transdecoder-gff3 ${transdecoder_gff3} \
        --output-gtf ${assembly_id}_transcript_mag.gtf \
        --output-map ${assembly_id}_transcript_mag_map.tsv \
        > ${assembly_id}_create_transcript_mag_gtf.log 2>&1

    cat <<-END_VERSIONS > ${assembly_id}_versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '%s\tcreate_transcript_mag_gtf\tCDS\t1\t39\t.\t+\t0\tgene_id "fake_MAG_protein_1"; transcript_id "%s"; orf_id "%s"; mag_id "fake_MAG_protein_1"; diamond_target "fake_MAG_protein_1"; pident "99.0"; evalue "1e-50"; bitscore "200";\n' \
        "${assembly_id}_fake_transcript_1" \
        "${assembly_id}_fake_transcript_1" \
        "${assembly_id}_fake_orf_1" \
        > ${assembly_id}_transcript_mag.gtf

    printf 'qseqid\ttranscript_id\tgene_id\tmag_id\tsseqid\tpident\tevalue\tbitscore\n%s\t%s\tfake_MAG_protein_1\tfake_MAG_protein_1\tfake_MAG_protein_1\t99.0\t1e-50\t200\n' \
        "${assembly_id}_fake_orf_1" \
        "${assembly_id}_fake_transcript_1" \
        > ${assembly_id}_transcript_mag_map.tsv

    touch ${assembly_id}_create_transcript_mag_gtf.log

    cat <<-END_VERSIONS > ${assembly_id}_versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
