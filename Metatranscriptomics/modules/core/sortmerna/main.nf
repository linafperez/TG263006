process SORTMERNA {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/919d9c8f5f2c3221a94efe96b81bde0c953c13ebb0a1eca6b690b90666006cad/data' :
        'community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e' }"
    publishDir { "${params.outdir}/rrna_and_host_removal/sortmerna/${meta.id}" }, mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path rrna_db_dir
    val rrna_idx_dir
    val outdir_name

    output:
    tuple val(meta), path("${meta.id}_non_rRNA_R1.fastq.gz"), path("${meta.id}_non_rRNA_R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_rRNA_R1.fastq.gz"), path("${meta.id}_rRNA_R2.fastq.gz"), emit: rrna
    tuple val(meta), path("${meta.id}_sortmerna.log"), emit: log
    tuple val(meta), path("${meta.id}_sortmerna_stdout.log"), emit: stdout_log
    tuple val(meta), path("${meta.id}_versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read1 = reads[0]
    def read2 = reads[1]

    """
    mkdir -p sortmerna_work

    sortmerna \
        ${args} \
        --ref ${rrna_db_dir}/smr_v4.3_fast_db.fasta \
        --ref ${rrna_db_dir}/smr_v4.3_default_db.fasta \
        --ref ${rrna_db_dir}/smr_v4.3_sensitive_db.fasta \
        --ref ${rrna_db_dir}/smr_v4.3_sensitive_db_rfam_seeds.fasta \
        --reads ${read1} \
        --reads ${read2} \
        --workdir sortmerna_work \
        --idx-dir ${rrna_idx_dir} \
        --aligned ${meta.id}_rRNA \
        --other ${meta.id}_non_rRNA \
        --threads ${task.cpus} \
        > ${meta.id}_sortmerna_stdout.log 2>&1

    native_log=""

    while IFS= read -r candidate; do
        if grep -q 'Minimal SW score based on E-value' "\$candidate"; then
            native_log="\$candidate"
            break
        fi
    done < <(
        find . \
            -maxdepth 4 \
            -type f \
            -name '*.log' \
            ! -name '${meta.id}_sortmerna_stdout.log' \
            -print
    )

    if [ -z "\$native_log" ]; then
        echo "ERROR: SortMeRNA statistics log was not generated" >&2
        echo "Log files found:" >&2
        find . -maxdepth 4 -type f -name '*.log' -print >&2
        exit 1
    fi

    cp "\$native_log" ${meta.id}_sortmerna.log

    mv ${meta.id}_non_rRNA_fwd.fq.gz ${meta.id}_non_rRNA_R1.fastq.gz
    mv ${meta.id}_non_rRNA_rev.fq.gz ${meta.id}_non_rRNA_R2.fastq.gz

    mv ${meta.id}_rRNA_fwd.fq.gz ${meta.id}_rRNA_R1.fastq.gz
    mv ${meta.id}_rRNA_rev.fq.gz ${meta.id}_rRNA_R2.fastq.gz

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        sortmerna: \$(sortmerna --version 2>&1 | grep -Eo '[0-9]+[.][0-9]+([.][0-9]+)?([A-Za-z0-9._+-]+)?' | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_non_rRNA_R1.fastq.gz
    touch ${meta.id}_non_rRNA_R2.fastq.gz
    touch ${meta.id}_rRNA_R1.fastq.gz
    touch ${meta.id}_rRNA_R2.fastq.gz
    touch ${meta.id}_sortmerna.log
    touch ${meta.id}_sortmerna_stdout.log

    cat <<-END_VERSIONS > ${meta.id}_versions.yml
    "${task.process}":
        sortmerna: 4.3.7
    END_VERSIONS
    """
}
