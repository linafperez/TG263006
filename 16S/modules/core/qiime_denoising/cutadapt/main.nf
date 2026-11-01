nextflow.enable.dsl = 2

process CUTADAPT_TRIM {

    tag "QIIME2 Cutadapt"

    label "qiime"

    publishDir "${params.outdir}/Qiime/cutadapt",
        mode: "copy"

    input:
    path demux
    path primers

    output:
    path "cutadapt_sequences.qza", emit: cutadapt_sequences
    path "cutadapt.log", emit: log

    script:
    """
    set -euo pipefail

    FRONT_F=\$(awk '
        /^>forward/ {getline; print; exit}
    ' ${primers})

    FRONT_R=\$(awk '
        /^>reverse/ {getline; print; exit}
    ' ${primers})

    if [[ -z "\${FRONT_F}" || -z "\${FRONT_R}" ]]; then
        echo "ERROR: primers FASTA must contain >forward and >reverse sequences." >&2
        exit 1
    fi

    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${demux} \
        --p-front-f "\${FRONT_F}" \
        --p-front-r "\${FRONT_R}" \
        --o-trimmed-sequences cutadapt_sequences.qza \
        --verbose \
        2>&1 | tee cutadapt.log
    """
}

workflow QIIME_CUTADAPT {

    take:
    demux
    primers

    main:

    CUTADAPT_TRIM(
        demux,
        primers
    )

    emit:

    cutadapt_sequences =
        CUTADAPT_TRIM.out.cutadapt_sequences
}
