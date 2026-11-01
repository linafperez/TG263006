nextflow.enable.dsl = 2

process SILVA_EXTRACTION {

    tag "SILVA amplicon extraction"

    label "qiime"

    publishDir "${params.outdir}/Qiime/silva/extraction",
        mode: "copy"

    input:
    path sequences
    path primers

    output:
    path "reference_reads.qza", emit: reference_reads

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

    qiime feature-classifier extract-reads \
        --i-sequences ${sequences} \
        --p-f-primer "\${FRONT_F}" \
        --p-r-primer "\${FRONT_R}" \
        --p-min-length ${params.silva_min_length} \
        --p-max-length ${params.silva_max_length} \
        --p-identity ${params.silva_extract_identity} \
        --p-n-jobs ${task.cpus} \
        --o-reads reference_reads.qza
    """
}
