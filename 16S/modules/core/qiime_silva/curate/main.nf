process SILVA_CURATION {

    tag "SILVA curation"

    label "qiime"

    publishDir "${params.outdir}/Qiime/silva/curation",
        mode: "copy"

    input:
    path sequences
    path taxonomy

    output:
    path "silva-sequences-clean.qza", emit: clean_sequences
    path "silva-sequences-derep.qza", emit: dereplicated_sequences
    path "silva-taxonomy-derep.qza", emit: dereplicated_taxonomy

    script:

    """
    set -euo pipefail

    qiime rescript cull-seqs \\
        --i-sequences ${sequences} \\
        --o-clean-sequences silva-sequences-clean.qza

    qiime rescript dereplicate \\
        --i-sequences silva-sequences-clean.qza \\
        --i-taxa ${taxonomy} \\
        --p-mode ${params.silva_derep_mode} \\
        --o-dereplicated-sequences silva-sequences-derep.qza \\
        --o-dereplicated-taxa silva-taxonomy-derep.qza
    """
}
