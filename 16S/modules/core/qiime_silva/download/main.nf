process SILVA_DOWNLOAD {

    tag "SILVA ${params.silva_version}"

    label "qiime"

    publishDir "${params.outdir}/Qiime/silva/download",
        mode: "copy"

    output:
    path "silva-sequences.qza", emit: sequences
    path "silva-taxonomy.qza", emit: taxonomy

    script:

    def species_labels = params.silva_include_species_labels ?
        "--p-include-species-labels" :
        "--p-no-include-species-labels"

    """
    qiime rescript get-silva-data \\
        --p-version ${params.silva_version} \\
        --p-target ${params.silva_target} \\
        ${species_labels} \\
        --o-silva-sequences silva-sequences.qza \\
        --o-silva-taxonomy silva-taxonomy.qza
    """
}
