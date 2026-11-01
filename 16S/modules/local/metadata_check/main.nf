nextflow.enable.dsl = 2


process METADATA_CHECK {

    tag "Metadata validation"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/Preprocessing/metadata",
        mode: "copy"

    input:
    path metadata
    path samplesheet

    output:
    path "validated_metadata.tsv", emit: validated_metadata
    path "metadata_validation.txt", emit: report

    script:

    """
    python ${projectDir}/bin/check_metadata.py \\
        --metadata ${metadata} \\
        --samplesheet ${samplesheet} \\
        --output validated_metadata.tsv \\
        --report metadata_validation.txt
    """
}
