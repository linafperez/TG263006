process SAMPLESHEET_CHECK {

    tag "$samplesheet"

    conda "${moduleDir}/environment.yml"

    input:
    path samplesheet

    output:
    path "samplesheet.valid.csv", emit: validated_samplesheet

    script:
    """
    python3 ${projectDir}/bin/check_samplesheet.py \
        "$samplesheet" \
        --base-dir "${launchDir}" \
        --output samplesheet.valid.csv
    """
}
