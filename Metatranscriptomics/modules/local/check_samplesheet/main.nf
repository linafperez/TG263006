process CHECK_SAMPLESHEET {

    tag "samplesheet"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    input:
    path samplesheet

    output:
    path "samplesheet.valid.csv", emit: samplesheet
    path "check_samplesheet.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    WORKDIR="\$PWD"
    INPUT="\$(readlink -f ${samplesheet})"

    cd "${projectDir}"

    check_samplesheet.py \
        --input "\$INPUT" \
        --output "\$WORKDIR/samplesheet.valid.csv" \
        > "\$WORKDIR/check_samplesheet.log" 2>&1

    cd "\$WORKDIR"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${samplesheet} samplesheet.valid.csv
    touch check_samplesheet.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}

