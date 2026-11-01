process VERSIONS_MULTIQC {

    tag "pipeline_versions"
    label 'process_low'

    input:
    path version_files, stageAs: 'versions??/*'

    output:
    path "pipeline_mqc_versions.yml", emit: yaml

    script:
    '''
    set -euo pipefail

    awk '
    function canonical(key, lower) {
        lower = tolower(key)

        if (lower == "fastqc") return "FastQC"
        if (lower == "fastp") return "fastp"
        if (lower == "sortmerna") return "SortMeRNA"
        if (lower == "star") return "STAR"
        if (lower == "spades") return "rnaSPAdes"
        if (lower == "rnaquast") return "rnaQUAST"
        if (lower == "transdecoder") return "TransDecoder"
        if (lower == "diamond") return "DIAMOND"
        if (lower == "salmon") return "Salmon"
        if (lower == "python") return "Python"
        if (lower == "r") return "R"
        if (lower == "tximport") return "tximport"
        if (lower == "deseq2") return "DESeq2"
        if (lower == "apeglm") return "apeglm"
        if (lower == "pcatools") return "PCAtools"
        if (lower == "ggplot2") return "ggplot2"
        if (lower == "pheatmap") return "pheatmap"

        return key
    }

    {
        line = $0

        sub(/^[[:space:]]+/, "", line)
        sub(/[[:space:]]+$/, "", line)

        if (line == "") next
        if (line ~ /^#/) next
        if (line == "END_VERSIONS") next
        if (line ~ /:$/) next

        pos = index(line, ":")
        if (pos == 0) next

        key = substr(line, 1, pos - 1)
        value = substr(line, pos + 1)

        sub(/^[[:space:]]+/, "", key)
        sub(/[[:space:]]+$/, "", key)

        sub(/^[[:space:]]+/, "", value)
        sub(/[[:space:]]+$/, "", value)

        if (key !~ /^[A-Za-z0-9_.+-]+$/) next

        lower_value = tolower(value)

        if (lower_value == "stub") next
        if (lower_value == "unknown") next
        if (lower_value == "na") next
        if (lower_value == "n/a") next
        if (lower_value == "none") next

        if (match(value, /[0-9][0-9]*[.][0-9][0-9A-Za-z._+-]*/)) {
            version = substr(value, RSTART, RLENGTH)
            print canonical(key) "|" version
        }
    }
    ' versions*/* > parsed_versions.txt

    echo 'MultiQC|1.25.2' >> parsed_versions.txt

    sort -u parsed_versions.txt > parsed_versions_sorted.txt

    awk '
    BEGIN {
        FS = "|"
    }

    !seen[$1]++ {
        print $1 ": " $2
    }
    ' parsed_versions_sorted.txt > pipeline_mqc_versions.yml

    echo "Generated versions:"
    cat pipeline_mqc_versions.yml
    '''
}
