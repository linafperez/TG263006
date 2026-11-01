process IMPORT {

    tag "QIIME2 Import"
    label "process_medium"
    label "qiime"

    publishDir "${params.outdir}/Qiime/import",
        mode: "copy"

    input:
    val samples

    output:
    path "manifest.tsv", emit: manifest
    path "demux.qza", emit: demux

    script:

    def manifestRows = samples.collect { sample ->

        def id = sample[0]
        def r1 = sample[1]
        def r2 = sample[2]

        "${id}\t${r1}\t${r2}"

    }.join('\n')

    """
    cat > manifest.tsv << EOF
sample-id\tforward-absolute-filepath\treverse-absolute-filepath
${manifestRows}
EOF

    echo "===== GENERATED QIIME2 MANIFEST ====="
    cat manifest.tsv
    echo "====================================="

    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-path manifest.tsv \\
        --input-format PairedEndFastqManifestPhred33V2 \\
        --output-path demux.qza
    """
}
