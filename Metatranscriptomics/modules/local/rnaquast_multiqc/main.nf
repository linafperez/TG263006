process RNAQUAST_MULTIQC {

    tag "${assembly_id}"
    label 'process_low'

    input:
    path report_dir
    val assembly_id
    val outdir_name

    output:
    path "${assembly_id}_rnaquast_summary_mqc.yaml", emit: yaml

    script:
    """
    short_report="${report_dir}/short_report.tsv"

    if [ ! -s "\$short_report" ]; then
        echo "ERROR: rnaQUAST short_report.tsv not found or empty: \$short_report" >&2
        exit 1
    fi

    transcripts=\$(awk -F '\\t' '\$1 == "Transcripts" { print \$2; exit }' "\$short_report" | tr -d '\\r')

    transcripts_500=\$(awk -F '\\t' '\$1 == "Transcripts > 500 bp" { print \$2; exit }' "\$short_report" | tr -d '\\r')

    transcripts_1000=\$(awk -F '\\t' '\$1 == "Transcripts > 1000 bp" { print \$2; exit }' "\$short_report" | tr -d '\\r')

    for value in "\$transcripts" "\$transcripts_500" "\$transcripts_1000"; do
        case "\$value" in
            ''|*[!0-9]*)
                echo "ERROR: invalid rnaQUAST metric value: '\$value'" >&2
                exit 1
                ;;
        esac
    done

    cat > ${assembly_id}_rnaquast_summary_mqc.yaml << EOF_MQC
id: "rnaquast_assembly_summary"
section_name: "rnaQUAST Assembly Summary"
description: "Transcript assembly metrics reported by rnaQUAST."
plot_type: "table"

pconfig:
  id: "rnaquast_assembly_summary_table"
  title: "rnaQUAST Assembly Summary"
  col1_header: "Assembly"

headers:
  transcripts:
    title: "Transcripts"
  transcripts_gt_500:
    title: "Transcripts > 500 bp"
  transcripts_gt_1000:
    title: "Transcripts > 1000 bp"

data:
  ${assembly_id}:
    transcripts: \$transcripts
    transcripts_gt_500: \$transcripts_500
    transcripts_gt_1000: \$transcripts_1000
EOF_MQC
    """

    stub:
    """
    cat > ${assembly_id}_rnaquast_summary_mqc.yaml << EOF_MQC
id: "rnaquast_assembly_summary"
section_name: "rnaQUAST Assembly Summary"
description: "Transcript assembly metrics reported by rnaQUAST."
plot_type: "table"

pconfig:
  id: "rnaquast_assembly_summary_table"
  title: "rnaQUAST Assembly Summary"
  col1_header: "Assembly"

headers:
  transcripts:
    title: "Transcripts"
  transcripts_gt_500:
    title: "Transcripts > 500 bp"
  transcripts_gt_1000:
    title: "Transcripts > 1000 bp"

data:
  ${assembly_id}:
    transcripts: 1
    transcripts_gt_500: 0
    transcripts_gt_1000: 0
EOF_MQC
    """
}
