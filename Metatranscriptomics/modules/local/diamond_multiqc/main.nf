process DIAMOND_MULTIQC {

    tag "${assembly_id}"
    label 'process_low'

    input:
    path hits
    path best_hits
    val assembly_id

    output:
    path "${assembly_id}_diamond_summary_mqc.yaml", emit: yaml

    script:
    """
    alignment_records=\$(awk '
        NF {
            first = tolower(\$1)
            if (first != "qseqid" && first != "query" && first != "query_id")
                n++
        }
        END {print n + 0}
    ' ${hits})

    queries_aligned=\$(awk '
        NF {
            first = tolower(\$1)
            if (first != "qseqid" && first != "query" && first != "query_id")
                print \$1
        }
    ' ${hits} | sort -u | wc -l | tr -d ' ')

    best_hit_records=\$(awk '
        NF {
            first = tolower(\$1)
            if (first != "qseqid" && first != "query" && first != "query_id")
                n++
        }
        END {print n + 0}
    ' ${best_hits})

    cat > ${assembly_id}_diamond_summary_mqc.yaml << EOF_MQC
id: "diamond_mapping_summary"
section_name: "DIAMOND - Transcript-to-MAG Mapping"
description: "Summary of transcript-to-MAG protein alignments and best hits retained by the pipeline."
plot_type: "table"
pconfig:
  id: "diamond_mapping_summary_table"
  title: "DIAMOND Mapping Summary"
  col1_header: "Assembly"
headers:
  queries_aligned:
    title: "Queries aligned"
  alignment_records:
    title: "Alignment records"
  best_hit_records:
    title: "Best hits retained"
data:
  "${assembly_id}":
    queries_aligned: \${queries_aligned}
    alignment_records: \${alignment_records}
    best_hit_records: \${best_hit_records}
EOF_MQC
    """

    stub:
    """
    cat > ${assembly_id}_diamond_summary_mqc.yaml << EOF_MQC
id: "diamond_mapping_summary"
section_name: "DIAMOND - Transcript-to-MAG Mapping"
plot_type: "table"
pconfig:
  id: "diamond_mapping_summary_table"
  title: "DIAMOND Mapping Summary"
  col1_header: "Assembly"
headers:
  queries_aligned:
    title: "Queries aligned"
  alignment_records:
    title: "Alignment records"
  best_hit_records:
    title: "Best hits retained"
data:
  "${assembly_id}":
    queries_aligned: 1
    alignment_records: 1
    best_hit_records: 1
EOF_MQC
    """
}
