process TRANSDECODER_MULTIQC {

    tag "${assembly_id}"
    label 'process_low'

    input:
    path pep
    path complete_pep
    val assembly_id

    output:
    path "${assembly_id}_transdecoder_summary_mqc.yaml", emit: yaml

    script:
    """
    predicted_orfs=\$(grep -c '^>' ${pep} || true)
    complete_orfs=\$(grep -c '^>' ${complete_pep} || true)

    awk '
        /^>/ {
            if (seen) print len
            seen = 1
            len = 0
            next
        }
        {
            gsub(/[[:space:]]/, "")
            gsub(/\\*/, "")
            len += length(\$0)
        }
        END {
            if (seen) print len
        }
    ' ${complete_pep} > complete_protein_lengths.txt

    total_aa=\$(awk '{s += \$1} END {print s + 0}' complete_protein_lengths.txt)

    longest_protein=\$(sort -nr complete_protein_lengths.txt | head -n 1)
    longest_protein=\${longest_protein:-0}

    retention_pct=\$(awk \
        -v complete="\$complete_orfs" \
        -v predicted="\$predicted_orfs" \
        'BEGIN {
            if (predicted > 0)
                printf "%.2f", (complete * 100.0 / predicted)
            else
                printf "0.00"
        }'
    )

    cat > ${assembly_id}_transdecoder_summary_mqc.yaml << EOF_MQC
id: "transdecoder_orf_summary"
section_name: "TransDecoder - ORF Prediction"
description: "Summary of predicted proteins and the complete-ORF subset retained by the pipeline."
plot_type: "table"
pconfig:
  id: "transdecoder_orf_summary_table"
  title: "TransDecoder ORF Summary"
  col1_header: "Assembly"
headers:
  predicted_orfs:
    title: "Predicted ORFs"
  complete_orfs:
    title: "Complete ORFs"
  retention_pct:
    title: "Complete (%)"
    suffix: "%"
  total_aa:
    title: "Total amino acids"
  longest_protein:
    title: "Longest protein"
data:
  "${assembly_id}":
    predicted_orfs: \${predicted_orfs}
    complete_orfs: \${complete_orfs}
    retention_pct: \${retention_pct}
    total_aa: \${total_aa}
    longest_protein: \${longest_protein}
EOF_MQC
    """

    stub:
    """
    cat > ${assembly_id}_transdecoder_summary_mqc.yaml << EOF_MQC
id: "transdecoder_orf_summary"
section_name: "TransDecoder - ORF Prediction"
plot_type: "table"
pconfig:
  id: "transdecoder_orf_summary_table"
  title: "TransDecoder ORF Summary"
  col1_header: "Assembly"
headers:
  predicted_orfs:
    title: "Predicted ORFs"
  complete_orfs:
    title: "Complete ORFs"
  retention_pct:
    title: "Complete (%)"
  total_aa:
    title: "Total amino acids"
  longest_protein:
    title: "Longest protein"
data:
  "${assembly_id}":
    predicted_orfs: 1
    complete_orfs: 1
    retention_pct: 100
    total_aa: 100
    longest_protein: 100
EOF_MQC
    """
}
