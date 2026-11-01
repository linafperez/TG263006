process EXPRESSION_MULTIQC {

    tag "${comparison_id}"
    label 'process_low'

    input:
    path summary
    path mag_summary
    val comparison_id

    output:
    path "${comparison_id}_expression_summary_mqc.yaml", emit: yaml

    script:
    """
    total_features=\$(awk '\$1 == "total_features" {print \$2}' ${summary})
    features_with_padj=\$(awk '\$1 == "features_with_padj" {print \$2}' ${summary})
    significant=\$(awk '\$1 == "significant_padj_lt_0.05" {print \$2}' ${summary})
    upregulated=\$(awk '\$1 == "upregulated_log2fc_gt_0" {print \$2}' ${summary})
    downregulated=\$(awk '\$1 == "downregulated_log2fc_lt_0" {print \$2}' ${summary})

    total_features=\${total_features:-0}
    features_with_padj=\${features_with_padj:-0}
    significant=\${significant:-0}
    upregulated=\${upregulated:-0}
    downregulated=\${downregulated:-0}

    mags_represented=\$(awk 'NR > 1 && NF {n++} END {print n + 0}' ${mag_summary})
    mags_significant=\$(awk 'NR > 1 && NF && \$3 + 0 > 0 {n++} END {print n + 0}' ${mag_summary})

    cat > ${comparison_id}_expression_summary_mqc.yaml << EOF_MQC
id: "expression_analysis_summary"
section_name: "Expression Analysis - DESeq2"
description: "Summary of differential expression and MAG-level results produced by the DESeq2 analysis."
plot_type: "table"
pconfig:
  id: "expression_analysis_summary_table"
  title: "Differential Expression Summary"
  col1_header: "Comparison"
headers:
  total_features:
    title: "Features tested"
  features_with_padj:
    title: "Features with padj"
  significant:
    title: "Significant"
    description: "Features with adjusted p-value < 0.05"
  upregulated:
    title: "Up-regulated"
  downregulated:
    title: "Down-regulated"
  mags_represented:
    title: "MAGs represented"
  mags_significant:
    title: "MAGs with significant features"
data:
  "${comparison_id}":
    total_features: \${total_features}
    features_with_padj: \${features_with_padj}
    significant: \${significant}
    upregulated: \${upregulated}
    downregulated: \${downregulated}
    mags_represented: \${mags_represented}
    mags_significant: \${mags_significant}
EOF_MQC
    """

    stub:
    """
    cat > ${comparison_id}_expression_summary_mqc.yaml << EOF_MQC
id: "expression_analysis_summary"
section_name: "Expression Analysis - DESeq2"
plot_type: "table"
pconfig:
  id: "expression_analysis_summary_table"
  title: "Differential Expression Summary"
  col1_header: "Comparison"
headers:
  total_features:
    title: "Features tested"
  features_with_padj:
    title: "Features with padj"
  significant:
    title: "Significant"
  upregulated:
    title: "Up-regulated"
  downregulated:
    title: "Down-regulated"
  mags_represented:
    title: "MAGs represented"
  mags_significant:
    title: "MAGs with significant features"
data:
  "${comparison_id}":
    total_features: 1
    features_with_padj: 0
    significant: 0
    upregulated: 0
    downregulated: 0
    mags_represented: 1
    mags_significant: 0
EOF_MQC
    """
}
