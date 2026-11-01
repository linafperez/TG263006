process MULTIQC {

    tag "multiqc"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.2--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/reporting_qc/multiqc" }, mode: 'copy'

    input:
    path report_files

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
    path "multiqc_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    cat > .multiqc_config.yml <<'EOF_MQC'
sample_names_replace_regex: true
sample_names_replace_exact: true

# Pipeline supplies a curated *_mqc_versions.yml file.
# Avoid duplicate / incorrectly parsed versions from module logs.
disable_version_detection: true
sample_names_replace:
  '(.+)_fastqc_raw_(R[12])[.]fastq[.]gz': '\\1_\\2'
  '(.+)_fastqc_clean_(R[12])[.]fastq[.]gz': '\\1_\\2'
  '(.+)_R1[.]fastq[.]gz_.+_R2[.]fastq[.]gz': '\\1'
  '(.+)_R[12][.]trimmed[.]fastq[.]gz': '\\1'
  '(.+)_Log[.]final[.]out': '\\1'
  '[.]': 'coassembly'

custom_data:

  expression_analysis_summary:
    section_name: "Summary"
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"
    parent_description: "Differential expression analysis and diagnostic plots for the configured case-control comparison."

  pca_plot:
    section_name: "PCA"
    description: "Principal component analysis of VST-normalized expression."
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"

  sample_distance_heatmap:
    section_name: "Sample Distance Heatmap"
    description: "Pairwise distances between samples using VST-normalized expression."
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"

  dispersion_plot:
    section_name: "Dispersion Plot"
    description: "DESeq2 dispersion estimates used by the negative-binomial model."
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"

  deg_heatmap:
    section_name: "DEG Heatmap"
    description: "Expression heatmap for differentially expressed features."
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"

  volcano_plot:
    section_name: "Volcano Plot"
    description: "Log2 fold change versus FDR-adjusted statistical significance."
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"

  ma_plot:
    section_name: "MA Plot"
    description: "Unshrunken and shrunken log2 fold changes across mean expression."
    parent_id: "expression_analysis"
    parent_name: "Expression Analysis - DESeq2"

custom_content:
  order:
    - spades_assembly_summary
    - rnaquast_assembly_summary
    - transdecoder_orf_summary
    - diamond_mapping_summary
    - expression_analysis_summary
    - pca_plot
    - sample_distance_heatmap
    - dispersion_plot
    - deg_heatmap
    - volcano_plot
    - ma_plot

module_order:

  - fastqc:
      name: "FastQC - Raw Reads"
      anchor: "fastqc_raw"
      info: "Quality assessment of raw sequencing reads before preprocessing."
      path_filters:
        - "*_fastqc_raw_R*_fastqc.zip"

  - fastp:
      name: "fastp - Read Preprocessing"
      anchor: "fastp"

  - fastqc:
      name: "FastQC - Clean Reads"
      anchor: "fastqc_clean"
      info: "Quality assessment of reads after preprocessing with fastp."
      path_filters:
        - "*_fastqc_clean_R*_fastqc.zip"

  - sortmerna:
      name: "SortMeRNA - rRNA Removal"
      anchor: "sortmerna"

  - star:
      name: "STAR - Host Read Removal"
      anchor: "star"

  - custom_content

  - salmon:
      name: "Salmon - Transcript Quantification"
      anchor: "salmon"

report_section_order:
  fastqc_raw:
    order: 1100

  fastp:
    order: 1000

  fastqc_clean:
    order: 900

  sortmerna:
    order: 800

  star:
    order: 700

  spades_assembly_summary:
    order: 600

  rnaquast_assembly_summary:
    order: 500

  transdecoder_orf_summary:
    order: 400

  diamond_mapping_summary:
    order: 300

  salmon:
    order: 200

  expression_analysis:
    order: 100
EOF_MQC

    multiqc \
        ${args} \
        --config .multiqc_config.yml \
        . \
        --outdir .

    echo '"MULTIQC":' > multiqc_versions.yml
    echo "    multiqc: \$(multiqc --version | sed 's/multiqc, version //')" >> multiqc_versions.yml
    """

    stub:
    """
    mkdir -p multiqc_data

    cat > multiqc_report.html << EOF
<html>
  <body>
    <h1>Fake MultiQC report</h1>
  </body>
</html>
EOF

    cat <<-END_VERSIONS > multiqc_versions.yml
    "${task.process}":
        multiqc: 1.25.2
    END_VERSIONS
    """
}
