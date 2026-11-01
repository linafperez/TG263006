process EXPRESSION_SUMMARY {

    tag "${outdir_name}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"

    publishDir {
        "${params.outdir}/expression_analysis/summary/${outdir_name}"
    }, mode: 'copy'

    input:
    path deseq_results
    path transcript_mag_map
    val outdir_name

    output:
    path "deseq2_mag_mapped_results.tsv", emit: mag_mapped_results
    path "deseq2_significant_results.tsv", emit: significant_results
    path "expression_summary.tsv", emit: summary
    path "mag_expression_summary.tsv", emit: mag_summary
    path "expression_summary_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_summary_tables.R \
        --deseq-results ${deseq_results} \
        --mag-map ${transcript_mag_map} \
        --output-mag-mapped deseq2_mag_mapped_results.tsv \
        --output-significant deseq2_significant_results.tsv \
        --output-summary expression_summary.tsv \
        --output-mag-summary mag_expression_summary.tsv

    {
        echo '"EXPRESSION_SUMMARY":'
        echo "    R: \$(Rscript -e 'cat(as.character(getRversion()))')"
        echo "    dplyr: \$(Rscript -e 'cat(as.character(packageVersion(\"dplyr\")))')"
        echo "    readr: \$(Rscript -e 'cat(as.character(packageVersion(\"readr\")))')"
    } > expression_summary_versions.yml
    """

    stub:
    """
    cat > deseq2_mag_mapped_results.tsv <<'STUB'
gene_idbaseMeanlog2FoldChangelfcSEstatpvaluepadjqseqidtranscript_idmag_idsseqidpidentevaluebitscore
fake_MAG_protein_1450.3NANANANAfake_orffake_transcriptfake_MAGfake_protein991e-50200
STUB

    cat > deseq2_significant_results.tsv <<'STUB'
gene_idbaseMeanlog2FoldChangelfcSEstatpvaluepadj
STUB

    cat > expression_summary.tsv <<'STUB'
metricvalue
total_features1
features_with_padj0
significant_padj_lt_0.050
upregulated_log2fc_gt_00
downregulated_log2fc_lt_00
STUB

    cat > mag_expression_summary.tsv <<'STUB'
mag_idtotal_featuressignificant_featuresmean_log2FoldChange
fake_MAG100.3
STUB

    cat > expression_summary_versions.yml <<'STUB'
"EXPRESSION_SUMMARY":
    R: stub
    dplyr: stub
    readr: stub
STUB
    """
}
