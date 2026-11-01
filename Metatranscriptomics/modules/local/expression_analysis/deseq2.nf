process DESEQ2_PD_VS_HC {

    tag "${case_label}_vs_${control_label}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    publishDir {
        "${params.outdir}/expression_analysis/deseq2/${outdir_name}"
    }, mode: 'copy'

    input:
    path tximport_rds
    path samplesheet
    val condition_column
    val case_label
    val control_label
    val outdir_name

    output:
    path "deseq2_dataset.rds", emit: deseq2_rds
    path "deseq2_pd_vs_hc.tsv", emit: deseq2_results
    path "deseq2_pd_vs_hc_unshrunken.tsv", emit: deseq2_unshrunken
    path "deseq2_qc_summary.tsv", emit: qc_summary
    path "vst_counts.tsv", emit: vst_counts

    path "pca_plot.pdf", emit: pca
    path "sample_distance_heatmap.pdf", emit: distance_heatmap
    path "dispersion_plot.pdf", emit: dispersion
    path "deg_heatmap.pdf", emit: deg_heatmap
    path "ma_plot.pdf", emit: ma_plot
    path "volcano_plot.pdf", emit: volcano

    path "pca_plot_mqc.png", emit: pca_mqc
    path "sample_distance_heatmap_mqc.png", emit: distance_heatmap_mqc
    path "dispersion_plot_mqc.png", emit: dispersion_mqc
    path "deg_heatmap_mqc.png", emit: deg_heatmap_mqc
    path "ma_plot_mqc.png", emit: ma_plot_mqc
    path "volcano_plot_mqc.png", emit: volcano_mqc

    path "deseq2_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_deseq2_pd_vs_hc.R \
        --tximport-rds ${tximport_rds} \
        --samplesheet ${samplesheet} \
        --sample-column sample_id \
        --condition-column ${condition_column} \
        --case ${case_label} \
        --control ${control_label} \
        --min-sample-counts 100000 \
        --min-gene-counts 10 \
        --output-results deseq2_pd_vs_hc.tsv \
        --output-unshrunken deseq2_pd_vs_hc_unshrunken.tsv \
        --output-rds deseq2_dataset.rds \
        --output-vst vst_counts.tsv \
        --output-pca pca_plot.pdf \
        --output-distance-heatmap sample_distance_heatmap.pdf \
        --output-dispersion dispersion_plot.pdf \
        --output-ma ma_plot.pdf \
        --output-volcano volcano_plot.pdf \
        --output-deg-heatmap deg_heatmap.pdf \
        --output-qc-summary deseq2_qc_summary.tsv

    for plot in \
        pca_plot \
        sample_distance_heatmap \
        dispersion_plot \
        deg_heatmap \
        ma_plot \
        volcano_plot
    do
        pdftoppm \
            -png \
            -singlefile \
            -r 150 \
            "\${plot}.pdf" \
            "\${plot}_mqc"
    done

    {
        echo '"DESEQ2_PD_VS_HC":'
        echo "    R: \$(Rscript -e 'cat(as.character(getRversion()))')"
        echo "    DESeq2: \$(Rscript -e 'cat(as.character(packageVersion(\"DESeq2\")))')"
        echo "    apeglm: \$(Rscript -e 'cat(as.character(packageVersion(\"apeglm\")))')"
        echo "    PCAtools: \$(Rscript -e 'cat(as.character(packageVersion(\"PCAtools\")))')"
        echo "    ggplot2: \$(Rscript -e 'cat(as.character(packageVersion(\"ggplot2\")))')"
        echo "    pheatmap: \$(Rscript -e 'cat(as.character(packageVersion(\"pheatmap\")))')"
    } > deseq2_versions.yml
    """

    stub:
    """
    touch deseq2_dataset.rds

    cat > deseq2_pd_vs_hc_unshrunken.tsv <<'STUB'
gene_idbaseMeanlog2FoldChangelfcSEstatpvaluepadj
fake_MAG_protein_1450.3NANANANA
STUB

    cat > deseq2_pd_vs_hc.tsv <<'STUB'
gene_idbaseMeanlog2FoldChangelfcSEstatpvaluepadj
fake_MAG_protein_1450.3NANANANA
STUB

    cat > deseq2_qc_summary.tsv <<'STUB'
sample_idtotal_countskeep_samplecondition
PD_01150000TRUEPD
HC_01140000TRUEHC
STUB

    cat > vst_counts.tsv <<'STUB'
gene_idPD_01HC_01
fake_MAG_protein_15.65.3
STUB

    touch \
        pca_plot.pdf \
        sample_distance_heatmap.pdf \
        dispersion_plot.pdf \
        deg_heatmap.pdf \
        ma_plot.pdf \
        volcano_plot.pdf

    MQC_PNG_BASE64="iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="

    for f in \
        pca_plot_mqc.png \
        sample_distance_heatmap_mqc.png \
        dispersion_plot_mqc.png \
        deg_heatmap_mqc.png \
        ma_plot_mqc.png \
        volcano_plot_mqc.png
    do
        printf "%s" "\$MQC_PNG_BASE64" | base64 -d > "\$f"
    done

    cat > deseq2_versions.yml <<'STUB'
"DESEQ2_PD_VS_HC":
    R: stub
    DESeq2: stub
    apeglm: stub
    PCAtools: stub
    ggplot2: stub
    pheatmap: stub
STUB
    """
}
