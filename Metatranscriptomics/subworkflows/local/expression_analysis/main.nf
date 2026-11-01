include { TXIMPORT } from '../../../modules/local/expression_analysis/tximport.nf'
include { DESEQ2_PD_VS_HC } from '../../../modules/local/expression_analysis/deseq2.nf'
include { EXPRESSION_SUMMARY } from '../../../modules/local/expression_analysis/summary.nf'


workflow RUN_EXPRESSION_ANALYSIS {

    take:
    salmon_quant_dirs_ch
    transcript_mag_map
    samplesheet_ch
    condition_column
    case_label
    control_label
    outdir_name

    main:

    TXIMPORT(
        salmon_quant_dirs_ch,
        transcript_mag_map,
        outdir_name
    )

    DESEQ2_PD_VS_HC(
        TXIMPORT.out.tximport_rds,
        samplesheet_ch,
        condition_column,
        case_label,
        control_label,
        outdir_name
    )

    EXPRESSION_SUMMARY(
        DESEQ2_PD_VS_HC.out.deseq2_results,
        transcript_mag_map,
        outdir_name
    )

    expression_versions_ch = Channel
        .empty()
        .mix(TXIMPORT.out.versions)
        .mix(DESEQ2_PD_VS_HC.out.versions)
        .mix(EXPRESSION_SUMMARY.out.versions)

    emit:

    tximport_rds = TXIMPORT.out.tximport_rds
    counts       = TXIMPORT.out.counts
    tpm          = TXIMPORT.out.tpm
    length       = TXIMPORT.out.length

    deseq2_rds       = DESEQ2_PD_VS_HC.out.deseq2_rds
    deseq2_results   = DESEQ2_PD_VS_HC.out.deseq2_results
    qc_summary       = DESEQ2_PD_VS_HC.out.qc_summary
    vst_counts       = DESEQ2_PD_VS_HC.out.vst_counts

    pca              = DESEQ2_PD_VS_HC.out.pca
    distance_heatmap = DESEQ2_PD_VS_HC.out.distance_heatmap
    dispersion       = DESEQ2_PD_VS_HC.out.dispersion
    deg_heatmap      = DESEQ2_PD_VS_HC.out.deg_heatmap
    volcano          = DESEQ2_PD_VS_HC.out.volcano
    ma_plot          = DESEQ2_PD_VS_HC.out.ma_plot

    pca_mqc              = DESEQ2_PD_VS_HC.out.pca_mqc
    distance_heatmap_mqc = DESEQ2_PD_VS_HC.out.distance_heatmap_mqc
    dispersion_mqc       = DESEQ2_PD_VS_HC.out.dispersion_mqc
    deg_heatmap_mqc      = DESEQ2_PD_VS_HC.out.deg_heatmap_mqc
    ma_plot_mqc          = DESEQ2_PD_VS_HC.out.ma_plot_mqc
    volcano_mqc          = DESEQ2_PD_VS_HC.out.volcano_mqc

    mag_mapped_results  = EXPRESSION_SUMMARY.out.mag_mapped_results
    significant_results = EXPRESSION_SUMMARY.out.significant_results
    summary             = EXPRESSION_SUMMARY.out.summary
    mag_summary         = EXPRESSION_SUMMARY.out.mag_summary

    versions = expression_versions_ch
}
