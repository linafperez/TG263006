include { EXPRESSION_MULTIQC } from '../modules/local/expression_multiqc/main.nf'
include { SPADES_MULTIQC } from '../modules/local/spades_multiqc/main.nf'
include { TRANSDECODER_MULTIQC } from '../modules/local/transdecoder_multiqc/main.nf'
include { DIAMOND_MULTIQC } from '../modules/local/diamond_multiqc/main.nf'
include { INPUT_PREPARATION } from '../subworkflows/local/input_preparation/main.nf'

include { QUALITY_CONTROL } from '../subworkflows/core/quality_control/main.nf'
include { RRNA_AND_HOST_REMOVAL } from '../subworkflows/core/rrna_and_host_removal/main.nf'
include { ASSEMBLY } from '../subworkflows/core/assembly/main.nf'
include { TRANSCRIPT_PREDICTION } from '../subworkflows/core/transcript_prediction/main.nf'
include { MAG_ANNOTATION } from '../subworkflows/core/mag_annotation/main.nf'
include { QUANTIFICATION } from '../subworkflows/core/quantification/main.nf'
include { REPORTING_QC } from '../subworkflows/core/reporting_qc/main.nf'
include { VERSIONS_MULTIQC } from '../modules/local/versions_multiqc/main.nf'

include { RUN_EXPRESSION_ANALYSIS } from '../subworkflows/local/expression_analysis/main.nf'


workflow METATRANSCRIPTOMICS {

    take:
    samplesheet
    rrna_db_dir
    star_index
    mag_protein_dir
    condition_column
    case_label
    control_label
    expression_outdir_name
    expression_report_name

    main:

    /*
     * 1. Validar samplesheet
     */
    INPUT_PREPARATION(
        samplesheet
    )

    /*
     * 2. Crear canal de reads desde samplesheet validado
     */
    reads_ch = INPUT_PREPARATION.out.samplesheet
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                condition: row.condition,
                sra_id: row.sra_id
            ]

            tuple(
                meta,
                [
                    file(row.fastq_1),
                    file(row.fastq_2)
                ]
            )
        }

    /*
     * 3. Quality control:
     * FastQC raw → fastp → FastQC clean
     */

    QUALITY_CONTROL(
        reads_ch
    )

    /*
     * 4. Remoción rRNA y humana:
     * SortMeRNA → STAR
     */

    rrna_idx_dir = rrna_db_dir.resolve('idx')

    RRNA_AND_HOST_REMOVAL(
        QUALITY_CONTROL.out.trimmed_reads,
        rrna_db_dir,
        rrna_idx_dir,
        star_index
    )

    /*
     * 5. Assembly:
     * rnaSPAdes → QUAST/rnaQUAST
     */
    ASSEMBLY(
        RRNA_AND_HOST_REMOVAL.out.non_human_reads
    )

    SPADES_MULTIQC(
        ASSEMBLY.out.transcripts,
        "coassembly"
    )

    /*
     * 6. Transcript prediction:
     * TransDecoder → filtro de ORFs completos
     */
    TRANSCRIPT_PREDICTION(
        ASSEMBLY.out.transcripts
    )

    TRANSDECODER_MULTIQC(
        TRANSCRIPT_PREDICTION.out.pep,
        TRANSCRIPT_PREDICTION.out.complete_pep,
        "coassembly"
    )

    /*
     * 7. MAG annotation:
     * MAG protein DB → DIAMOND → best hits → transcript-MAG map/GTF
     */
    MAG_ANNOTATION(
        TRANSCRIPT_PREDICTION.out.complete_pep,
        TRANSCRIPT_PREDICTION.out.complete_gff3,
        mag_protein_dir
    )

    DIAMOND_MULTIQC(
        MAG_ANNOTATION.out.diamond_hits,
        MAG_ANNOTATION.out.diamond_best_hits,
        "coassembly"
    )

    /*
     * 8. Quantification:
     * Salmon index → Salmon quant
     */
    QUANTIFICATION(
        ASSEMBLY.out.transcripts,
        RRNA_AND_HOST_REMOVAL.out.non_human_reads
    )

    /*
     * 9. Expression analysis:
     * tximport → DESeq2 → plots/tables
     */
    RUN_EXPRESSION_ANALYSIS(
        QUANTIFICATION.out.quant_dirs,
        MAG_ANNOTATION.out.transcript_mag_map,
        INPUT_PREPARATION.out.samplesheet,
        condition_column,
        case_label,
        control_label,
        expression_outdir_name
    )

    EXPRESSION_MULTIQC(
        RUN_EXPRESSION_ANALYSIS.out.summary,
        RUN_EXPRESSION_ANALYSIS.out.mag_summary,
        expression_report_name
    )

    /*
     * 10. MultiQC final
     */
    version_files_ch = Channel
        .empty()
        .mix(QUALITY_CONTROL.out.fastqc_raw_versions.map { meta, file -> file })
        .mix(QUALITY_CONTROL.out.fastp_versions.map { meta, file -> file })
        .mix(QUALITY_CONTROL.out.fastqc_clean_versions.map { meta, file -> file })
        .mix(RRNA_AND_HOST_REMOVAL.out.sortmerna_versions.map { meta, file -> file })
        .mix(RRNA_AND_HOST_REMOVAL.out.star_versions.map { meta, file -> file })
        .mix(ASSEMBLY.out.spades_versions)
        .mix(ASSEMBLY.out.quast_versions)
        .mix(TRANSCRIPT_PREDICTION.out.transdecoder_versions)
        .mix(TRANSCRIPT_PREDICTION.out.filter_versions)
        .mix(MAG_ANNOTATION.out.prepare_versions)
        .mix(MAG_ANNOTATION.out.diamond_versions)
        .mix(MAG_ANNOTATION.out.filter_versions)
        .mix(MAG_ANNOTATION.out.gtf_versions)
        .mix(QUANTIFICATION.out.salmon_index_versions)
        .mix(QUANTIFICATION.out.salmon_quant_versions.map { meta, file -> file })
        .mix(RUN_EXPRESSION_ANALYSIS.out.versions)
        .collect()

    VERSIONS_MULTIQC(version_files_ch)

    multiqc_input_ch = Channel
        .empty()
        .mix(QUALITY_CONTROL.out.fastqc_raw_zip.map { meta, reports -> reports }.flatten())
        .mix(QUALITY_CONTROL.out.fastqc_clean_zip.map { meta, reports -> reports }.flatten())
        .mix(QUALITY_CONTROL.out.fastp_json.map { meta, file -> file })
        .mix(RRNA_AND_HOST_REMOVAL.out.sortmerna_log.map { meta, file -> file })
        .mix(RRNA_AND_HOST_REMOVAL.out.star_log.map { meta, file -> file })
        .mix(SPADES_MULTIQC.out.yaml)
        .mix(ASSEMBLY.out.rnaquast_multiqc)
        .mix(TRANSDECODER_MULTIQC.out.yaml)
        .mix(DIAMOND_MULTIQC.out.yaml)
        .mix(QUANTIFICATION.out.quant_dirs.flatten())
        .mix(EXPRESSION_MULTIQC.out.yaml)
        .mix(RUN_EXPRESSION_ANALYSIS.out.pca_mqc)
        .mix(RUN_EXPRESSION_ANALYSIS.out.distance_heatmap_mqc)
        .mix(RUN_EXPRESSION_ANALYSIS.out.dispersion_mqc)
        .mix(RUN_EXPRESSION_ANALYSIS.out.deg_heatmap_mqc)
        .mix(RUN_EXPRESSION_ANALYSIS.out.ma_plot_mqc)
        .mix(RUN_EXPRESSION_ANALYSIS.out.volcano_mqc)
        .mix(VERSIONS_MULTIQC.out.yaml)
        .collect()

    REPORTING_QC(
        multiqc_input_ch
    )

    emit:
    validated_samplesheet = INPUT_PREPARATION.out.samplesheet

    trimmed_reads   = QUALITY_CONTROL.out.trimmed_reads
    non_human_reads = RRNA_AND_HOST_REMOVAL.out.non_human_reads
    transcripts     = ASSEMBLY.out.transcripts

    predicted_pep  = TRANSCRIPT_PREDICTION.out.pep
    predicted_cds  = TRANSCRIPT_PREDICTION.out.cds
    predicted_gff3 = TRANSCRIPT_PREDICTION.out.gff3

    transcript_mag_map = MAG_ANNOTATION.out.transcript_mag_map
    transcript_mag_gtf = MAG_ANNOTATION.out.transcript_mag_gtf
    diamond_best_hits  = MAG_ANNOTATION.out.diamond_best_hits

    salmon_quant_dirs = QUANTIFICATION.out.quant_dirs

    expression_summary   = RUN_EXPRESSION_ANALYSIS.out.summary
    deseq2_results       = RUN_EXPRESSION_ANALYSIS.out.deseq2_results
    mag_mapped_results   = RUN_EXPRESSION_ANALYSIS.out.mag_mapped_results
    significant_results  = RUN_EXPRESSION_ANALYSIS.out.significant_results
    mag_summary          = RUN_EXPRESSION_ANALYSIS.out.mag_summary

    multiqc_report = REPORTING_QC.out.html
}
