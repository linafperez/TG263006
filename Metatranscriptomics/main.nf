nextflow.enable.dsl = 2

/*
 * Main pipeline:
 * Metatranscriptomics PD vs HC
 */

include { METATRANSCRIPTOMICS } from './workflows/metatranscriptomics.nf'


workflow {

    /*
     * Input files and reference directories.
     *
     * Default parameter values are defined in:
     *   conf/base.config
     *   conf/genomes.config
     */

    samplesheet_ch = Channel.value(
        file(
            params.input,
            checkIfExists: true
        )
    )

    rrna_db_dir_ch = Channel.value(
        file(
            params.rrna_db_dir,
            checkIfExists: true
        )
    )

    star_index_ch = Channel.value(
        file(
            params.star_index,
            checkIfExists: true
        )
    )

    mag_protein_dir_ch = Channel.value(
        file(
            params.mag_protein_dir,
            checkIfExists: true
        )
    )

    condition_column_ch = Channel.value(params.condition_column)
    case_label_ch = Channel.value(params.case_label)
    control_label_ch = Channel.value(params.control_label)
    expression_outdir_name_ch = Channel.value(params.expression_outdir_name)
    expression_report_name_ch = Channel.value(params.expression_report_name)

    METATRANSCRIPTOMICS(
        samplesheet_ch,
        rrna_db_dir_ch,
        star_index_ch,
        mag_protein_dir_ch,
        condition_column_ch,
        case_label_ch,
        control_label_ch,
        expression_outdir_name_ch,
        expression_report_name_ch
    )
}
