nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE SAMPLES SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Input:
        samplesheet.csv

    Steps:
        SAMPLESHEET_CHECK
        PARSE_SAMPLESHEET
        GROUP_LANES
        MERGE_LANES
        FASTQC

    Outputs:
        validated_samplesheet
        merged paired-end FASTQ
        FastQC HTML reports
        FastQC ZIP reports

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    SAMPLESHEET_CHECK
} from '../../../modules/local/samplesheet_check'


include {
    PARSE_SAMPLESHEET
} from '../../../modules/local/parse_samplesheet'


include {
    GROUP_LANES
} from '../../../modules/local/group_lanes'


include {
    MERGE_LANES
} from '../../../modules/local/merge_lanes'


include {
    FASTQC
} from '../../../modules/core/fastqc'


workflow PREPARE_SAMPLES {

    take:
    samplesheet


    main:

    /*
     * Validate samplesheet structure and FASTQ paths.
     */

    SAMPLESHEET_CHECK(
        samplesheet
    )


    /*
     * Convert validated samplesheet rows into:
     *
     * tuple(
     *     sample,
     *     fastq_1,
     *     fastq_2
     * )
     */

    PARSE_SAMPLESHEET(
        SAMPLESHEET_CHECK.out.validated_samplesheet
    )


    /*
     * Group multiple sequencing lanes belonging
     * to the same biological sample.
     */

    GROUP_LANES(
        PARSE_SAMPLESHEET.out.parsed
    )


    /*
     * Merge lanes when necessary.
     *
     * Single-lane samples are linked without
     * recompressing the FASTQ files.
     */

    MERGE_LANES(
        GROUP_LANES.out.grouped
    )


    /*
     * Initial FASTQ quality control.
     */

    FASTQC(
        MERGE_LANES.out.merged
    )


    emit:

    validated_samplesheet =
        SAMPLESHEET_CHECK.out.validated_samplesheet

    merged =
        MERGE_LANES.out.merged

    fastqc_html =
        FASTQC.out.html

    fastqc_zip =
        FASTQC.out.zip
}
