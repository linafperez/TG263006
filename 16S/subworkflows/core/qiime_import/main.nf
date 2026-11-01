nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QIIME2 IMPORT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Input:
        merged paired-end FASTQ samples

    Steps:
        IMPORT
        DEMUX_SUMMARIZE

    Outputs:
        manifest
        demux
        demux_summary

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    IMPORT
} from '../../../modules/core/qiime_import/import'


include {
    DEMUX_SUMMARIZE
} from '../../../modules/core/qiime_import/demux'


workflow QIIME_IMPORT {

    take:
    samples


    main:

    IMPORT(
        samples
    )


    DEMUX_SUMMARIZE(
        IMPORT.out.demux
    )


    emit:

    manifest =
        IMPORT.out.manifest

    demux =
        IMPORT.out.demux

    demux_summary =
        DEMUX_SUMMARIZE.out.summary
}
