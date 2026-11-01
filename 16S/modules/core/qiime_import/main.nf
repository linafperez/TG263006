nextflow.enable.dsl=2

include {
    IMPORT
} from './import'

workflow QIIME_IMPORT {

    take:
    samples

    main:

    IMPORT(samples)

    emit:

    manifest = IMPORT.out.manifest
    demux    = IMPORT.out.demux
}
