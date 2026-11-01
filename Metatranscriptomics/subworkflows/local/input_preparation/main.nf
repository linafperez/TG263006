include { CHECK_SAMPLESHEET } from '../../../modules/local/check_samplesheet/main.nf'

workflow INPUT_PREPARATION {

    take:
    samplesheet

    main:

    CHECK_SAMPLESHEET(
        samplesheet
    )

    emit:
    samplesheet            = CHECK_SAMPLESHEET.out.samplesheet
    check_samplesheet_log  = CHECK_SAMPLESHEET.out.log
    versions               = CHECK_SAMPLESHEET.out.versions
}
