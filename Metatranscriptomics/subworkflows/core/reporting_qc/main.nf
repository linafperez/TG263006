include { MULTIQC } from '../../../modules/core/multiqc/main.nf'

workflow REPORTING_QC {

    take:
    multiqc_input_ch

    main:

    MULTIQC(
        multiqc_input_ch
    )

    emit:
    html     = MULTIQC.out.html
    data     = MULTIQC.out.data
    versions = MULTIQC.out.versions
}
