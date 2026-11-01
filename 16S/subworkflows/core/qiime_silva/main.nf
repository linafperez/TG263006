nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QIIME2 SILVA SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    SILVA can be obtained in two ways:

    1. Automatic mode
       If no local SILVA artifacts are provided:
           SILVA_DOWNLOAD
                ↓
           SILVA_CURATION

    2. Local mode
       If both params.silva_sequences and params.silva_taxonomy
       are provided:
           local SILVA QIIME2 artifacts
                ↓
           SILVA_CURATION

    Outputs:
        clean_sequences
        dereplicated_sequences
        dereplicated_taxonomy

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include {
    SILVA_DOWNLOAD
} from '../../../modules/core/qiime_silva/download'


include {
    SILVA_CURATION
} from '../../../modules/core/qiime_silva/curate'


workflow QIIME_SILVA {


    main:

    if (params.silva_sequences && params.silva_taxonomy) {

        log.info "SILVA source: local QIIME2 artifacts"

        silva_sequences_ch = Channel.fromPath(
            params.silva_sequences,
            checkIfExists: true
        )

        silva_taxonomy_ch = Channel.fromPath(
            params.silva_taxonomy,
            checkIfExists: true
        )

    } else {

        log.info "SILVA source: automatic download"

        SILVA_DOWNLOAD()

        silva_sequences_ch =
            SILVA_DOWNLOAD.out.sequences

        silva_taxonomy_ch =
            SILVA_DOWNLOAD.out.taxonomy
    }


    SILVA_CURATION(
        silva_sequences_ch,
        silva_taxonomy_ch
    )


    emit:

    clean_sequences =
        SILVA_CURATION.out.clean_sequences

    dereplicated_sequences =
        SILVA_CURATION.out.dereplicated_sequences

    dereplicated_taxonomy =
        SILVA_CURATION.out.dereplicated_taxonomy
}
