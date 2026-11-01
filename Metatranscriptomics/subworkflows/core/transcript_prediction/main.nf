include { TRANSDECODER } from '../../../modules/core/transdecoder/main.nf'
include { FILTER_TRANSDECODER_COMPLETE_ISOFORMS } from '../../../modules/local/filter_transdecoder_complete_isoforms/main.nf'

workflow TRANSCRIPT_PREDICTION {

    take:
    transcripts

    main:

    TRANSDECODER(
        transcripts,
        "coassembly",
        "coassembly"
    )

    FILTER_TRANSDECODER_COMPLETE_ISOFORMS(
        TRANSDECODER.out.pep,
        TRANSDECODER.out.cds,
        TRANSDECODER.out.gff3,
        "coassembly",
        "coassembly"
    )

    emit:

    /*
     * Raw TransDecoder outputs.
     *
     * These are used by MAG annotation / DIAMOND to match
     * the workflow used in the reference scripts.
     */
    pep  = TRANSDECODER.out.pep
    cds  = TRANSDECODER.out.cds
    gff3 = TRANSDECODER.out.gff3

    /*
     * Complete-CDS subset retained as an explicit output.
     */
    complete_pep  = FILTER_TRANSDECODER_COMPLETE_ISOFORMS.out.pep
    complete_cds  = FILTER_TRANSDECODER_COMPLETE_ISOFORMS.out.cds
    complete_gff3 = FILTER_TRANSDECODER_COMPLETE_ISOFORMS.out.gff3

    transdecoder_versions = TRANSDECODER.out.versions
    filter_log            = FILTER_TRANSDECODER_COMPLETE_ISOFORMS.out.log
    filter_versions       = FILTER_TRANSDECODER_COMPLETE_ISOFORMS.out.versions
}
