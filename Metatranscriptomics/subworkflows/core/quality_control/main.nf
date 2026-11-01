include { FASTQC as FASTQC_RAW } from '../../../modules/core/fastqc/main.nf'
include { FASTQC as FASTQC_CLEAN } from '../../../modules/core/fastqc/main.nf'
include { FASTP } from '../../../modules/core/fastp/main.nf'

workflow QUALITY_CONTROL {

    take:
    reads_ch

    main:

    FASTQC_RAW(
        reads_ch,
        "fastqc_raw"
    )

    FASTP(
        reads_ch,
        "fastp_trimmed"
    )

    trimmed_reads_ch = FASTP.out.reads
        .map { meta, read1, read2 ->
            tuple(
                meta,
                [
                    read1,
                    read2
                ]
            )
        }

    FASTQC_CLEAN(
        trimmed_reads_ch,
        "fastqc_clean"
    )

    emit:
    fastqc_raw_versions   = FASTQC_RAW.out.versions
    fastp_versions        = FASTP.out.versions
    fastqc_clean_versions = FASTQC_CLEAN.out.versions

    trimmed_reads = trimmed_reads_ch

    fastqc_raw_html   = FASTQC_RAW.out.html
    fastqc_clean_html = FASTQC_CLEAN.out.html
    
    fastqc_raw_zip    = FASTQC_RAW.out.zip
    fastqc_clean_zip  = FASTQC_CLEAN.out.zip

    fastp_html        = FASTP.out.html
    fastp_json        = FASTP.out.json
}
