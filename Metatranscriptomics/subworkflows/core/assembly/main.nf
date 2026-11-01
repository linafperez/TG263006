include { SPADES_RNA } from '../../../modules/core/spades_rna/main.nf'
include { QUAST } from '../../../modules/core/quast/main.nf'
include { RNAQUAST_MULTIQC } from '../../../modules/local/rnaquast_multiqc/main.nf'

workflow ASSEMBLY {

    take:
    non_human_reads_ch

    main:

    spades_reads_ch = non_human_reads_ch
        .map { meta, read1, read2 ->
            [
                read1,
                read2
            ]
        }
        .flatten()
        .collect()

    SPADES_RNA(
        spades_reads_ch,
        "coassembly",
        "coassembly"
    )

    QUAST(
        SPADES_RNA.out.transcripts,
        "coassembly",
        "coassembly"
    )

    RNAQUAST_MULTIQC(
        QUAST.out.report_dir,
        "coassembly",
        "coassembly"
    )

    emit:
    transcripts = SPADES_RNA.out.transcripts

    spades_log   = SPADES_RNA.out.log
    quast_report = QUAST.out.report_dir
    rnaquast_multiqc = RNAQUAST_MULTIQC.out.yaml

    spades_versions = SPADES_RNA.out.versions
    quast_versions  = QUAST.out.versions
}
