workflow PARSE_SAMPLESHEET {

    take:
    samplesheet

    main:

    parsed = samplesheet
        .splitCsv(header: true)
        .map { row ->

            tuple(
                row.sample,
                file(row.fastq_1),
                file(row.fastq_2)
            )

        }

    emit:

    parsed
}
