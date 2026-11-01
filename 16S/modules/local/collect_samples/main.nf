workflow COLLECT_SAMPLES {

    take:
    samples

    main:

    collected = samples.collect()

    emit:

    collected
}
