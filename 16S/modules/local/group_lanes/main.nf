workflow GROUP_LANES {

    take:
    parsed

    main:

    grouped = parsed
        .groupTuple()

   
    emit:

    grouped

}
