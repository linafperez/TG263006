process MERGE_LANES {

    tag "$sample"

    publishDir "${params.outdir}/Preprocessing/merged_fastq",
        mode: "copy"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample),
          path("${sample}_R1.fastq.gz"),
          path("${sample}_R2.fastq.gz"),
          emit: merged

    script:

    def r1 = reads1.collect { it.getName() }.join(" ")
    def r2 = reads2.collect { it.getName() }.join(" ")

    def singleLane = (reads1.size() == 1)

    """
    set -euo pipefail

    if ${singleLane}; then

        ln -s ${reads1[0].getName()} ${sample}_R1.fastq.gz
        ln -s ${reads2[0].getName()} ${sample}_R2.fastq.gz

    else

        cat ${r1} > ${sample}_R1.fastq.gz
        cat ${r2} > ${sample}_R2.fastq.gz

    fi
    """
}
