process PREPARE_MAG_PROTEIN_DATABASE {

    tag "MAGs_protein_database"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.10--h43eeafb_0' :
        'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0' }"

    publishDir { "${params.outdir}/mag_annotation/prepare_mag_protein_database" }, mode: 'copy'

    input:
    path protein_dir
    val outdir_name

    output:
    path "MAGs_protein.faa", emit: fasta
    path "MAGs_protein.dmnd", emit: db
    path "prepare_mag_protein_database.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    find -L ${protein_dir} \
        -type f \
        -name '*.proteins.faa' \
        | sort > mag_protein_files.txt

    TOTAL_FILES=\$(wc -l < mag_protein_files.txt)

    if [ "\${TOTAL_FILES}" -eq 0 ]; then
        echo "ERROR: No *.proteins.faa files found in ${protein_dir}" >&2
        exit 1
    fi

    : > MAGs_protein.faa

    while IFS= read -r faa; do
        cat "\${faa}" >> MAGs_protein.faa
    done < mag_protein_files.txt

    TOTAL_SEQUENCES=\$(grep -c '^>' MAGs_protein.faa || true)

    if [ "\${TOTAL_SEQUENCES}" -eq 0 ]; then
        echo "ERROR: No protein sequences found" >&2
        exit 1
    fi

    {
        echo "Found \${TOTAL_FILES} MAG protein FASTA files"
        echo "Wrote \${TOTAL_SEQUENCES} protein sequences"
        echo "Output FASTA: MAGs_protein.faa"
    } > prepare_mag_protein_database.log

    diamond makedb \
        --in MAGs_protein.faa \
        --db MAGs_protein \
        >> prepare_mag_protein_database.log 2>&1

    printf '"${task.process}":\n    diamond: %s\n' \
        "\$(diamond version | head -n 1 | sed 's/diamond version //')" \
        > versions.yml
    """

    stub:
    """
    printf '>fake_MAG_protein_1\nMRTLASLLALAAVAA\n' > MAGs_protein.faa

    touch MAGs_protein.dmnd
    touch prepare_mag_protein_database.log

    printf '"${task.process}":\n    diamond: 2.1.10\n' > versions.yml
    """
}
