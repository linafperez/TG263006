process TXIMPORT {

    tag "${outdir_name}"
    label 'process_high'

    conda "${moduleDir}/environment.yaml"

    publishDir {
        "${params.outdir}/expression_analysis/tximport/${outdir_name}"
    }, mode: 'copy'

    input:
    path salmon_quant_dirs
    path transcript_mag_map
    val outdir_name

    output:
    path "salmon_dirs.txt", emit: salmon_dirs
    path "tximport_results.rds", emit: tximport_rds
    path "tximport_counts.tsv", emit: counts
    path "tximport_tpm.tsv", emit: tpm
    path "tximport_length.tsv", emit: length
    path "tximport_versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf "%s\\n" ${salmon_quant_dirs} > salmon_dirs.txt

    run_tximport_salmon.R \
        --salmon-dirs-file salmon_dirs.txt \
        --tx2gene ${transcript_mag_map} \
        --output-rds tximport_results.rds \
        --counts-output tximport_counts.tsv \
        --abundance-output tximport_tpm.tsv \
        --length-output tximport_length.tsv

    {
        echo '"TXIMPORT":'
        echo "    R: \$(Rscript -e 'cat(as.character(getRversion()))')"
        echo "    tximport: \$(Rscript -e 'cat(as.character(packageVersion(\"tximport\")))')"
    } > tximport_versions.yml
    """

    stub:
    """
    printf "%s\\n" PD_01 HC_01 > salmon_dirs.txt

    touch tximport_results.rds

    cat > tximport_counts.tsv <<'STUB_COUNTS'
gene_idPD_01HC_01
fake_MAG_protein_15040
STUB_COUNTS

    cat > tximport_tpm.tsv <<'STUB_TPM'
gene_idPD_01HC_01
fake_MAG_protein_110.58.1
STUB_TPM

    cat > tximport_length.tsv <<'STUB_LENGTH'
gene_idPD_01HC_01
fake_MAG_protein_1800800
STUB_LENGTH

    cat > tximport_versions.yml <<'STUB_VERSIONS'
"TXIMPORT":
    R: stub
    tximport: stub
STUB_VERSIONS
    """
}
