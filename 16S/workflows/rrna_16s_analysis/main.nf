nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    16S rRNA ANALYSIS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    General-purpose workflow for paired-end 16S rRNA amplicon analysis.

    Architecture:

        LOCAL
        PREPARE_SAMPLES
              ↓
        CORE
        QIIME_IMPORT
              ↓
        QIIME_DENOISING
              ↓
        QIIME_SILVA
              ↓
        QIIME_TAXONOMY
              ↓
        LOCAL
        ABUNDANCE
          ├── ALPHA_DIVERSITY
          └── BETA_DIVERSITY

    Sample metadata is optional and validated once before being passed
    to downstream R subworkflows.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
===============================================================================
    CORE SUBWORKFLOWS
===============================================================================
*/


include {
    QIIME_IMPORT
} from '../../subworkflows/core/qiime_import'


include {
    QIIME_DENOISING
} from '../../subworkflows/core/qiime_denoising'


include {
    QIIME_SILVA
} from '../../subworkflows/core/qiime_silva'


include {
    QIIME_TAXONOMY
} from '../../subworkflows/core/qiime_taxonomy'


/*
===============================================================================
    LOCAL SUBWORKFLOWS
===============================================================================
*/


include {
    PREPARE_SAMPLES
} from '../../subworkflows/local/prepare_samples'


include {
    ABUNDANCE
} from '../../subworkflows/local/abundance'


include {
    ALPHA_DIVERSITY
} from '../../subworkflows/local/alpha_diversity'


include {
    BETA_DIVERSITY
} from '../../subworkflows/local/beta_diversity'


/*
===============================================================================
    SHARED UTILITIES
===============================================================================
*/


include {
    METADATA_CHECK
} from '../../modules/local/metadata_check'


include {
    MULTIQC
} from '../../modules/core/multiqc'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow RRNA_16S_ANALYSIS {


    /*
     * ========================================================================
     * INPUT CHANNELS
     * ========================================================================
     */


    samplesheet_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )


    primers_ch = Channel.fromPath(
        params.primers,
        checkIfExists: true
    )


    /*
     * Metadata is optional.
     *
     * When no metadata is provided, an empty channel is created.
     * METADATA_CHECK will therefore not execute and its output channel
     * remains empty.
     */

    metadata_input_ch = params.metadata
        ? Channel.fromPath(
            params.metadata,
            checkIfExists: true
        )
        : Channel.empty()


    /*
     * ========================================================================
     * 1. LOCAL — PREPARE SAMPLES
     * ========================================================================
     *
     * samplesheet validation
     * parsing
     * lane grouping
     * lane merging
     * FastQC
     */


    PREPARE_SAMPLES(
        samplesheet_ch
    )


    /*
     * ========================================================================
     * OPTIONAL METADATA VALIDATION
     * ========================================================================
     *
     * Metadata is validated once and the same validated file is distributed
     * to downstream analyses.
     *
     * If metadata_input_ch is empty, this process does not execute.
     */


    METADATA_CHECK(
        metadata_input_ch,
        PREPARE_SAMPLES.out.validated_samplesheet
    )


    /*
     * ========================================================================
     * 2. CORE — QIIME IMPORT
     * ========================================================================
     *
     * merged FASTQ
     *      ↓
     * QIIME2 import
     *      ↓
     * demux summary
     *
     * IMPORT requires the complete collection of samples to construct
     * the manifest, therefore collect(flat:false) is intentional.
     */


    QIIME_IMPORT(
        PREPARE_SAMPLES.out.merged.collect(
            flat: false
        )
    )


    /*
     * ========================================================================
     * 3. CORE — QIIME DENOISING
     * ========================================================================
     *
     * demux.qza
     *      ↓
     * Cutadapt
     *      ↓
     * post-Cutadapt quality summary
     *      ↓
     * DADA2
     *      ↓
     * denoising summaries
     */


    QIIME_DENOISING(
        QIIME_IMPORT.out.demux,
        primers_ch
    )


    /*
     * ========================================================================
     * 4. CORE — QIIME SILVA
     * ========================================================================
     *
     * Current implementation:
     *
     * SILVA download
     *      ↓
     * curation
     *      ↓
     * dereplication
     *
     * Support for a user-provided local SILVA database will be added
     * separately.
     */


    QIIME_SILVA()


    /*
     * ========================================================================
     * 5. CORE — QIIME TAXONOMY
     * ========================================================================
     *
     * SILVA dereplicated sequences
     *      ↓
     * amplicon extraction
     *      ↓
     * classifier training
     *
     * representative ASVs
     *      ↓
     * taxonomic classification
     *      ↓
     * taxonomic filtering
     *      ↓
     * export for R
     */


    QIIME_TAXONOMY(
        QIIME_SILVA.out.dereplicated_sequences,
        QIIME_SILVA.out.dereplicated_taxonomy,
        QIIME_DENOISING.out.table,
        QIIME_DENOISING.out.rep_seqs,
        primers_ch
    )


    /*
     * ========================================================================
     * 6. LOCAL — ABUNDANCE
     * ========================================================================
     *
     * counts_raw.tsv
     *      ↓
     * prevalence filtering
     *      ↓
     * filtered_counts.tsv
     *      ↓
     * relative abundance
     *
     * relative abundance + taxonomy
     *      ↓
     * taxonomic composition
     *
     * Grouped composition is optional when metadata and
     * composition_group_column are provided.
     */


    ABUNDANCE(
        QIIME_TAXONOMY.out.counts_raw,
        QIIME_TAXONOMY.out.taxonomy_table,
        METADATA_CHECK.out.validated_metadata
    )


    /*
     * ========================================================================
     * 7. LOCAL — ALPHA DIVERSITY
     * ========================================================================
     *
     * IMPORTANT:
     *
     * Alpha diversity uses taxonomically filtered raw ASV counts
     * BEFORE prevalence filtering.
     */


    ALPHA_DIVERSITY(
        ABUNDANCE.out.counts_raw,
        METADATA_CHECK.out.validated_metadata
    )


    /*
     * ========================================================================
     * 8. LOCAL — BETA DIVERSITY
     * ========================================================================
     *
     * Beta diversity uses prevalence-filtered ASV counts.
     *
     * filtered_counts.tsv
     *      ↓
     * multiplicative zero replacement
     *      ↓
     * CLR
     *      ↓
     * Aitchison distance
     *      ↓
     * PCoA
     *
     * Optional:
     *      grouped PCoA
     *      PERMANOVA
     *      PERMDISP
     */


    BETA_DIVERSITY(
        ABUNDANCE.out.filtered_counts,
        METADATA_CHECK.out.validated_metadata
    )


    /*
     * ========================================================================
     * 9. QC — MULTIQC
     * ========================================================================
     *
     * MultiQC currently aggregates FastQC reports.
     */

    ch_multiqc_files = Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(
        PREPARE_SAMPLES.out.fastqc_zip
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        QIIME_DENOISING.out.cutadapt_log
    )

    MULTIQC(
        ch_multiqc_files.collect()
    )


    emit:

    multiqc_report =
        MULTIQC.out.report

}
