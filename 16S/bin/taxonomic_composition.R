#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        "Usage: taxonomic_composition.R ",
        "<relative_abundance.tsv> ",
        "<taxonomy.tsv> ",
        "<taxonomic_rank> ",
        "<rare_threshold>"
    )
}

abundance_file <- args[1]
taxonomy_file <- args[2]
taxonomic_rank <- args[3]
rare_threshold <- as.numeric(args[4])


# ============================================================
# VALIDATE PARAMETERS
# ============================================================

valid_ranks <- c(
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
)

if (!(taxonomic_rank %in% valid_ranks)) {
    stop(
        "Invalid taxonomic rank: ", taxonomic_rank,
        ". Allowed values: ",
        paste(valid_ranks, collapse = ", ")
    )
}

if (
    is.na(rare_threshold) ||
    rare_threshold < 0 ||
    rare_threshold >= 1
) {
    stop(
        "composition_rare_threshold must be >= 0 and < 1"
    )
}


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})


# ============================================================
# READ RELATIVE ABUNDANCE
# ============================================================

abundance <- read.delim(
    abundance_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)

if (ncol(abundance) < 2) {
    stop(
        "Relative abundance table must contain FeatureID ",
        "and at least one sample."
    )
}

if (!("FeatureID" %in% colnames(abundance))) {
    stop(
        "Relative abundance table must contain a FeatureID column."
    )
}

if (anyDuplicated(abundance$FeatureID)) {
    stop("Duplicated FeatureIDs detected in abundance table.")
}

sample_columns <- setdiff(
    colnames(abundance),
    "FeatureID"
)

abundance_matrix <- as.matrix(
    abundance[, sample_columns, drop = FALSE]
)

storage.mode(abundance_matrix) <- "numeric"

if (any(is.na(abundance_matrix))) {
    stop(
        "Relative abundance table contains NA or non-numeric values."
    )
}

if (any(abundance_matrix < 0)) {
    stop(
        "Relative abundance table contains negative values."
    )
}


# ============================================================
# READ TAXONOMY
# ============================================================

taxonomy <- read.delim(
    taxonomy_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)

required_taxonomy_columns <- c(
    "Feature ID",
    "Taxon"
)

if (!all(
    required_taxonomy_columns %in% colnames(taxonomy)
)) {
    stop(
        "Taxonomy table must contain columns: Feature ID and Taxon"
    )
}

if (anyDuplicated(taxonomy[["Feature ID"]])) {
    stop("Duplicated Feature IDs detected in taxonomy.")
}


# ============================================================
# CHECK FEATURE CORRESPONDENCE
# ============================================================

missing_taxonomy <- setdiff(
    abundance$FeatureID,
    taxonomy[["Feature ID"]]
)

if (length(missing_taxonomy) > 0) {
    stop(
        "Features without taxonomy detected: ",
        paste(head(missing_taxonomy, 10), collapse = ", ")
    )
}


# ============================================================
# TAXONOMIC RANK
# ============================================================

rank_positions <- c(
    Domain = 1,
    Phylum = 2,
    Class = 3,
    Order = 4,
    Family = 5,
    Genus = 6,
    Species = 7
)

rank_position <- rank_positions[[taxonomic_rank]]


extract_rank <- function(taxon_string, position) {

    if (
        is.na(taxon_string) ||
        taxon_string == ""
    ) {
        return("Unclassified")
    }

    levels <- strsplit(
        taxon_string,
        ";",
        fixed = TRUE
    )[[1]]

    levels <- trimws(levels)

    if (length(levels) < position) {
        return("Unclassified")
    }

    taxon <- levels[position]

    taxon <- sub(
        "^[A-Za-z]__",
        "",
        taxon
    )

    taxon <- trimws(taxon)

    if (
        is.na(taxon) ||
        taxon == "" ||
        tolower(taxon) %in% c(
            "uncultured",
            "unclassified",
            "unknown"
        )
    ) {
        return("Unclassified")
    }

    taxon
}


taxonomy$SelectedTaxon <- vapply(
    taxonomy$Taxon,
    extract_rank,
    FUN.VALUE = character(1),
    position = rank_position
)


# ============================================================
# JOIN ABUNDANCE + TAXONOMY
# ============================================================

combined <- abundance %>%
    left_join(
        taxonomy %>%
            select(
                `Feature ID`,
                SelectedTaxon
            ),
        by = c(
            "FeatureID" = "Feature ID"
        )
    )

if (any(is.na(combined$SelectedTaxon))) {
    stop(
        "Taxonomy join generated missing taxonomic assignments."
    )
}


# ============================================================
# LONG FORMAT
# ============================================================

composition_long <- combined %>%
    pivot_longer(
        cols = all_of(sample_columns),
        names_to = "SampleID",
        values_to = "RelativeAbundance"
    )


# ============================================================
# COLLAPSE ASVs TO SELECTED TAXONOMIC RANK
# ============================================================

composition_taxa <- composition_long %>%
    group_by(
        SampleID,
        SelectedTaxon
    ) %>%
    summarise(
        RelativeAbundance = sum(RelativeAbundance),
        .groups = "drop"
    )


# ============================================================
# IDENTIFY RARE TAXA
#
# A taxon is considered rare when its MEAN relative abundance
# across all samples is below the configured threshold.
# ============================================================

taxon_means <- composition_taxa %>%
    group_by(SelectedTaxon) %>%
    summarise(
        MeanRelativeAbundance = mean(RelativeAbundance),
        .groups = "drop"
    )

rare_taxa <- taxon_means %>%
    filter(
        MeanRelativeAbundance < rare_threshold
    ) %>%
    pull(SelectedTaxon)


# ============================================================
# GROUP RARE TAXA AS "Other"
# ============================================================

composition_final <- composition_taxa %>%
    mutate(
        Taxon = ifelse(
            SelectedTaxon %in% rare_taxa,
            "Other",
            SelectedTaxon
        )
    ) %>%
    group_by(
        SampleID,
        Taxon
    ) %>%
    summarise(
        RelativeAbundance = sum(RelativeAbundance),
        .groups = "drop"
    )


# ============================================================
# VALIDATE COMPOSITION SUMS
# ============================================================

sample_sums <- composition_final %>%
    group_by(SampleID) %>%
    summarise(
        Total = sum(RelativeAbundance),
        .groups = "drop"
    )

if (any(abs(sample_sums$Total - 1) > 1e-8)) {
    stop(
        "Taxonomic composition does not sum to 1 for all samples."
    )
}


# ============================================================
# OUTPUT TABLES
# ============================================================

write.table(
    composition_final,
    file = "taxonomic_composition.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    taxon_means,
    file = "taxon_mean_abundance.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# SAMPLE-LEVEL STACKED BARPLOT
# ============================================================

p <- ggplot(
    composition_final,
    aes(
        x = SampleID,
        y = RelativeAbundance,
        fill = Taxon
    )
) +
    geom_col(
        width = 0.9
    ) +
    scale_y_continuous(
        labels = scales::percent
    ) +
    labs(
        x = "Sample",
        y = "Relative abundance",
        fill = taxonomic_rank,
        title = paste(
            "Taxonomic composition at",
            taxonomic_rank,
            "level"
        )
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        legend.position = "right"
    )


ggsave(
    "taxonomic_composition_by_sample.pdf",
    plot = p,
    width = 12,
    height = 7
)

ggsave(
    "taxonomic_composition_by_sample.png",
    plot = p,
    width = 12,
    height = 7,
    dpi = 300
)


# ============================================================
# SUMMARY
# ============================================================

n_taxa_before <- length(
    unique(composition_taxa$SelectedTaxon)
)

n_rare <- length(
    unique(rare_taxa)
)

n_taxa_after <- length(
    unique(composition_final$Taxon)
)

summary_lines <- c(
    "TAXONOMIC COMPOSITION SUMMARY",
    "============================================================",
    paste("Samples:", length(sample_columns)),
    paste("Input ASVs:", nrow(abundance)),
    paste("Taxonomic rank:", taxonomic_rank),
    paste("Rare threshold:", rare_threshold),
    paste(
        "Rare definition: mean relative abundance <",
        rare_threshold
    ),
    paste("Taxa before rare grouping:", n_taxa_before),
    paste("Taxa grouped as rare:", n_rare),
    paste(
        "Displayed categories after grouping:",
        n_taxa_after
    )
)

writeLines(
    summary_lines,
    "taxonomic_composition_summary.txt"
)

cat(
    paste(summary_lines, collapse = "\n"),
    "\n"
)
