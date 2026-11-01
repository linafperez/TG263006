#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: taxonomic_composition_grouped.R ",
        "<taxonomic_composition.tsv> ",
        "<metadata.tsv> ",
        "<group_column>"
    )
}

composition_file <- args[1]
metadata_file <- args[2]
group_column <- args[3]


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})


# ============================================================
# READ COMPOSITION
# ============================================================

composition <- read.delim(
    composition_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_columns <- c(
    "SampleID",
    "Taxon",
    "RelativeAbundance"
)

if (!all(
    required_columns %in% colnames(composition)
)) {
    stop(
        "Composition table must contain: ",
        paste(required_columns, collapse = ", ")
    )
}


# ============================================================
# READ METADATA
# ============================================================

metadata <- read.delim(
    metadata_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)


# Remove optional QIIME2 #q2:types row

if (nrow(metadata) > 0) {

    first_value <- as.character(
        metadata[1, 1]
    )

    if (
        !is.na(first_value) &&
        grepl("^#q2:", first_value)
    ) {
        metadata <- metadata[-1, , drop = FALSE]
    }
}


# ============================================================
# NORMALIZE SAMPLE ID COLUMN
# ============================================================

colnames(metadata)[1] <- "SampleID"


# ============================================================
# VALIDATE GROUP COLUMN
# ============================================================

if (!(group_column %in% colnames(metadata))) {
    stop(
        "Requested composition_group_column '",
        group_column,
        "' was not found in metadata. Available columns: ",
        paste(colnames(metadata), collapse = ", ")
    )
}


# ============================================================
# JOIN
# ============================================================

plot_data <- composition %>%
    left_join(
        metadata,
        by = "SampleID"
    )

if (any(is.na(plot_data[[group_column]]))) {

    missing_samples <- unique(
        plot_data$SampleID[
            is.na(plot_data[[group_column]])
        ]
    )

    stop(
        "Missing group metadata for samples: ",
        paste(missing_samples, collapse = ", ")
    )
}


# ============================================================
# GROUP MEAN COMPOSITION
# ============================================================

grouped <- plot_data %>%
    group_by(
        .data[[group_column]],
        Taxon
    ) %>%
    summarise(
        MeanRelativeAbundance = mean(RelativeAbundance),
        .groups = "drop"
    )

colnames(grouped)[1] <- "Group"


# ============================================================
# VALIDATE SUMS
# ============================================================

group_sums <- grouped %>%
    group_by(Group) %>%
    summarise(
        Total = sum(MeanRelativeAbundance),
        .groups = "drop"
    )

if (any(abs(group_sums$Total - 1) > 1e-8)) {
    stop(
        "Grouped taxonomic composition does not sum to 1."
    )
}


# ============================================================
# TABLE
# ============================================================

write.table(
    grouped,
    file = "taxonomic_composition_by_group.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
    grouped,
    aes(
        x = Group,
        y = MeanRelativeAbundance,
        fill = Taxon
    )
) +
    geom_col(
        width = 0.8
    ) +
    scale_y_continuous(
        labels = scales::percent
    ) +
    labs(
        x = group_column,
        y = "Mean relative abundance",
        fill = "Taxon",
        title = "Mean taxonomic composition by group"
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
    "taxonomic_composition_by_group.pdf",
    plot = p,
    width = 10,
    height = 7
)

ggsave(
    "taxonomic_composition_by_group.png",
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
)
