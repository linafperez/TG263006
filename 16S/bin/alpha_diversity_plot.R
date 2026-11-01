#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: alpha_diversity_plot.R ",
        "<alpha_diversity_long.tsv> ",
        "<metadata.tsv> ",
        "<group_column>"
    )
}

alpha_file <- args[1]
metadata_file <- args[2]
group_column <- args[3]


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(ggplot2)
})


# ============================================================
# READ ALPHA DIVERSITY
# ============================================================

alpha <- read.delim(
    alpha_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


required_alpha <- c(
    "SampleID",
    "Metric",
    "Value"
)

if (!all(required_alpha %in% colnames(alpha))) {
    stop(
        "Alpha diversity table must contain: ",
        paste(required_alpha, collapse = ", ")
    )
}


# ============================================================
# READ QIIME2 METADATA
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

    if (!is.na(first_value) &&
        grepl("^#q2:", first_value)) {

        metadata <- metadata[-1, , drop = FALSE]
    }
}


# ============================================================
# NORMALIZE SAMPLE ID COLUMN
# ============================================================

metadata_id_column <- colnames(metadata)[1]

colnames(metadata)[1] <- "SampleID"


# ============================================================
# VALIDATE GROUP COLUMN
# ============================================================

if (!(group_column %in% colnames(metadata))) {
    stop(
        "Requested alpha_group_column '",
        group_column,
        "' was not found in metadata. Available columns: ",
        paste(colnames(metadata), collapse = ", ")
    )
}


# ============================================================
# MERGE
# ============================================================

plot_data <- merge(
    alpha,
    metadata,
    by = "SampleID",
    all.x = TRUE,
    sort = FALSE
)


if (any(is.na(plot_data[[group_column]]))) {

    missing_samples <- unique(
        plot_data$SampleID[
            is.na(plot_data[[group_column]])
        ]
    )

    stop(
        "Missing metadata values for samples: ",
        paste(missing_samples, collapse = ", ")
    )
}


# ============================================================
# PLOT
# ============================================================

plot_data$GroupForPlot <- factor(
    plot_data[[group_column]]
)

p <- ggplot(
    plot_data,
    aes(
        x = GroupForPlot,
        y = Value
    )
) +
    geom_boxplot(
        outlier.shape = NA
    ) +
    geom_jitter(
        width = 0.12,
        height = 0
    ) +
    facet_wrap(
        ~ Metric,
        scales = "free_y"
    ) +
    labs(
        x = group_column,
        y = "Alpha diversity",
        title = "Alpha diversity by group"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    )


# ============================================================
# OUTPUT
# ============================================================

ggsave(
    "alpha_diversity_by_group.pdf",
    plot = p,
    width = 10,
    height = 7
)

ggsave(
    "alpha_diversity_by_group.png",
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
)
