#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        "Usage: beta_diversity_plot.R ",
        "<pcoa_coordinates.tsv> ",
        "<pcoa_variance.tsv> ",
        "<metadata.tsv> ",
        "<group_column>"
    )
}

coordinates_file <- args[1]
variance_file <- args[2]
metadata_file <- args[3]
group_column <- args[4]


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(ggplot2)
})


# ============================================================
# READ PCoA
# ============================================================

coordinates <- read.delim(
    coordinates_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

variance <- read.delim(
    variance_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_coordinates <- c(
    "SampleID",
    "PCoA1",
    "PCoA2"
)

if (!all(
    required_coordinates %in%
        colnames(coordinates)
)) {
    stop(
        "PCoA coordinates must contain: ",
        paste(
            required_coordinates,
            collapse = ", "
        )
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
        grepl(
            "^#q2:",
            first_value
        )
    ) {
        metadata <- metadata[
            -1,
            ,
            drop = FALSE
        ]
    }
}


# ============================================================
# SAMPLE ID
# ============================================================

colnames(metadata)[1] <- "SampleID"

if (!(group_column %in% colnames(metadata))) {
    stop(
        "Requested beta_group_column '",
        group_column,
        "' was not found in metadata. Available columns: ",
        paste(
            colnames(metadata),
            collapse = ", "
        )
    )
}


# ============================================================
# MATCH SAMPLES
# ============================================================

idx <- match(
    coordinates$SampleID,
    metadata$SampleID
)

if (any(is.na(idx))) {

    missing <- coordinates$SampleID[
        is.na(idx)
    ]

    stop(
        "PCoA samples missing from metadata: ",
        paste(
            missing,
            collapse = ", "
        )
    )
}

plot_data <- coordinates

plot_data$Group <- metadata[
    idx,
    group_column
]


if (any(is.na(plot_data$Group))) {
    stop(
        "Missing values detected in beta grouping column."
    )
}


# ============================================================
# AXIS LABELS
# ============================================================

get_variance <- function(axis_name) {

    value <- variance$ProportionExplained[
        variance$Axis == axis_name
    ]

    if (
        length(value) == 1 &&
        is.finite(value)
    ) {
        return(value)
    }

    NA_real_
}


pc1_var <- get_variance("PCoA1")
pc2_var <- get_variance("PCoA2")

x_label <- if (!is.na(pc1_var)) {
    sprintf(
        "PCoA1 (%.1f%%)",
        100 * pc1_var
    )
} else {
    "PCoA1"
}

y_label <- if (!is.na(pc2_var)) {
    sprintf(
        "PCoA2 (%.1f%%)",
        100 * pc2_var
    )
} else {
    "PCoA2"
}


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
    plot_data,
    aes(
        x = PCoA1,
        y = PCoA2,
        colour = Group
    )
) +
    geom_point(
        size = 4
    ) +
    geom_text(
        aes(label = SampleID),
        vjust = -0.8,
        check_overlap = TRUE,
        show.legend = FALSE
    ) +
    labs(
        x = x_label,
        y = y_label,
        colour = group_column,
        title = paste(
            "Aitchison PCoA by",
            group_column
        )
    ) +
    theme_bw()


ggsave(
    "aitchison_pcoa_by_group.pdf",
    plot = p,
    width = 8,
    height = 6
)

ggsave(
    "aitchison_pcoa_by_group.png",
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
)
