#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        "Usage: beta_diversity.R ",
        "<filtered_counts.tsv> ",
        "<zero_method> ",
        "<zero_fraction> ",
        "<zero_threshold>",
        "<min_positive_samples>"
    )
}

counts_file <- args[1]
zero_method <- args[2]
zero_fraction <- as.numeric(args[3])
zero_threshold <- as.numeric(args[4])

min_positive_samples <- as.integer(args[5])

if (
    is.na(min_positive_samples) ||
    min_positive_samples < 1
) {
    stop(
        "beta_min_positive_samples must be >= 1"
    )
}


# ============================================================
# PARAMETERS
# ============================================================

supported_zero_methods <- c(
    "CZM",
    "GBM",
    "SQ",
    "BL"
)

if (!(zero_method %in% supported_zero_methods)) {
    stop(
        "Unsupported beta_zero_method: ",
        zero_method,
        ". Supported methods: ",
        paste(supported_zero_methods, collapse = ", ")
    )
}

if (
    is.na(zero_fraction) ||
    zero_fraction <= 0 ||
    zero_fraction >= 1
) {
    stop("beta_zero_fraction must be > 0 and < 1")
}

if (
    is.na(zero_threshold) ||
    zero_threshold <= 0
) {
    stop("beta_zero_threshold must be > 0")
}


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(zCompositions)
    library(ggplot2)
})


# ============================================================
# READ COUNTS
# ============================================================

counts_df <- read.delim(
    counts_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)

if (!("FeatureID" %in% colnames(counts_df))) {
    stop(
        "Count table must contain a FeatureID column."
    )
}

if (ncol(counts_df) < 3) {
    stop(
        "Beta diversity requires at least two sample columns."
    )
}

feature_ids <- counts_df$FeatureID

if (
    any(is.na(feature_ids)) ||
    any(feature_ids == "")
) {
    stop("Empty FeatureIDs detected.")
}

if (anyDuplicated(feature_ids)) {
    stop("Duplicated FeatureIDs detected.")
}


# ============================================================
# FEATURE x SAMPLE MATRIX
# ============================================================

counts <- as.matrix(
    counts_df[, setdiff(
        colnames(counts_df),
        "FeatureID"
    ), drop = FALSE]
)

storage.mode(counts) <- "numeric"

rownames(counts) <- feature_ids


# ============================================================
# VALIDATE COUNTS
# ============================================================

if (any(is.na(counts))) {
    stop(
        "Count table contains NA or non-numeric values."
    )
}

if (any(counts < 0)) {
    stop("Count table contains negative values.")
}

if (any(counts %% 1 != 0)) {
    stop(
        "Beta diversity expects integer count data."
    )
}

# ============================================================
# FEATURE ELIGIBILITY FOR BETA DIVERSITY
# ============================================================

feature_totals <- rowSums(counts)

nonzero_feature_mask <- feature_totals > 0

counts <- counts[
    nonzero_feature_mask,
    ,
    drop = FALSE
]


feature_positive_samples <- rowSums(
    counts > 0
)

feature_keep_mask <- (
    feature_positive_samples >=
    min_positive_samples
)


feature_filter_report <- data.frame(
    FeatureID = rownames(counts),
    PositiveSamples = feature_positive_samples,
    Status = ifelse(
        feature_keep_mask,
        "retained",
        "removed_low_positive_samples"
    ),
    stringsAsFactors = FALSE
)


write.table(
    feature_filter_report,
    file = "beta_feature_filter_report.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


counts <- counts[
    feature_keep_mask,
    ,
    drop = FALSE
]


if (nrow(counts) == 0) {
    stop(
        "Beta-diversity filtering removed all features. ",
        "Review beta_min_positive_samples."
    )
}


# ============================================================
# SAMPLE ELIGIBILITY AFTER FEATURE FILTERING
# ============================================================

sample_totals <- colSums(counts)

zero_samples <- names(sample_totals)[
    sample_totals == 0
]


excluded_samples_report <- data.frame(
    SampleID = zero_samples,
    Reason = rep(
        "zero_total_after_beta_feature_filter",
        length(zero_samples)
    ),
    stringsAsFactors = FALSE
)


write.table(
    excluded_samples_report,
    file = "beta_excluded_samples.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


if (length(zero_samples) > 0) {

    counts <- counts[
        ,
        sample_totals > 0,
        drop = FALSE
    ]
}


if (ncol(counts) < 3) {
    stop(
        "Fewer than three samples remain after ",
        "beta-diversity eligibility filtering."
    )
}


# ============================================================
# TRANSPOSE: SAMPLE x FEATURE
# ============================================================

community <- t(counts)

n_samples <- nrow(community)
n_features <- ncol(community)

if (n_samples < 3) {
    stop(
        "At least three samples are required ",
        "for a two-dimensional PCoA."
    )
}


# ============================================================
# ZERO STATISTICS
# ============================================================

n_values <- length(community)
n_zeros <- sum(community == 0)
zero_fraction_observed <- n_zeros / n_values


# ============================================================
# MULTIPLICATIVE ZERO REPLACEMENT
# ============================================================

if (n_zeros > 0) {

    if (zero_method == "CZM") {

        zero_replaced <- zCompositions::cmultRepl(
            community,
            label = 0,
            method = zero_method,
            output = "prop",
            frac = zero_fraction,
            threshold = zero_threshold,
            z.delete = FALSE,
            suppress.print = TRUE
        )

    } else {

        zero_replaced <- zCompositions::cmultRepl(
            community,
            label = 0,
            method = zero_method,
            output = "prop",
            frac = zero_fraction,
            z.delete = FALSE,
            suppress.print = TRUE
        )
    }

    zero_replaced <- as.matrix(
        zero_replaced
    )

} else {

    zero_replaced <- community /
        rowSums(community)
}


if (!identical(dim(zero_replaced), dim(community))) {
    stop(
        "Zero replacement changed matrix dimensions. ",
        "This should not occur when z.delete = FALSE."
    )
}


rownames(zero_replaced) <- rownames(community)
colnames(zero_replaced) <- colnames(community)


# ============================================================
# VALIDATE ZERO REPLACEMENT
# ============================================================

if (any(!is.finite(zero_replaced))) {
    stop(
        "Zero replacement produced non-finite values."
    )
}

if (any(zero_replaced <= 0)) {
    stop(
        "Zero replacement left zero or negative values."
    )
}

composition_sums <- rowSums(
    zero_replaced
)

if (
    any(abs(composition_sums - 1) > 1e-8)
) {
    stop(
        "Zero-replaced compositions do not sum to 1."
    )
}


# ============================================================
# CLR TRANSFORMATION
# ============================================================

log_composition <- log(
    zero_replaced
)

clr <- log_composition -
    rowMeans(log_composition)


# ============================================================
# VALIDATE CLR
# ============================================================

if (any(!is.finite(clr))) {
    stop(
        "CLR transformation produced non-finite values."
    )
}

clr_sums <- rowSums(clr)

if (
    any(abs(clr_sums) > 1e-8)
) {
    stop(
        "CLR rows do not sum approximately to zero."
    )
}


# ============================================================
# AITCHISON DISTANCE
#
# Euclidean distance on CLR coordinates
# ============================================================

aitchison <- dist(
    clr,
    method = "euclidean"
)

distance_matrix <- as.matrix(
    aitchison
)


# ============================================================
# PCoA
# ============================================================

k <- min(
    2,
    n_samples - 1
)

pcoa <- cmdscale(
    aitchison,
    k = k,
    eig = TRUE,
    add = FALSE
)

coordinates <- as.data.frame(
    pcoa$points,
    check.names = FALSE
)

axis_names <- paste0(
    "PCoA",
    seq_len(ncol(coordinates))
)

colnames(coordinates) <- axis_names

coordinates$SampleID <- rownames(
    coordinates
)

coordinates <- coordinates[
    ,
    c(
        "SampleID",
        axis_names
    ),
    drop = FALSE
]


# ============================================================
# VARIANCE EXPLAINED
# ============================================================

eigenvalues <- pcoa$eig

positive_eigenvalues <- eigenvalues[
    eigenvalues > 0
]

positive_total <- sum(
    positive_eigenvalues
)

variance_explained <- if (
    positive_total > 0
) {
    positive_eigenvalues /
        positive_total
} else {
    rep(
        NA_real_,
        length(positive_eigenvalues)
    )
}

variance_table <- data.frame(
    Axis = paste0(
        "PCoA",
        seq_along(positive_eigenvalues)
    ),
    Eigenvalue = positive_eigenvalues,
    ProportionExplained = variance_explained,
    stringsAsFactors = FALSE
)


# ============================================================
# WRITE MATRIX FUNCTION
# ============================================================

write_sample_matrix <- function(
    x,
    filename
) {

    output <- data.frame(
        SampleID = rownames(x),
        x,
        check.names = FALSE
    )

    write.table(
        output,
        file = filename,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


# ============================================================
# OUTPUT NUMERIC TABLES
# ============================================================

write_sample_matrix(
    zero_replaced,
    "zero_replaced_composition.tsv"
)

write_sample_matrix(
    clr,
    "clr_matrix.tsv"
)

write_sample_matrix(
    distance_matrix,
    "aitchison_distance.tsv"
)

write.table(
    coordinates,
    file = "pcoa_coordinates.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    variance_table,
    file = "pcoa_variance.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# UNGROUPED PCoA PLOT
# ============================================================

if (ncol(coordinates) < 3) {
    stop(
        "Two PCoA axes were not generated."
    )
}

pc1_var <- variance_table$ProportionExplained[
    variance_table$Axis == "PCoA1"
]

pc2_var <- variance_table$ProportionExplained[
    variance_table$Axis == "PCoA2"
]

x_label <- if (
    length(pc1_var) == 1 &&
    !is.na(pc1_var)
) {
    sprintf(
        "PCoA1 (%.1f%%)",
        100 * pc1_var
    )
} else {
    "PCoA1"
}

y_label <- if (
    length(pc2_var) == 1 &&
    !is.na(pc2_var)
) {
    sprintf(
        "PCoA2 (%.1f%%)",
        100 * pc2_var
    )
} else {
    "PCoA2"
}


p <- ggplot(
    coordinates,
    aes(
        x = PCoA1,
        y = PCoA2
    )
) +
    geom_point(
        size = 3
    ) +
    geom_text(
        aes(label = SampleID),
        vjust = -0.8,
        check_overlap = TRUE
    ) +
    labs(
        x = x_label,
        y = y_label,
        title = "Aitchison PCoA"
    ) +
    theme_bw()


ggsave(
    "aitchison_pcoa.pdf",
    plot = p,
    width = 8,
    height = 6
)

ggsave(
    "aitchison_pcoa.png",
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
)


# ============================================================
# SUMMARY
# ============================================================

summary_lines <- c(
    "BETA DIVERSITY SUMMARY",
    "============================================================",
    paste("Samples:", n_samples),
    paste("Features:", n_features),
    paste("Observed zeros:", n_zeros),
    paste(
        "Observed zero fraction:",
        signif(
            zero_fraction_observed,
            6
        )
    ),
    paste(
        "Zero replacement method:",
        zero_method
    ),
    paste(
        "Zero replacement fraction:",
        zero_fraction
    ),
    paste(
        "Zero replacement threshold:",
        zero_threshold
    ),
    "",
    "Transformation:",
    "Multiplicative zero replacement -> CLR",
    "",
    "Distance:",
    "Euclidean distance on CLR coordinates (Aitchison)",
    "",
    "Ordination:",
    "Principal Coordinates Analysis (PCoA)"
)

writeLines(
    summary_lines,
    "beta_diversity_summary.txt"
)

cat(
    paste(summary_lines, collapse = "\n"),
    "\n"
)
