#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop(
        "Usage: abundance_preprocessing.R ",
        "<counts_raw.tsv> <prevalence_threshold>"
    )
}

counts_file <- args[1]
prevalence_threshold <- as.numeric(args[2])

if (is.na(prevalence_threshold) ||
    prevalence_threshold < 0 ||
    prevalence_threshold > 1) {
    stop("prevalence_threshold must be between 0 and 1")
}


# ============================================================
# READ COUNT TABLE
# ============================================================

counts_df <- read.delim(
    counts_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)

if (ncol(counts_df) < 2) {
    stop("Count table must contain one feature column and >=1 sample column")
}

feature_ids <- counts_df[[1]]

if (any(is.na(feature_ids)) || any(feature_ids == "")) {
    stop("Feature table contains empty feature IDs")
}

if (anyDuplicated(feature_ids)) {
    stop("Feature table contains duplicated feature IDs")
}

count_matrix <- as.matrix(counts_df[, -1, drop = FALSE])

storage.mode(count_matrix) <- "numeric"

rownames(count_matrix) <- feature_ids


# ============================================================
# VALIDATE COUNTS
# ============================================================

if (any(is.na(count_matrix))) {
    stop("Count table contains NA or non-numeric values")
}

if (any(count_matrix < 0)) {
    stop("Count table contains negative abundances")
}

n_features_initial <- nrow(count_matrix)
n_samples <- ncol(count_matrix)

if (n_samples == 0) {
    stop("Count table contains zero samples")
}


# ============================================================
# REMOVE ZERO-TOTAL FEATURES
# ============================================================

total_abundance <- rowSums(count_matrix)

nonzero_mask <- total_abundance > 0

count_matrix_nonzero <- count_matrix[
    nonzero_mask,
    ,
    drop = FALSE
]


# ============================================================
# PREVALENCE
# ============================================================

prevalence <- rowSums(count_matrix_nonzero > 0) / n_samples

keep_mask <- prevalence > prevalence_threshold

filtered_counts <- count_matrix_nonzero[
    keep_mask,
    ,
    drop = FALSE
]

if (nrow(filtered_counts) == 0) {
    stop(
        "Prevalence filtering removed all features. ",
        "Review prevalence_threshold."
    )
}


# ============================================================
# VALIDATE SAMPLE TOTALS AFTER FILTERING
# ============================================================

sample_totals <- colSums(filtered_counts)

if (any(sample_totals == 0)) {

    zero_samples <- names(sample_totals)[sample_totals == 0]

    stop(
        "One or more samples have zero total abundance after filtering: ",
        paste(zero_samples, collapse = ", ")
    )
}


# ============================================================
# RELATIVE ABUNDANCE
# ============================================================

relative_abundance <- sweep(
    filtered_counts,
    2,
    sample_totals,
    FUN = "/"
)


# ============================================================
# PREVALENCE REPORT
# ============================================================

all_total_abundance <- rowSums(count_matrix)

all_prevalence <- rowSums(count_matrix > 0) / n_samples

status <- ifelse(
    all_total_abundance == 0,
    "removed_zero_total",
    ifelse(
        all_prevalence <= prevalence_threshold,
        "removed_low_prevalence",
        "retained"
    )
)

prevalence_report <- data.frame(
    FeatureID = rownames(count_matrix),
    TotalAbundance = all_total_abundance,
    Prevalence = all_prevalence,
    Status = status,
    stringsAsFactors = FALSE
)


# ============================================================
# WRITE TABLE FUNCTION
# ============================================================

write_matrix <- function(x, filename) {

    out <- data.frame(
        FeatureID = rownames(x),
        x,
        check.names = FALSE
    )

    write.table(
        out,
        file = filename,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


# ============================================================
# OUTPUTS
# ============================================================

write_matrix(
    filtered_counts,
    "filtered_counts.tsv"
)

write_matrix(
    relative_abundance,
    "relative_abundance.tsv"
)

write.table(
    prevalence_report,
    file = "prevalence_report.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# R OBJECT
# ============================================================

abundance_data <- list(
    counts = filtered_counts,
    relative_abundance = relative_abundance,
    prevalence = prevalence_report,
    prevalence_threshold = prevalence_threshold
)

saveRDS(
    abundance_data,
    file = "abundance_data.rds"
)


# ============================================================
# SUMMARY
# ============================================================

n_zero_removed <- sum(all_total_abundance == 0)

n_low_prevalence_removed <- sum(
    all_total_abundance > 0 &
    all_prevalence <= prevalence_threshold
)

n_retained <- nrow(filtered_counts)

summary_lines <- c(
    "ABUNDANCE PREPROCESSING SUMMARY",
    "============================================================",
    paste("Samples:", n_samples),
    paste("Initial features:", n_features_initial),
    paste("Prevalence threshold:", prevalence_threshold),
    paste("Retention rule: prevalence >", prevalence_threshold),
    paste("Zero-total features removed:", n_zero_removed),
    paste(
        "Low-prevalence features removed:",
        n_low_prevalence_removed
    ),
    paste("Features retained:", n_retained)
)

writeLines(
    summary_lines,
    "abundance_preprocessing_summary.txt"
)

cat(paste(summary_lines, collapse = "\n"))
cat("\n")
