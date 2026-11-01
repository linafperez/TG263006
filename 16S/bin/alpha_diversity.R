#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: alpha_diversity.R <counts_raw.tsv>")
}

counts_file <- args[1]


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(vegan)
})


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
    stop(
        "Count table must contain one FeatureID column ",
        "and at least one sample column."
    )
}

feature_ids <- counts_df[[1]]

if (any(is.na(feature_ids)) || any(feature_ids == "")) {
    stop("Count table contains empty feature IDs.")
}

if (anyDuplicated(feature_ids)) {
    stop("Count table contains duplicated feature IDs.")
}


# ============================================================
# BUILD FEATURE x SAMPLE MATRIX
# ============================================================

counts <- as.matrix(
    counts_df[, -1, drop = FALSE]
)

storage.mode(counts) <- "numeric"

rownames(counts) <- feature_ids


# ============================================================
# VALIDATE COUNTS
# ============================================================

if (any(is.na(counts))) {
    stop("Count table contains NA or non-numeric values.")
}

if (any(counts < 0)) {
    stop("Count table contains negative abundances.")
}

if (any(counts %% 1 != 0)) {
    stop(
        "Alpha diversity expects raw integer counts, ",
        "but non-integer values were detected."
    )
}

sample_totals <- colSums(counts)

if (any(sample_totals == 0)) {

    zero_samples <- names(sample_totals)[
        sample_totals == 0
    ]

    stop(
        "Samples with zero total abundance detected: ",
        paste(zero_samples, collapse = ", ")
    )
}


# ============================================================
# TRANSPOSE:
# vegan expects samples x features
# ============================================================

community <- t(counts)


# ============================================================
# OBSERVED RICHNESS
# ============================================================

observed <- rowSums(
    community > 0
)


# ============================================================
# CHAO1
# ============================================================

richness_estimates <- vegan::estimateR(
    community
)

chao1 <- richness_estimates[
    "S.chao1",
]


# ============================================================
# SHANNON
# ============================================================

shannon <- vegan::diversity(
    community,
    index = "shannon"
)


# ============================================================
# SIMPSON
# ============================================================

simpson <- vegan::diversity(
    community,
    index = "simpson"
)


# ============================================================
# RESULT TABLE
# ============================================================

sample_ids <- rownames(community)

alpha <- data.frame(
    SampleID = sample_ids,
    Observed = as.numeric(observed),
    Chao1 = as.numeric(chao1),
    Shannon = as.numeric(shannon),
    Simpson = as.numeric(simpson),
    stringsAsFactors = FALSE,
    check.names = FALSE
)


# ============================================================
# VALIDATE RESULTS
# ============================================================

if (any(is.na(alpha$Observed))) {
    stop("NA values detected in Observed richness.")
}

if (any(is.na(alpha$Chao1))) {
    warning("NA values detected in Chao1.")
}

if (any(is.na(alpha$Shannon))) {
    stop("NA values detected in Shannon diversity.")
}

if (any(is.na(alpha$Simpson))) {
    stop("NA values detected in Simpson diversity.")
}


# ============================================================
# WRITE WIDE TABLE
# ============================================================

write.table(
    alpha,
    file = "alpha_diversity.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# LONG FORMAT
# ============================================================

metrics <- c(
    "Observed",
    "Chao1",
    "Shannon",
    "Simpson"
)

alpha_long <- do.call(
    rbind,
    lapply(
        metrics,
        function(metric) {
            data.frame(
                SampleID = alpha$SampleID,
                Metric = metric,
                Value = alpha[[metric]],
                stringsAsFactors = FALSE
            )
        }
    )
)

write.table(
    alpha_long,
    file = "alpha_diversity_long.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# SUMMARY
# ============================================================

summary_lines <- c(
    "ALPHA DIVERSITY SUMMARY",
    "============================================================",
    paste("Samples:", nrow(alpha)),
    paste("Features:", ncol(community)),
    "",
    "Input:",
    "Taxonomically filtered raw ASV counts.",
    "",
    "Metrics:",
    "Observed = number of detected ASVs",
    "Chao1 = richness estimator",
    "Shannon = Shannon diversity index",
    "Simpson = Simpson diversity index",
    "",
    "No prevalence filtering or relative-abundance transformation",
    "was applied before alpha-diversity calculation."
)

writeLines(
    summary_lines,
    "alpha_diversity_summary.txt"
)

cat(
    paste(summary_lines, collapse = "\n"),
    "\n"
)
