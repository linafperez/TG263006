#!/usr/bin/env Rscript


# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        "Usage: permanova.R ",
        "<aitchison_distance.tsv> ",
        "<metadata.tsv> ",
        "<group_column> ",
        "<permutations>"
    )
}

distance_file <- args[1]
metadata_file <- args[2]
group_column <- args[3]
n_permutations <- as.integer(args[4])


# ============================================================
# VALIDATE PARAMETERS
# ============================================================

if (
    is.na(n_permutations) ||
    n_permutations < 99
) {
    stop(
        "permanova_permutations must be >= 99"
    )
}


# ============================================================
# DEPENDENCIES
# ============================================================

suppressPackageStartupMessages({
    library(vegan)
})


# ============================================================
# READ DISTANCE MATRIX
# ============================================================

distance_df <- read.delim(
    distance_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!("SampleID" %in% colnames(distance_df))) {
    stop(
        "Distance matrix must contain SampleID."
    )
}

sample_ids <- distance_df$SampleID

distance_matrix <- as.matrix(
    distance_df[
        ,
        setdiff(
            colnames(distance_df),
            "SampleID"
        ),
        drop = FALSE
    ]
)

storage.mode(distance_matrix) <- "numeric"

rownames(distance_matrix) <- sample_ids


# ============================================================
# VALIDATE DISTANCE MATRIX
# ============================================================

if (
    nrow(distance_matrix) !=
    ncol(distance_matrix)
) {
    stop(
        "Distance matrix must be square."
    )
}

if (
    !identical(
        rownames(distance_matrix),
        colnames(distance_matrix)
    )
) {
    stop(
        "Distance matrix row/column sample IDs do not match."
    )
}

if (any(!is.finite(distance_matrix))) {
    stop(
        "Distance matrix contains non-finite values."
    )
}

if (
    max(
        abs(
            distance_matrix -
            t(distance_matrix)
        )
    ) > 1e-8
) {
    stop(
        "Distance matrix is not symmetric."
    )
}

if (
    any(
        abs(diag(distance_matrix)) >
            1e-8
    )
) {
    stop(
        "Distance matrix diagonal is not zero."
    )
}


distance_object <- as.dist(
    distance_matrix
)


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


colnames(metadata)[1] <- "SampleID"


if (!(group_column %in% colnames(metadata))) {
    stop(
        "Requested permanova_column '",
        group_column,
        "' was not found in metadata."
    )
}


# ============================================================
# ALIGN METADATA
# ============================================================

idx <- match(
    sample_ids,
    metadata$SampleID
)

if (any(is.na(idx))) {

    missing <- sample_ids[
        is.na(idx)
    ]

    stop(
        "Distance-matrix samples missing from metadata: ",
        paste(
            missing,
            collapse = ", "
        )
    )
}


analysis_metadata <- metadata[
    idx,
    ,
    drop = FALSE
]

group <- analysis_metadata[[group_column]]

if (any(is.na(group)) || any(group == "")) {
    stop(
        "Missing values detected in PERMANOVA grouping variable."
    )
}

group <- factor(group)

if (nlevels(group) < 2) {
    stop(
        "PERMANOVA requires at least two groups."
    )
}

group_counts <- table(group)

if (any(group_counts < 2)) {

    problematic <- names(
        group_counts[
            group_counts < 2
        ]
    )

    stop(
        "PERMANOVA/PERMDISP requires at least two samples ",
        "per group in this pipeline. Groups failing this rule: ",
        paste(
            problematic,
            collapse = ", "
        )
    )
}


analysis_metadata$PERMANOVA_GROUP <- group


# ============================================================
# PERMANOVA
# ============================================================

set.seed(1)

permanova <- vegan::adonis2(
    distance_object ~ PERMANOVA_GROUP,
    data = analysis_metadata,
    permutations = n_permutations
)

permanova_table <- as.data.frame(
    permanova
)

permanova_table$Term <- rownames(
    permanova_table
)

rownames(permanova_table) <- NULL

permanova_table <- permanova_table[
    ,
    c(
        "Term",
        setdiff(
            colnames(permanova_table),
            "Term"
        )
    ),
    drop = FALSE
]

write.table(
    permanova_table,
    file = "permanova_results.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PERMDISP
# ============================================================

dispersion <- vegan::betadisper(
    distance_object,
    group,
    type = "median"
)

set.seed(1)

dispersion_test <- permutest(
    dispersion,
    permutations = n_permutations
)

permdisp_table <- as.data.frame(
    dispersion_test$tab
)

permdisp_table$Term <- rownames(
    permdisp_table
)

rownames(permdisp_table) <- NULL

permdisp_table <- permdisp_table[
    ,
    c(
        "Term",
        setdiff(
            colnames(permdisp_table),
            "Term"
        )
    ),
    drop = FALSE
]

write.table(
    permdisp_table,
    file = "permdisp_results.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# DISTANCE TO GROUP MEDIAN
# ============================================================

dispersion_by_sample <- data.frame(
    SampleID = sample_ids,
    Group = as.character(group),
    DistanceToGroupMedian =
        dispersion$distances,
    stringsAsFactors = FALSE
)

write.table(
    dispersion_by_sample,
    file = "dispersion_by_sample.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# SUMMARY
# ============================================================

permanova_p <- permanova_table[
    permanova_table$Term ==
        "PERMANOVA_GROUP",
    "Pr(>F)"
]

permanova_r2 <- permanova_table[
    permanova_table$Term ==
        "PERMANOVA_GROUP",
    "R2"
]

permdisp_p_column <- grep(
    "^Pr\\(>F\\)$",
    colnames(permdisp_table),
    value = TRUE
)

permdisp_p <- if (
    length(permdisp_p_column) == 1
) {
    permdisp_table[
        1,
        permdisp_p_column
    ]
} else {
    NA_real_
}


summary_lines <- c(
    "PERMANOVA / PERMDISP SUMMARY",
    "============================================================",
    paste(
        "Grouping column:",
        group_column
    ),
    paste(
        "Groups:",
        paste(
            names(group_counts),
            group_counts,
            sep = "=",
            collapse = ", "
        )
    ),
    paste(
        "Permutations:",
        n_permutations
    ),
    "",
    paste(
        "PERMANOVA R2:",
        permanova_r2
    ),
    paste(
        "PERMANOVA p-value:",
        permanova_p
    ),
    paste(
        "PERMDISP p-value:",
        permdisp_p
    )
)

writeLines(
    summary_lines,
    "permanova_summary.txt"
)

cat(
    paste(summary_lines, collapse = "\n"),
    "\n"
)
