#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})


get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)

  if (is.na(idx)) {
    return(default)
  }

  if (idx == length(args)) {
    stop(paste("Missing value for", flag))
  }

  args[[idx + 1]]
}


deseq_results <- get_arg(
  "--deseq-results"
)

mag_map_file <- get_arg(
  "--mag-map"
)

output_mag_mapped <- get_arg(
  "--output-mag-mapped",
  "deseq2_mag_mapped_results.tsv"
)

output_significant <- get_arg(
  "--output-significant",
  "deseq2_significant_results.tsv"
)

output_summary <- get_arg(
  "--output-summary",
  "expression_summary.tsv"
)

output_mag_summary <- get_arg(
  "--output-mag-summary",
  "mag_expression_summary.tsv"
)


if (is.null(deseq_results)) {
  stop("Missing --deseq-results")
}

if (is.null(mag_map_file)) {
  stop("Missing --mag-map")
}


res <- read_tsv(
  deseq_results,
  show_col_types = FALSE
)

mag_map <- read_tsv(
  mag_map_file,
  show_col_types = FALSE
)


required_result_columns <- c(
  "gene_id",
  "log2FoldChange",
  "padj"
)

if (!all(
  required_result_columns %in% colnames(res)
)) {
  stop(
    paste(
      "DESeq2 results must contain:",
      paste(
        required_result_columns,
        collapse = ", "
      )
    )
  )
}


required_map_columns <- c(
  "gene_id",
  "mag_id"
)

if (!all(
  required_map_columns %in% colnames(mag_map)
)) {
  stop(
    "MAG map must contain columns: gene_id and mag_id"
  )
}


gene_to_mag <- mag_map %>%
  select(
    gene_id,
    mag_id
  ) %>%
  distinct()


conflicting_genes <- gene_to_mag %>%
  count(gene_id) %>%
  filter(n > 1)

if (nrow(conflicting_genes) > 0) {
  stop(
    paste(
      "Some gene IDs map to multiple MAGs. Examples:",
      paste(
        head(
          conflicting_genes$gene_id,
          10
        ),
        collapse = ", "
      )
    )
  )
}


mag_mapped <- res %>%
  left_join(
    gene_to_mag,
    by = "gene_id"
  )


significant <- mag_mapped %>%
  filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj <= 0.05,
    abs(log2FoldChange) > 1.5
  )


summary_table <- data.frame(
  metric = c(
    "total_genes",
    "genes_with_mag_mapping",
    "genes_without_mag_mapping",
    "genes_with_padj",
    "significant_genes",
    "significant_upregulated",
    "significant_downregulated"
  ),
  value = c(
    nrow(mag_mapped),
    sum(!is.na(mag_mapped$mag_id)),
    sum(is.na(mag_mapped$mag_id)),
    sum(!is.na(mag_mapped$padj)),
    nrow(significant),
    sum(
      significant$log2FoldChange > 1.5
    ),
    sum(
      significant$log2FoldChange < -1.5
    )
  )
)


mag_summary <- mag_mapped %>%
  filter(
    !is.na(mag_id)
  ) %>%
  group_by(
    mag_id
  ) %>%
  summarise(
    total_genes = n(),
    significant_genes = sum(
      !is.na(padj) &
      !is.na(log2FoldChange) &
      padj <= 0.05 &
      abs(log2FoldChange) > 1.5
    ),
    mean_log2FoldChange = if (
      all(is.na(log2FoldChange))
    ) {
      NA_real_
    } else {
      mean(
        log2FoldChange,
        na.rm = TRUE
      )
    },
    .groups = "drop"
  )


write_tsv(
  mag_mapped,
  output_mag_mapped
)

write_tsv(
  significant,
  output_significant
)

write_tsv(
  summary_table,
  output_summary
)

write_tsv(
  mag_summary,
  output_mag_summary
)


message("Expression summary tables completed")
message(
  "Total genes: ",
  nrow(mag_mapped)
)
message(
  "Genes mapped to MAGs: ",
  sum(!is.na(mag_mapped$mag_id))
)
message(
  "Significant genes: ",
  nrow(significant)
)
