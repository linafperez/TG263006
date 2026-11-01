#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
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

salmon_dirs_file <- get_arg("--salmon-dirs-file")
tx2gene_file <- get_arg("--tx2gene")

output_rds <- get_arg(
  "--output-rds",
  "tximport_results.rds"
)

counts_output <- get_arg(
  "--counts-output",
  "tximport_counts.tsv"
)

abundance_output <- get_arg(
  "--abundance-output",
  "tximport_tpm.tsv"
)

length_output <- get_arg(
  "--length-output",
  "tximport_length.tsv"
)

if (is.null(salmon_dirs_file)) {
  stop("Missing --salmon-dirs-file")
}

if (is.null(tx2gene_file)) {
  stop("Missing --tx2gene")
}

salmon_dirs <- readLines(salmon_dirs_file)
salmon_dirs <- salmon_dirs[nchar(salmon_dirs) > 0]

sample_ids <- basename(salmon_dirs)

if (anyDuplicated(sample_ids) > 0) {
  stop(
    paste(
      "Duplicate Salmon sample IDs:",
      paste(
        unique(sample_ids[duplicated(sample_ids)]),
        collapse = ", "
      )
    )
  )
}

files <- file.path(salmon_dirs, "quant.sf")
names(files) <- sample_ids

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing quant.sf files:",
      paste(missing_files, collapse = ", ")
    )
  )
}

message("Found Salmon quant files:")
print(files)

tx2gene_raw <- read_tsv(
  tx2gene_file,
  show_col_types = FALSE
)

required_columns <- c(
  "transcript_id",
  "gene_id"
)

if (!all(required_columns %in% colnames(tx2gene_raw))) {
  stop(
    "tx2gene map must contain columns: transcript_id and gene_id"
  )
}

# Same scientific input used by the professor:
# transcript_id -> gene_id.
# Do not deduplicate or otherwise transform the mapping here.
tx2gene <- tx2gene_raw %>%
  select(transcript_id, gene_id)

message(
  "tx2gene mapping rows: ",
  nrow(tx2gene)
)

txi <- tximport(
  files = files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  importer = function(path) {
    readr::read_tsv(
      path,
      show_col_types = FALSE,
      progress = FALSE,
      num_threads = 1
    )
  }
)

saveRDS(
  txi,
  output_rds
)

write_tsv(
  as.data.frame(txi$counts) %>%
    tibble::rownames_to_column("gene_id"),
  counts_output
)

write_tsv(
  as.data.frame(txi$abundance) %>%
    tibble::rownames_to_column("gene_id"),
  abundance_output
)

write_tsv(
  as.data.frame(txi$length) %>%
    tibble::rownames_to_column("gene_id"),
  length_output
)

message("tximport completed")
message("Saved: ", output_rds)
