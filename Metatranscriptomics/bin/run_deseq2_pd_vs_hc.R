#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(apeglm)
  library(PCAtools)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(RColorBrewer)
  library(statmod)
  library(readr)
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

tximport_rds <- get_arg("--tximport-rds")
samplesheet <- get_arg("--samplesheet")

sample_column <- get_arg(
  "--sample-column",
  "sample_id"
)

condition_column <- get_arg(
  "--condition-column",
  "condition"
)

case_label <- get_arg(
  "--case",
  "PD"
)

control_label <- get_arg(
  "--control",
  "HC"
)

min_sample_counts <- as.numeric(
  get_arg("--min-sample-counts", "100000")
)

min_gene_counts <- as.numeric(
  get_arg("--min-gene-counts", "10")
)

output_results <- get_arg(
  "--output-results",
  "deseq2_pd_vs_hc.tsv"
)

output_unshrunken <- get_arg(
  "--output-unshrunken",
  "deseq2_pd_vs_hc_unshrunken.tsv"
)

output_rds <- get_arg(
  "--output-rds",
  "deseq2_dataset.rds"
)

output_vst <- get_arg(
  "--output-vst",
  "vst_counts.tsv"
)

output_pca <- get_arg(
  "--output-pca",
  "pca_plot.pdf"
)

output_distance_heatmap <- get_arg(
  "--output-distance-heatmap",
  "sample_distance_heatmap.pdf"
)

output_dispersion <- get_arg(
  "--output-dispersion",
  "dispersion_plot.pdf"
)

output_ma <- get_arg(
  "--output-ma",
  "ma_plot.pdf"
)

output_volcano <- get_arg(
  "--output-volcano",
  "volcano_plot.pdf"
)

output_deg_heatmap <- get_arg(
  "--output-deg-heatmap",
  "deg_heatmap.pdf"
)

output_qc_summary <- get_arg(
  "--output-qc-summary",
  "deseq2_qc_summary.tsv"
)

if (is.null(tximport_rds)) {
  stop("Missing --tximport-rds")
}

if (is.null(samplesheet)) {
  stop("Missing --samplesheet")
}

txi <- readRDS(tximport_rds)

metadata <- if (
  grepl(
    "\\.csv$",
    samplesheet,
    ignore.case = TRUE
  )
) {
  read_csv(
    samplesheet,
    show_col_types = FALSE
  )
} else {
  read_tsv(
    samplesheet,
    show_col_types = FALSE
  )
}

if (!sample_column %in% colnames(metadata)) {
  stop(
    paste(
      "Missing sample column:",
      sample_column
    )
  )
}

if (!condition_column %in% colnames(metadata)) {
  stop(
    paste(
      "Missing condition column:",
      condition_column
    )
  )
}

if (anyDuplicated(metadata[[sample_column]]) > 0) {
  stop("Duplicate sample IDs found in metadata")
}

quant_samples <- colnames(txi$counts)
metadata_samples <- metadata[[sample_column]]

missing_metadata <- setdiff(
  quant_samples,
  metadata_samples
)

missing_quant <- setdiff(
  metadata_samples,
  quant_samples
)

if (length(missing_metadata) > 0) {
  stop(
    paste(
      "Samples present in Salmon but missing from metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

if (length(missing_quant) > 0) {
  stop(
    paste(
      "Samples present in metadata but missing from Salmon:",
      paste(missing_quant, collapse = ", ")
    )
  )
}

unexpected_conditions <- setdiff(
  unique(metadata[[condition_column]]),
  c(control_label, case_label)
)

if (length(unexpected_conditions) > 0) {
  stop(
    paste(
      "Unexpected conditions in metadata:",
      paste(unexpected_conditions, collapse = ", ")
    )
  )
}

metadata <- metadata %>%
  mutate(
    condition_for_deseq = factor(
      .data[[condition_column]],
      levels = c(
        control_label,
        case_label
      )
    )
  )

rownames(metadata) <- metadata[[sample_column]]

metadata <- metadata[
  quant_samples,
  ,
  drop = FALSE
]

# Professor: remove samples with total counts < 100K.
sample_total_counts <- colSums(
  txi$counts
)

keep_samples <- (
  sample_total_counts >= min_sample_counts
)

qc_summary <- data.frame(
  sample_id = names(sample_total_counts),
  total_counts = as.numeric(sample_total_counts),
  keep_sample = keep_samples,
  condition = metadata[[condition_column]]
)

write_tsv(
  qc_summary,
  output_qc_summary
)

txi_filtered <- txi

txi_filtered$counts <- txi$counts[
  ,
  keep_samples,
  drop = FALSE
]

txi_filtered$abundance <- txi$abundance[
  ,
  keep_samples,
  drop = FALSE
]

txi_filtered$length <- txi$length[
  ,
  keep_samples,
  drop = FALSE
]

metadata <- metadata[
  keep_samples,
  ,
  drop = FALSE
]

group_counts <- table(
  metadata$condition_for_deseq
)

message("Samples per condition after filtering:")
print(group_counts)

if (any(group_counts < 2)) {
  stop(
    paste(
      "Not enough biological replicates for DESeq2.",
      "At least 2 samples per condition are required after filtering."
    )
  )
}

# Professor: DESeqDataSetFromTximport with design ~ Group.
dds <- DESeqDataSetFromTximport(
  txi_filtered,
  colData = metadata,
  design = ~ condition_for_deseq
)

# Professor: pre-filter genes whose total count is < 10.
keep_genes <- rowSums(
  counts(dds)
) >= min_gene_counts

message(
  "Genes retained after pre-filtering: ",
  sum(keep_genes),
  " / ",
  length(keep_genes)
)

dds <- dds[
  keep_genes,
]

if (nrow(dds) == 0) {
  stop(
    "No genes remained after DESeq2 pre-filtering"
  )
}

# Professor: run DESeq2.
dds <- DESeq(dds)

message("DESeq2 results names:")
print(resultsNames(dds))


# Additional PNG copies for embedding in MultiQC.
# The original PDF outputs are preserved unchanged.
mqc_png_path <- function(path) {
  if (grepl("\\.[^.]+$", path)) {
    sub("\\.[^.]+$", "_mqc.png", path)
  } else {
    paste0(path, "_mqc.png")
  }
}

output_pca_mqc <- mqc_png_path(output_pca)
output_distance_heatmap_mqc <- mqc_png_path(output_distance_heatmap)
output_dispersion_mqc <- mqc_png_path(output_dispersion)
output_deg_heatmap_mqc <- mqc_png_path(output_deg_heatmap)
output_ma_mqc <- mqc_png_path(output_ma)
output_volcano_mqc <- mqc_png_path(output_volcano)

# Professor: diagnostic plot of dispersion estimates.
pdf(
  output_dispersion,
  width = 7,
  height = 5
)

plotDispEsts(dds)

dev.off()

# MultiQC PNG copy of the same plot.
png(
  output_dispersion_mqc,
  width = 7,
  height = 5,
  units = "in",
  res = 160
)

plotDispEsts(dds)

dev.off()

# Professor: PD vs Healthy/control, BH correction.
res_unshrunken <- results(
  dds,
  contrast = c(
    "condition_for_deseq",
    case_label,
    control_label
  ),
  pAdjustMethod = "BH"
)

message("Unshrunken DESeq2 summary:")
print(
  summary(
    res_unshrunken,
    alpha = 0.05
  )
)

# apeglm requires the coefficient name rather than a contrast.
coef_name <- paste0(
  "condition_for_deseq_",
  case_label,
  "_vs_",
  control_label
)

if (!coef_name %in% resultsNames(dds)) {
  stop(
    paste(
      "Expected DESeq2 coefficient was not found:",
      coef_name,
      "Available:",
      paste(resultsNames(dds), collapse = ", ")
    )
  )
}

# Professor: apeglm LFC shrinkage.
res_shrunken <- lfcShrink(
  dds,
  coef = coef_name,
  type = "apeglm"
)

message("Shrunken DESeq2 summary:")
print(
  summary(
    res_shrunken,
    alpha = 0.05
  )
)

res_unshrunken_df <- as.data.frame(
  res_unshrunken
) %>%
  tibble::rownames_to_column(
    "gene_id"
  ) %>%
  arrange(
    padj,
    pvalue
  )

res_shrunken_df <- as.data.frame(
  res_shrunken
) %>%
  tibble::rownames_to_column(
    "gene_id"
  ) %>%
  arrange(
    padj,
    pvalue
  )

# Primary downstream table = professor's shrunken results.
write_tsv(
  res_shrunken_df,
  output_results
)

# Preserve the unshrunken results too.
write_tsv(
  res_unshrunken_df,
  output_unshrunken
)

saveRDS(
  dds,
  output_rds
)

# Professor: VST with blind = FALSE.
vst_data <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vst_data
)

write_tsv(
  as.data.frame(vst_matrix) %>%
    tibble::rownames_to_column("gene_id"),
  output_vst
)

# Professor: PCAtools PCA over VST counts.
pca_metadata <- as.data.frame(
  colData(vst_data)
)

pca_res <- pca(
  vst_matrix,
  metadata = pca_metadata
)

explained_var_percent <- round(
  100 * (
    pca_res$variance /
      sum(pca_res$variance)
  ),
  2
)

xlab <- paste0(
  "PC1 (",
  explained_var_percent[1],
  "%)"
)

ylab <- paste0(
  "PC2 (",
  explained_var_percent[2],
  "%)"
)

pca_colors <- setNames(
  c(
    "steelblue",
    "darkorange"
  ),
  c(
    control_label,
    case_label
  )
)

pca_shapes <- setNames(
  c(
    21,
    22
  ),
  c(
    control_label,
    case_label
  )
)

pdf(
  output_pca,
  width = 7,
  height = 5
)

biplot(
  pca_res,
  colby = "condition_for_deseq",
  colkey = pca_colors,
  shape = "condition_for_deseq",
  shapekey = pca_shapes,
  pointSize = 3,
  legendPosition = "right",
  title = "PCA of VST-normalized Expression",
  xlab = xlab,
  ylab = ylab
)

dev.off()

# MultiQC PNG copy of the same plot.
png(
  output_pca_mqc,
  width = 7,
  height = 5,
  units = "in",
  res = 160
)

biplot(
  pca_res,
  colby = "condition_for_deseq",
  colkey = pca_colors,
  shape = "condition_for_deseq",
  shapekey = pca_shapes,
  pointSize = 3,
  legendPosition = "right",
  title = "PCA of VST-normalized Expression",
  xlab = xlab,
  ylab = ylab
)

dev.off()

# Professor: pairwise sample distances from VST counts.
sample_distances <- dist(
  t(vst_matrix)
)

sample_distance_matrix <- as.matrix(
  sample_distances
)

rownames(sample_distance_matrix) <- colnames(
  vst_data
)

colnames(sample_distance_matrix) <- colnames(
  vst_data
)

annotation_col <- data.frame(
  Group = colData(vst_data)$condition_for_deseq,
  row.names = colnames(vst_data)
)

annotation_colors <- list(
  Group = setNames(
    c(
      "steelblue",
      "darkorange"
    ),
    c(
      control_label,
      case_label
    )
  )
)

pdf(
  output_distance_heatmap,
  width = 7,
  height = 6
)

pheatmap(
  sample_distance_matrix,
  clustering_distance_rows = sample_distances,
  clustering_distance_cols = sample_distances,
  color = colorRampPalette(
    rev(
      brewer.pal(
        9,
        "Blues"
      )
    )
  )(255),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Heatmap of Sample Distances (VST-normalized)",
  fontsize_row = 8,
  fontsize_col = 8
)

dev.off()

# MultiQC PNG copy of the same plot.
png(
  output_distance_heatmap_mqc,
  width = 7,
  height = 6,
  units = "in",
  res = 160
)

pheatmap(
  sample_distance_matrix,
  clustering_distance_rows = sample_distances,
  clustering_distance_cols = sample_distances,
  color = colorRampPalette(
    rev(
      brewer.pal(
        9,
        "Blues"
      )
    )
  )(255),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Heatmap of Sample Distances (VST-normalized)",
  fontsize_row = 8,
  fontsize_col = 8
)

dev.off()

# Professor: DEG heatmap selected from apeglm-shrunken results.
#
# The professor replaces row names with MAG_gene labels after functional
# annotation. Functional annotation is intentionally decoupled from this
# analysis run, so gene_id is retained here. The selected genes, VST values,
# Z-score normalization and clustering inputs remain the same.
deg_df <- res_shrunken_df %>%
  filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < 0.05,
    abs(log2FoldChange) > 1.5
  )

if (nrow(deg_df) == 0) {
  stop(
    paste(
      "No significant DEGs found with padj < 0.05 and",
      "|log2FoldChange| > 1.5"
    )
  )
}

n_top <- 50

up_genes <- deg_df %>%
  filter(log2FoldChange > 0) %>%
  arrange(padj) %>%
  slice_head(n = n_top / 2) %>%
  pull(gene_id)

down_genes <- deg_df %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj) %>%
  slice_head(n = n_top / 2) %>%
  pull(gene_id)

top_genes <- unique(
  c(
    up_genes,
    down_genes
  )
)

message(
  "DEG heatmap: selected ",
  length(up_genes),
  " upregulated and ",
  length(down_genes),
  " downregulated genes"
)

vst_deg_matrix <- vst_matrix[
  top_genes,
  ,
  drop = FALSE
]

zscore_data <- t(
  scale(
    t(vst_deg_matrix)
  )
)

zscore_data <- zscore_data[
  complete.cases(zscore_data),
  ,
  drop = FALSE
]

if (nrow(zscore_data) == 0) {
  stop(
    "No non-constant DEG rows remained after Z-score normalization"
  )
}

deg_annotation_col <- data.frame(
  Group = colData(dds)$condition_for_deseq,
  row.names = colnames(dds)
)

deg_annotation_colors <- list(
  Group = setNames(
    c(
      "#fe9929",
      "#8c6bb1"
    ),
    c(
      control_label,
      case_label
    )
  )
)

deg_heatmap_colors <- colorRampPalette(
  brewer.pal(
    11,
    "RdYlBu"
  )
)(100)

pdf(
  output_deg_heatmap,
  width = 10,
  height = 8
)

pheatmap(
  zscore_data,
  color = deg_heatmap_colors,
  annotation_col = deg_annotation_col,
  annotation_colors = deg_annotation_colors,
  cluster_rows = FALSE,
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = TRUE,
  show_colnames = ncol(zscore_data) <= 30,
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  scale = "none"
)

dev.off()

# MultiQC PNG copy of the same plot.
png(
  output_deg_heatmap_mqc,
  width = 10,
  height = 8,
  units = "in",
  res = 160
)

pheatmap(
  zscore_data,
  color = deg_heatmap_colors,
  annotation_col = deg_annotation_col,
  annotation_colors = deg_annotation_colors,
  cluster_rows = FALSE,
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = TRUE,
  show_colnames = ncol(zscore_data) <= 30,
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  scale = "none"
)

dev.off()

# Professor: MA plot of both unshrunken and shrunken LFC.
pdf(
  output_ma,
  width = 14,
  height = 5
)

par(
  mfrow = c(
    1,
    2
  )
)

plotMA(
  res_unshrunken,
  main = paste(
    "Unshrunken:",
    case_label,
    "vs",
    control_label
  ),
  ylim = c(
    -5,
    5
  )
)

abline(
  h = c(
    -1.5,
    1.5
  ),
  col = "dodgerblue",
  lwd = 2
)

plotMA(
  res_shrunken,
  main = paste(
    "Shrunken:",
    case_label,
    "vs",
    control_label
  ),
  ylim = c(
    -5,
    5
  )
)

abline(
  h = c(
    -1.5,
    1.5
  ),
  col = "dodgerblue",
  lwd = 2
)

dev.off()

# MultiQC PNG copy of the same plot.
png(
  output_ma_mqc,
  width = 14,
  height = 5,
  units = "in",
  res = 160
)

par(
  mfrow = c(
    1,
    2
  )
)

plotMA(
  res_unshrunken,
  main = paste(
    "Unshrunken:",
    case_label,
    "vs",
    control_label
  ),
  ylim = c(
    -5,
    5
  )
)

abline(
  h = c(
    -1.5,
    1.5
  ),
  col = "dodgerblue",
  lwd = 2
)

plotMA(
  res_shrunken,
  main = paste(
    "Shrunken:",
    case_label,
    "vs",
    control_label
  ),
  ylim = c(
    -5,
    5
  )
)

abline(
  h = c(
    -1.5,
    1.5
  ),
  col = "dodgerblue",
  lwd = 2
)

dev.off()

# Professor's volcano is based on the unshrunken results.
volcano_data <- res_unshrunken_df %>%
  mutate(
    log10_padj = -log10(padj),
    significance = case_when(
      padj < 0.05 &
        log2FoldChange > 0 ~ paste(case_label, "Enriched"),

      padj < 0.05 &
        log2FoldChange < 0 ~ paste(control_label, "Enriched"),

      TRUE ~ "Not Significant"
    ),
    label = ifelse(
      padj < 0.05 &
        abs(log2FoldChange) > 1.5,
      gene_id,
      NA
    )
  )

volcano <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = log10_padj
  )
) +
  geom_point(
    aes(
      color = significance,
      fill = significance
    ),
    alpha = 0.7,
    size = 2,
    shape = 21
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = c(
      -1.5,
      1.5
    ),
    linetype = "dotted",
    alpha = 0.3
  ) +
  geom_text_repel(
    aes(label = label),
    na.rm = TRUE,
    max.overlaps = 20,
    size = 3
  ) +
  theme_minimal(
    base_size = 12
  ) +
  labs(
    x = paste0(
      "Log2 Fold Change (",
      case_label,
      " vs ",
      control_label,
      ")"
    ),
    y = "-log10(FDR-adjusted p-value)",
    color = "Significance",
    fill = "Significance"
  )

ggsave(
  output_volcano,
  plot = volcano,
  width = 10,
  height = 8,
  dpi = 300
)

# MultiQC PNG copy of the volcano plot.
ggsave(
  output_volcano_mqc,
  plot = volcano,
  width = 10,
  height = 8,
  dpi = 300
)

message("DESeq2 completed")
message("Shrunken results: ", output_results)
message("Unshrunken results: ", output_unshrunken)
