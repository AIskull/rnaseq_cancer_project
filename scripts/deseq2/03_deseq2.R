############################################################
# RNA-seq Differential Expression Analysis using DESeq2
#
# Project: RNA-seq Cancer Analysis
# Script: 03_deseq2.R
#
# ------------------------------------------------------------
# Working directory setup (works for both interactive + source)
# ------------------------------------------------------------

set_project_root <- function(project_name = "rnaseq_cancer_project") {

  # Case 1: Running via source(".../03_deseq2.R") -> we can detect file path
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

  if (!is.null(script_file) && nzchar(script_file)) {
    script_dir <- dirname(normalizePath(script_file))
    project_root <- normalizePath(file.path(script_dir, "..", ".."))  # scripts/deseq2 -> project root
    setwd(project_root)
  }

  # Case 2: Interactive / VS Code line-by-line -> rely on current working dir
  # Validate we are in the project root (fail loudly if not)
  if (basename(getwd()) != project_name) {
    stop(
      "Wrong working directory.\n",
      "Expected project root folder named: ", project_name, "\n",
      "Current getwd(): ", getwd(), "\n\n",
      "Fix by running setwd('C:/Users/kurie/rnaseq_cancer_project') first, ",
      "or run source() from the project root."
    )
  }

  message("Working directory: ", getwd())
}

set_project_root()

# Description:
# This script performs differential gene expression analysis
# between tumor and normal samples using RNA-seq count data.
# The analysis is conducted using the DESeq2 Bioconductor
# package and includes normalization, statistical testing,
# and visualization of results.
#
# Input:
# - Gene-level count matrix (rows = genes, columns = samples)
# - Sample metadata table with condition labels
#
# Output:
# - Differential expression results table
# - MA plot
# - PCA plot
# - Heatmap of top differentially expressed genes
#
# Author: Kevin Arredondo
# Date: 2026-01-31
############################################################
# Set project root as working directory
# Set working directory to project root based on script location (robust)

# ------------------------------------------------------------
# Working directory setup (safe for interactive + source)
# ------------------------------------------------------------

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(script_file) && nzchar(script_file)) {
  # Running via source(): setwd to project root automatically
  script_dir <- dirname(normalizePath(script_file))
  project_root <- normalizePath(file.path(script_dir, "..", ".."))
  setwd(project_root)
}

message("Working directory: ", getwd())



# -------------------------------
# 1. Load required libraries
# -------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})



# -------------------------------
# 2. Load count data and metadata
# -------------------------------
# Read count matrix
counts <- read.table(
  file = "data/counts.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# Read sample metadata
metadata <- read.table(
  file = "data/metadata.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Use sample IDs as rownames (DESeq2 expects this)
rownames(metadata) <- metadata$sample_id

# Ensure sample order matches counts columns
metadata <- metadata[colnames(counts), , drop = FALSE]

# Sanity check
stopifnot(all(rownames(metadata) == colnames(counts)))



# -------------------------------
# 3. Create DESeq2 dataset object
# -------------------------------
# Convert condition to factor
metadata$condition <- factor(metadata$condition)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Inspect dataset
dds


# --------------------------------
# 4. Differential expression analysis
# NOTE:
# This project uses a very small toy dataset for demonstration purposes.
# With extremely small datasets, DESeq2's default dispersion trend fitting fails.
# We therefore use gene-wise dispersion estimates directly.

dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)

# Set final dispersions explicitly
dispersions(dds) <- mcols(dds)$dispGeneEst

# Run Wald test
dds <- nbinomWaldTest(dds)



# -------------------------------
# 5. Save results table
# -------------------------------
# Make sure output folders exist
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Get DE results (Tumor vs Normal)
res_tbl <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Convert to data.frame and add gene_id as a column
res_df <- as.data.frame(res_tbl)
res_df$gene_id <- rownames(res_df)

# Reorder columns so gene_id is first
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

# Save as TSV
write.table(
  res_df,
  file = "results/tables/deseq2_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Quick check
head(res_df)


# -------------------------------
# 6. Generate MA plot
# -------------------------------
# Make sure figures folder exists
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Save MA plot
png("results/figures/ma_plot.png", width = 1200, height = 900, res = 150)

plotMA(
  res_tbl,
  ylim = c(-6, 6),
  main = "MA Plot: Tumor vs Normal"
)

dev.off()

# Confirm file exists
file.exists("results/figures/ma_plot.png")

# --------------------------------
# 7. Generate Volcano Plot
# --------------------------------

library(ggplot2)

# Add significance labels
res_df$significance <- "Not Significant"
res_df$significance[
  res_df$padj < 0.05 & res_df$log2FoldChange > 0
] <- "Upregulated in Tumor"

res_df$significance[
  res_df$padj < 0.05 & res_df$log2FoldChange < 0
] <- "Downregulated in Tumor"

# Create volcano plot
volcano_plot <- ggplot(
  res_df,
  aes(x = log2FoldChange, y = -log10(padj))
) +
  geom_point(aes(color = significance), size = 3, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Upregulated in Tumor" = "red",
      "Downregulated in Tumor" = "blue",
      "Not Significant" = "gray"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot: Tumor vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Status"
  ) +
  theme_minimal(base_size = 14)

# Save volcano plot
ggsave(
  filename = "results/figures/volcano_plot.png",
  plot = volcano_plot,
  width = 10,
  height = 8,
  dpi = 150
)

# Display plot in R session
volcano_plot

# --------------------------------
# 8. Simple log2 transform of normalized counts (toy dataset-safe)
norm_counts <- counts(dds, normalized = TRUE)
log_norm_counts <- log2(norm_counts + 1)

# Inspect transformed values
head(log_norm_counts)




# -------------------------------
# 9. PCA visualization (toy dataset-safe)

pca <- prcomp(t(log_norm_counts), scale. = TRUE)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df$condition <- metadata$condition[
  match(pca_df$sample, metadata$sample_id)
]

pvar <- (pca$sdev^2) / sum(pca$sdev^2) * 100

p <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, show.legend = FALSE) +
  xlab(paste0("PC1: ", round(pvar[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(pvar[2], 1), "% variance")) +
  ggtitle("PCA (log2 normalized counts): Tumor vs Normal") +
  theme_minimal(base_size = 14)

ggsave(
  "results/figures/pca_plot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 150
)

p

# -------------------------------
# 10. Heatmap of top DE genes (toy dataset-safe, no pheatmap needed)
# ---------------------------------------------------------------

# Make sure figures folder exists
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Safety checks
stopifnot(exists("res_df"))
stopifnot(exists("log_norm_counts"))
stopifnot(is.matrix(log_norm_counts) || is.data.frame(log_norm_counts))

# Convert log_norm_counts to matrix if needed
log_norm_counts <- as.matrix(log_norm_counts)

# Pick genes to plot:
# Prefer smallest padj; if padj missing, fall back to largest |log2FC|
res_df_clean <- res_df[!is.na(res_df$gene_id), ]

if ("padj" %in% colnames(res_df_clean) && any(!is.na(res_df_clean$padj))) {
  res_df_clean <- res_df_clean[order(res_df_clean$padj), ]
} else {
  # fallback: rank by effect size
  res_df_clean <- res_df_clean[order(-abs(res_df_clean$log2FoldChange)), ]
}

top_n <- min(20, nrow(res_df_clean))
genes_to_plot <- res_df_clean$gene_id[1:top_n]

# Keep only genes that actually exist in the count matrix rownames
genes_to_plot <- intersect(genes_to_plot, rownames(log_norm_counts))

if (length(genes_to_plot) == 0) {
  stop("No genes from res_df matched rownames(log_norm_counts). Check gene_id names.")
}

heat_mat <- log_norm_counts[genes_to_plot, , drop = FALSE]

# Row z-score (standard heatmap style)
heat_scaled <- t(scale(t(heat_mat)))
heat_scaled[is.na(heat_scaled)] <- 0  # in case a row had 0 variance

# Build a sample order by condition (optional but nice)
# (Assumes metadata has: sample_id and condition)
sample_order <- metadata$sample_id[order(metadata$condition)]
sample_order <- intersect(sample_order, colnames(heat_scaled))
heat_scaled <- heat_scaled[, sample_order, drop = FALSE]

# Column side colors by condition
cond_vec <- metadata$condition[match(colnames(heat_scaled), metadata$sample_id)]
cond_fac <- factor(cond_vec)
cond_cols <- setNames(rainbow(length(levels(cond_fac))), levels(cond_fac))
col_side <- cond_cols[as.character(cond_fac)]

# Save heatmap (IMPORTANT: always close device)
out_file <- "results/figures/heatmap_top_genes.png"
png(out_file, width = 1200, height = 900, res = 150)
par(mar = c(8, 8, 4, 2))

heatmap(
  heat_scaled,
  Rowv = NA, Colv = NA,
  scale = "none",
  col = colorRampPalette(c("blue", "white", "red"))(101),
  margins = c(8, 8),
  ColSideColors = col_side,
  main = "Heatmap (row z-score): Top DE genes"
)

legend(
  "topright",
  legend = names(cond_cols),
  fill = cond_cols,
  title = "Condition",
  cex = 1
)

dev.off()

# Confirm
file.exists(out_file)


