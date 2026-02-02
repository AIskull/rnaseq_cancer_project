# RNA-seq Differential Expression Analysis (Tumor vs Normal)

## Project Overview

This project demonstrates a complete RNA-seq differential gene expression analysis
pipeline comparing **Tumor** and **Normal** samples using the **DESeq2** framework.

A small toy dataset is intentionally used to illustrate best practices for
normalization, dispersion estimation, statistical testing, and visualization
in a fully reproducible workflow. Although the dataset is small, the analysis
mirrors the structure and logic used in real RNA-seq studies.

---

## Goals

- Perform differential gene expression analysis (Tumor vs Normal)
- Identify significantly upregulated and downregulated genes
- Generate standard RNA-seq visualizations:
  - MA plot
  - Volcano plot
  - Principal Component Analysis (PCA)
  - Heatmap of top differentially expressed genes
- Build a clean, reproducible, and well-documented RNA-seq analysis pipeline

---

## Tools and Methods

- R / Bioconductor
  - DESeq2 (differential expression analysis)
- Visualization
  - ggplot2
  - pheatmap
- Development Environment
  - Visual Studio Code

---

## Project Structure

RNASEQ_CANCER_PROJECT/
├── data/
│   ├── counts.tsv
│   └── metadata.tsv
├── scripts/
│   └── deseq2/
│       └── 03_deseq2.R
├── results/
│   ├── tables/
│   │   └── deseq2_results.tsv
│   └── figures/
│       ├── ma_plot.png
│       ├── volcano_plot.png
│       ├── pca_plot.png
│       └── heatmap_top_genes.png
└── README.md

---

## How to Run

Requirements:
- R (≥ 4.2)
- Bioconductor
- Packages: DESeq2, ggplot2, pheatmap, dplyr, tibble

Install Bioconductor if needed:

install.packages("BiocManager")
BiocManager::install("DESeq2")

Install remaining packages:

install.packages(c("ggplot2", "pheatmap", "dplyr", "tibble"))

---

## Input Data

data/counts.tsv
- Gene-level raw RNA-seq count matrix
- Rows = genes
- Columns = samples

data/metadata.tsv
- sample_id
- condition (Tumor or Normal)

IMPORTANT:
Sample IDs must exactly match column names in counts.tsv.

---

## Run the Analysis

From the project root directory:

source("scripts/deseq2/03_deseq2.R")

Pipeline steps:
- Data loading and validation
- Normalization and dispersion estimation
- Differential expression testing (Wald test)
- Figure generation
- Results export

---

## Output Files

Tables:
results/tables/deseq2_results.tsv
- Log2 fold change
- Wald statistic
- Raw p-value
- Adjusted p-value (FDR)
- Significance labels

Figures:
results/figures/ma_plot.png
results/figures/volcano_plot.png
results/figures/pca_plot.png
results/figures/heatmap_top_genes.png

---

## Notes

- Toy dataset used for demonstration purposes
- Statistical results are not biologically interpretable
- Pipeline mirrors real RNA-seq workflows and scales easily

---

## Methods

Data Input:
Counts and metadata were provided as tab-delimited files.
Sample order was explicitly matched prior to analysis.

Differential Expression:
DESeq2 was used with size-factor normalization.
Gene-wise dispersion estimates were used due to small dataset size.
Wald test was applied for statistical testing.

Statistical Outputs:
- Log2 fold change
- Wald test statistic
- Raw p-value
- Benjamini–Hochberg adjusted p-value (FDR)

Gene categories:
- Upregulated in Tumor
- Downregulated in Tumor
- Not Significant

---

## Results Summary

MA plot, volcano plot, PCA, and heatmap show consistent separation
between Tumor and Normal samples. While this dataset is small,
the workflow demonstrates a complete, reproducible RNA-seq
differential expression pipeline.

---

## Reproducibility

All analyses are scripted.
Running the same script with the same inputs reproduces all results.

---

## Future Work

- Scale to real RNA-seq datasets
- Add QC with FastQC / MultiQC
- Integrate tximport
- Support batch correction and multifactor designs
- Add pathway enrichment analysis
- Containerize with Docker or Singularity
- Automate with Snakemake or Nextflow
