# =====================================================
# Single-Cell RNA-seq Quality Control Analysis
# Author: Amal
# Purpose:
# Import and process 10X Genomics scRNA-seq data
# using Seurat, generate initial QC metrics, and save
# results and figures for documentation.
# =====================================================


# -------------------------
# Load libraries
# -------------------------

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)


# -------------------------
# Create output folders
# -------------------------

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# -------------------------
# Create Seurat objects
# -------------------------

for (file in c(
  "ctrl_raw_feature_bc_matrix",
  "stim_raw_feature_bc_matrix"
)) {
  
  seurat_data <- Read10X(
    data.dir = paste0("data/", file)
  )
  
  seurat_obj <- CreateSeuratObject(
    counts = seurat_data,
    min.features = 100,
    project = file
  )
  
  assign(file, seurat_obj)
}


# -------------------------
# Merge datasets
# -------------------------

merged_seurat <- merge(
  x = ctrl_raw_feature_bc_matrix,
  y = stim_raw_feature_bc_matrix,
  add.cell.id = c("ctrl", "stim")
)

merged_seurat <- JoinLayers(merged_seurat)


# -------------------------
# Add sample metadata
# -------------------------

merged_seurat$sample <- ifelse(
  grepl("^ctrl_", colnames(merged_seurat)),
  "Control",
  "Stimulated"
)


# -------------------------
# Calculate mitochondrial percentage
# -------------------------

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
  merged_seurat,
  pattern = "^MT-"
)


# -------------------------
# Inspect metadata
# -------------------------

head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)


# -------------------------
# Save merged Seurat object
# -------------------------

saveRDS(
  merged_seurat,
  "results/merged_seurat.rds"
)


# -------------------------
# Save overall summary
# -------------------------

sink("results/summary.txt")

cat("Single-Cell RNA-seq QC Summary\n")
cat("================================\n\n")

cat("Number of genes:",
    nrow(merged_seurat), "\n")

cat("Number of cells:",
    ncol(merged_seurat), "\n")

cat("Average transcripts per cell:",
    mean(merged_seurat$nCount_RNA), "\n")

cat("Average genes per cell:",
    mean(merged_seurat$nFeature_RNA), "\n")

cat("Average mitochondrial percentage:",
    mean(merged_seurat$percent.mt), "\n")

sink()


# -------------------------
# Save sample-level summary table
# -------------------------

qc_summary <- merged_seurat@meta.data %>%
  group_by(sample) %>%
  summarise(
    cells = n(),
    average_transcripts = mean(nCount_RNA),
    average_genes = mean(nFeature_RNA),
    average_mitochondrial_percentage = mean(percent.mt)
  )

write.csv(
  qc_summary,
  "results/qc_summary_by_sample.csv",
  row.names = FALSE
)


# -------------------------
# QC violin plot by sample
# -------------------------

qc_violin_by_sample <- VlnPlot(
  merged_seurat,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3
)

ggsave(
  "figures/qc_violin_by_sample.png",
  plot = qc_violin_by_sample,
  width = 12,
  height = 5
)


# -------------------------
# Scatter plot: counts vs detected genes
# -------------------------

qc_counts_vs_features <- FeatureScatter(
  merged_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "sample"
)

ggsave(
  "figures/qc_counts_vs_features.png",
  plot = qc_counts_vs_features,
  width = 7,
  height = 5
)


# -------------------------
# Save metadata preview
# -------------------------

metadata_preview <- head(merged_seurat@meta.data, 20)

write.csv(
  metadata_preview,
  "results/metadata_preview.csv"
)


# -------------------------
# Print completion message
# -------------------------

cat("Analysis completed successfully.\n")
cat("Outputs saved in results/ and figures/ folders.\n")