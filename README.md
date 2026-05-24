# Single-Cell RNA-seq Quality Control Analysis

## Description

This project is a practical learning exercise focused on understanding the basics of single-cell RNA sequencing (scRNA-seq) analysis using R and Seurat. The exercise follows a workshop-based workflow and explores how single-cell datasets are imported, structured, and prepared for downstream quality control analysis.

## Purpose

The purpose of this exercise is to gain familiarity with commonly used bioinformatics workflows and tools while building a foundational understanding of how computational methods are applied to biological data.

## Workflow

The following steps were completed during this exercise:

1. Load required R packages
2. Import 10X Genomics count data
3. Create Seurat objects
4. Explore cell metadata
5. Load multiple datasets using a loop
6. Merge control and stimulated datasets
7. Prepare data for later quality control analysis

## Dataset

The exercise uses PBMC (Peripheral Blood Mononuclear Cell) single-cell RNA-seq data from control and interferon beta-stimulated samples in 10X Genomics format.

## Tools Used

- R
- RStudio
- Seurat
- SingleCellExperiment
- tidyverse
- Matrix

## Learning Outcomes

During this exercise I learned and explored:

- Sparse matrices and why they are used in scRNA-seq
- Seurat objects and object-based data structures in R
- Cell metadata and cell barcodes
- Reading and organizing 10X Genomics data
- Creating and merging datasets
- Basic R syntax and workflow structure

## Project Structure

```text
single_cell_rnaseq/
├── data
├── figures
├── report
├── results
├── scripts
│   └── quality_control.R
├── README.md
└── single_cell_rnaseq.Rproj