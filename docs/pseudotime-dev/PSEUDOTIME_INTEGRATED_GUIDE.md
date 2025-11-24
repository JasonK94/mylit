# Pseudotime Analysis Module Integrated Guide

This document is the integrated guide for the Pseudotime Analysis module. It describes the process from Trajectory Inference to Gene Dynamics analysis.

## 1. Introduction

### Purpose
Infers the continuity of cell differentiation paths or state changes (Trajectory) in single-cell data and analyzes how gene expression changes along a virtual time axis (Pseudotime).

### Key Tools
*   **Monocle3**: Widely used Trajectory Inference tool. Strong in inferring complex branching structures.
*   **Slingshot**: Flexible and powerful Trajectory Inference.
*   **GAM (Generalized Additive Models)**: Models non-linear gene expression patterns along Pseudotime.

## 2. Data Preparation

### Preprocessing
Before Pseudotime analysis, metadata cleaning and format conversion are essential.
The `preprocess_pseudotime_data()` function performs the following:
*   **sex**: Standardize to `M`, `F` Factor.
*   **g3**: Convert to Factor (`1`, `2`).
*   **Time variables**: Parse `icu_adm_dt`, `ia_start`, etc. as datetime format and calculate time differences (Duration).

### Key Metadata
*   `nih_change`: NIH score change (clinical indicator).
*   `sex`: Sex.
*   `g3`: Group variable.
*   `arrival_gcs_score`, `age`: Clinical covariates.

## 3. Functions & Usage

### Trajectory Inference
1.  **`run_monocle3_from_seurat(sobj)`**:
    *   Converts Seurat object to Monocle3 object.
    *   Performs Preprocessing, Reduction, Cluster Cells, Learn Graph, Order Cells.
    *   Root cell can be specified automatically or manually.

2.  **`run_slingshot_from_seurat(sobj)`**:
    *   Runs Slingshot from Seurat object.
    *   Calculates Lineages and Pseudotime.

### Gene Dynamics Analysis
1.  **`analyze_gene_dynamics()`**:
    *   Models gene expression changes along Pseudotime using GAM.
    *   Calculates p-values and visualizes.

2.  **`detect_expression_patterns()`**:
    *   Classifies gene expression patterns: Increasing, Decreasing, Oscillatory, Bimodal, etc.

3.  **`analyze_metadata_correlation_dynamics()`**:
    *   Analyzes changes in correlations between gene expression and metadata by Pseudotime windows.

## 4. Workflow Example

```r
# 1. Data preprocessing
source("scripts/pseudotime-dev/preprocess_data.R")
sobj_prep <- preprocess_pseudotime_data(
  input_file = "/data/user3/sobj/IS_scvi_251107_ds2500.qs",
  output_file = "/data/user3/sobj/IS_scvi_251107_prep.qs"
)

# 2. Basic analysis (Trajectory + Gene Dynamics)
source("scripts/pseudotime-dev/test_pseudotime_basic.R")
# Script calls run_monocle3_from_seurat() and analyze_gene_dynamics() internally
```

## 5. Appendix

### Script Locations
*   `scripts/pseudotime-dev/preprocess_data.R`: Preprocessing function.
*   `scripts/pseudotime-dev/test_pseudotime_basic.R`: Basic analysis test.
*   `scripts/pseudotime-dev/test_preprocessing.R`: Preprocessing test.

### Key Features (Target Genes)
`DDIT4`, `UTY`, `S100B`, `XIST`, `HLA-B`, `CCL4`, `HLA-C`, `TXNIP`, etc.

