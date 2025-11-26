# Analysis Module Integrated Guide

This document is the integrated guide for the `analysis` module. It describes the differential expression analysis (DEA) pipeline using NEBULA and muscat.

## 1. Introduction

### Purpose
The goal is to accurately detect differentially expressed genes (DEG) between groups in single-cell RNA sequencing data, accounting for subject-level variability.

### Key Features
*   **Mixed-Effects Model**: Models inter-subject variation (Random Effect) using the NEBULA package.
*   **Flexible Formula**: Supports complex experimental designs (interactions, covariate adjustment).
*   **Interactive/Batch Execution**: Provides various execution modes from lightweight testing to full data analysis.

### Formula 1
The Formula currently used for standard analysis is as follows:
```r
~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)
```
*   **g3**: Main group variable (Target)
*   **sex**: Sex (Covariate)
*   **anno3.scvi**: Cell type (Covariate)
*   **GEM**: Batch (treated as Fixed Effect)
*   **(1|GEM/hos_no)**: Nested Random Effect (patients are nested within GEM)

## 2. Workflow

### 2.1 Analysis Pipeline

1.  **Data Loading**: Load Seurat object (`.qs`)
2.  **Data Validation**: Check required columns (`g3`, `hos_no`, `anno3.scvi`, etc.) and missing values
3.  **Formula Parsing**: Separate Random Effects and Fixed Effects, construct design matrix
4.  **NEBULA Execution**: Fit Negative Binomial Mixed Model for each gene
5.  **Result Saving**: Save result summary and statistics in `.qs` format

### 2.2 Execution Methods

**1. Lightweight Test (Recommended)**
Quickly validate pipeline operation.
```bash
Rscript scripts/analysis/test_formula1_light.R
```

**2. Full Data Analysis**
Perform analysis on the entire dataset (may take several hours).
```bash
Rscript scripts/analysis/run_formula1_analysis.R
```

**3. Interactive Execution**
Call functions directly in an R session to perform analysis.
```r
source("scripts/analysis/run_formula1_interactive.R")
result <- run_formula1_analysis_interactive()
```

## 3. Methodology

### NEBULA Model Structure
NEBULA (Negative Binomial mixed-effects model) follows the following formula:

$$Y_{ij} \sim \text{NB}(\mu_{ij}, \phi)$$
$$\log(\mu_{ij}) = \log(\text{offset}_{ij}) + X_{ij}^T \beta + b_i$$

*   $Y_{ij}$: Count for gene $j$, cell $i$
*   $X_{ij}$: Fixed Effects (group, sex, batch, etc.)
*   $b_i$: Random Effect (subject-specific variation, $b_i \sim N(0, \sigma^2_b)$)

### Complete Separation Problem
*   **Phenomenon**: Cases where no data exists for certain Fixed Effect combinations (e.g., only g3=1 patients exist in a specific GEM).
*   **Problem**: The design matrix becomes singular, making coefficient estimation impossible.
*   **Solution**: 
    *   Complete separation between Random Effects and Fixed Effects is not a problem (since it's inter-subject comparison, not intra-subject variation).
    *   If complete separation occurs between Fixed Effects, remove the corresponding covariate or filter the dataset.

## 4. Troubleshooting

### 1. "No random effects in Formula"
*   Check if the Formula contains `(1|...)` syntax.
*   NEBULA requires at least one Random Effect.

### 2. "Design matrix is singular"
*   Complete separation problem has occurred. Simplify interaction terms or increase `min_count` to exclude sparse genes.

### 3. Memory shortage / Long execution time
*   Test first with `light_test = TRUE`.
*   Increase `min_count` to 20 or higher.
*   Consider dividing large datasets by cluster for analysis.

## 5. Appendix

### Key Parameters
*   `min_count`: Minimum number of expressing cells (default: 10). Genes expressed in fewer cells are excluded.
*   `remove_na_cells`: Whether to remove cells containing NA values (default: TRUE).
*   `layer`: Assay Layer to analyze (default: "counts").

### Result File Structure
*   `result$summary`: Per-gene statistics (p-value, logFC, etc.).
*   `result$formula`: Formula used.
*   `result$fixed_effects`: List of estimated Fixed Effect coefficients.

