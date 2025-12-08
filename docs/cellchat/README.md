
# CellChat Pipeline

This directory contains documentation for the CellChat analysis pipeline.

## Overview
The CellChat pipeline is designed to automate the analysis of cell-cell communication using the CellChat package.
It takes a Seurat object and metadata columns as input and generates:
- CellChat object with computed probabilities
- Visualization plots (Circle plots, Bubble plots)

## Usage

### Wrapper Function
The main wrapper function is `run_cellchat_analysis` located in `myR/R/cci_cellchat_wrapper.R`.

```r
source("myR/R/cci_cellchat_wrapper.R")

cellchat <- run_cellchat_analysis(
  sobj = sobj,
  group.by = "cell_type_column",
  species = "human",
  output_dir = "path/to/output"
)
```

### Arguments
- `sobj`: Seurat object.
- `group.by`: Character string. The metadata column to use for cell labels (e.g., cell types).
- `species`: Character string. "human" or "mouse". Default: "human".
- `db.use`: Character vector. Subset of CellChatDB to use (e.g., "Secreted Signaling"). If NULL, uses the entire database.
- `assay_name`: Character string. Assay to use. Default: "SCT".
- `output_dir`: Character string. Path to save results.
- `n_cores`: Integer. Number of cores for parallel processing. Default: 4.

## Outputs
- `cellchat_object.rds`: The computed CellChat object (also saved as .qs if qs package is available).
- `net_visual_circle.pdf`: Circle plots of interaction counts and weights.
- `net_visual_bubble.pdf`: Bubble plots of signaling pathways.

## Example
See `scripts/cellchat/test_cellchat.R` for a running example.
