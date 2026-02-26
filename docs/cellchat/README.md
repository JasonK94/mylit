
# CellChat Pipeline

This directory contains documentation for the CellChat analysis pipeline.

## Overview
The CellChat pipeline is designed to automate the analysis of cell-cell communication using the CellChat package.
It supports both "Proper" (sample-wise followed by aggregation) and "Pooled" analysis strategies.

## Usage

### Wrapper Function
The main wrapper function is `run_cellchat_analysis` located in `myR/R/cci_cellchat_wrapper.R`.

```r
source("myR/R/cci_cellchat_wrapper.R")

cellchat <- run_cellchat_analysis(
  input_data = sobj,
  group.by = "cell_type_column",
  species = "human",
  output_dir = "path/to/output",
  # Proper method:
  split.by = "hos_no",
  aggregate.by = "g3"
)
```

### CLI Usage
The primary entry point is `scripts/cellchat/run_cellchat_cli.R`.

```bash
Rscript scripts/cellchat/run_cellchat_cli.R \
  -i input.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  -o output_dir
```

### Arguments
- `input_data`: Seurat object or path.
- `group.by`: Character string. The metadata column to use for cell labels.
- `split.by`: Column to split samples (e.g. `hos_no`). Recommended for rigorous statistical analysis.
- `aggregate.by`: Column to aggregate/compare groups (e.g. `g3`).
- `species`: "human" or "mouse". Default: "human".
- `db.use`: Subset of CellChatDB (e.g., "Secreted Signaling").
- `output_dir`: Path to save results.
- `n_cores`: Number of cores.
- `prob.threshold`: P-value threshold for `identifyOverExpressedGenes` (default: 0.05).

## Version Compatibility
This pipeline is compatible with CellChat v2.
It includes automatic `updateCellChat()` calls to handle objects created with older versions.

## Outputs
- `cellchat.qs` (each sample/condition)
- `net_circle.pdf`: Interaction circle plot.
- `net_bubble.pdf`: Bubble plot.
- Merged objects for comparison.
