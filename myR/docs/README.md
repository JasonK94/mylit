# myR Package

## Overview

`myR` is a comprehensive R package for single-cell RNA sequencing (scRNA-seq) analysis, with specialized support for GeoMx spatial transcriptomics data. The package provides functions for differential expression analysis, pseudotime inference, pathway enrichment, cell-cell interaction analysis, and comprehensive visualization tools.

## Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install the package from local directory
devtools::install("path/to/myR", dependencies = TRUE)
```

## Key Features

- **Differential Expression Analysis**: Pseudobulk and LMM-based DEG analysis
- **Pathway Enrichment**: GO, KEGG, and GSEA analyses
- **Pseudotime Analysis**: Trajectory inference using Slingshot and Monocle3
- **Cell-Cell Interactions**: NicheNet-based interaction analysis
- **Module Scoring**: PROGENy and other signature scoring methods
- **GeoMx Analysis**: Specialized tools for NanoString GeoMx data

## Main Components

### Core Analysis Modules

- **pseudobulk_deg.R**: Pseudobulk differential expression analysis
- **pathway_analysis.R**: Pathway enrichment analysis (GO, KEGG, GSEA)
- **pseudotime.R**: Trajectory and pseudotime analysis
- **CCI.R**: Cell-cell interaction analysis
- **deg.R**: Module scoring and DEG analysis

### Utility Modules

- **plots.R**: Comprehensive visualization functions
- **utils.R**: Utility functions for data handling
- **gene_list.R**: Gene list management and operations
- **cluster_frequency.R**: Cluster proportion analysis

## Function Comparison Table

### Core Analysis Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `prepare_geomx_data` | count_file, metadata files | List (expression matrix, metadata, gene_info) | Prepares GeoMx data for analysis by extracting expression matrix and metadata |
| `q3_normalize` | expr_matrix, scaling_factor | normalized matrix (log2) | Normalizes expression by Q3 (75th percentile) |
| `find_deg_geomx` | norm_expr, metadata, group_var | DEG results data frame | Finds differentially expressed genes using limma/edgeR |
| `run_lmm_multiple_genes` | seurat_obj, genes, config | LMM results | Runs Linear Mixed Models for multiple genes |
| `find_response_differential_genes` | lmm_summary, pval_threshold | Significant genes data frame | Identifies genes with significant treatment response |

### Pseudobulk and DEG Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `run_pseudobulk_deg` | seurat_obj, analysis_level, cluster_group, condition_col | DEG results | Performs pseudobulk differential expression analysis |
| `prepare_pseudobulk_edgeR` | seurat_obj, cluster_group, sample_col, counts_assay | Pseudobulk matrix | Prepares data for edgeR pseudobulk analysis |
| `cluster_pseudobulk_deg` | sobj, cluster_group, condition_col, genes, ... | DEG results per cluster | Cluster-specific pseudobulk DEG analysis |
| `pseudobulk_linear_fit` | sobj, genes, sample_col, numeric_predictor, ... | Linear fit results | Fits linear models on pseudobulk data |

### Pathway Analysis Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `myGO` | DEG data frame, pathway_set, analysis_type | Named list of pathway results | Comprehensive GO/GSEA pathway enrichment analysis |
| `run_go_analysis` | genes, gene_type, ont, ... | GO results | Gene Ontology enrichment analysis |
| `run_kegg_analysis` | genes, background, pval_cutoff | KEGG results | KEGG pathway enrichment |
| `run_gsea_analysis` | ranked_genes, pathway_set, ... | GSEA results | Gene Set Enrichment Analysis |

### Module Scoring and Signatures

| Function Name | Input | Output | Function |
|---|---|---|---|
| `AddMultipleModuleScores` | seurat_object, gene_modules | Seurat object with module scores | Adds multiple module scores to Seurat object |
| `add_progeny_scores` | seurat_obj, organism, topn | Seurat object with PROGENy scores | Adds PROGENy pathway activity scores |
| `add_signature_enrichit` | seurat_obj, gene_sets, ... | Seurat object with signature scores | Adds gene signature enrichment scores |
| `find_gene_signature` | data, signature_type, ... | Gene signature list | Finds gene signatures from data |

### Visualization Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `plot_volcano` | lmm_summary, ... | Volcano plot | Creates volcano plot from LMM results |
| `PlotModuleScoreHeatmap` | seurat_object, features, assays | Heatmap | Plots module score heatmap |
| `plot_cluster_fractions` | sobj_metadata, cluster_col, ... | Cluster fraction plots | Visualizes cluster proportions |
| `plot_interaction_for_gene` | sobj, gene, patient, treatment, timepoint | Interaction plot | Plots gene expression interaction |
| `myhm_genes4` | sobj, features, group.by, ... | Heatmap | Creates gene expression heatmap |
| `cdf` | data, probability_col, ratio_col | Cumulative distribution plot | Plots cumulative distribution |

### Pseudotime Analysis Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `run_slingshot_from_seurat` | seurat_obj, cluster_col, reduced_dim_name | SingleCellExperiment object | Runs Slingshot trajectory inference |
| `analyze_gene_dynamics` | gene_id, cds_obj, condition_col_name, ... | Analysis results | Analyzes gene dynamics along pseudotime |
| `process_gene_list_dynamics` | gene_list, cds_obj, condition_col_name, ... | List of results | Processes multiple genes for dynamics analysis |
| `analyze_gene_dynamics_tradeSeq` | gene_id, cds_obj, condition_col_name, ... | tradeSeq results | Analyzes gene dynamics using tradeSeq |

### Cell-Cell Interaction Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| Functions in CCI.R | Multiple | NicheNet analysis results | Analyses cell-cell interactions using NicheNet |
| `run_nichenet_analysis` | sobj, sender, receiver, condition_col | Interaction results | Runs NicheNet ligand-receptor analysis |

## Usage Examples

### Pseudobulk DEG Analysis

```r
# Load the package
library(myR)

# Run pseudobulk DEG analysis
deg_results <- run_pseudobulk_deg(
  seurat_obj = seurat_object,
  analysis_level = "per_cluster",
  cluster_group = "cell_type",
  condition_col = "treatment"
)
```

### Pathway Enrichment Analysis

```r
# Run comprehensive pathway analysis
pathway_results <- myGO(
  DEG = deg_results,
  analysis_type = "ALL",
  fc_threshold = 0.25,
  pval_threshold = 0.05
)
```

### Module Scoring

```r
# Add PROGENy pathway scores
seurat_obj <- add_progeny_scores(
  seurat_obj = seurat_object,
  organism = "Human",
  topn = 100
)
```

## Dependencies

- Seurat
- edgeR
- limma
- slingshot
- monocle3
- NicheNet
- PROGENy
- clusterProfiler
- msigdbr
- dplyr
- ggplot2

## Documentation

For detailed documentation on each function, see the package help files:

```r
?run_pseudobulk_deg
?myGO
?add_progeny_scores
```

## Authors

- Development team

## License

See LICENSE file for details.

