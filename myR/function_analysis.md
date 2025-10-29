# Function Comparison Analysis

## Overview
This document compares functions from different commits and branches in the myR package to track changes and refactoring efforts.

## Key Functions Comparison

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

### Utility Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `convert_gene_ids` | genes, from, to | Converted gene IDs | Converts gene identifiers between formats |
| `get_barcode_mapping` | sample_names | barcode mapping | Maps barcode IDs to samples |
| `prepare_gene_lists` | deg_df, fc_threshold, ... | Gene lists (up/down) | Prepares gene lists from DEG results |
| `sort_samples` | samples | sorted samples | Sorts sample names by numeric pattern |

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

## Functions from Refactoring Branches (refactoring-v0.2 and refactor)

### New Signature Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `score_signature` | expr_data, signature, normalize | signature scores | Scores gene signatures from expression data |
| `CompareModuleScoringMethods` | seurat_object, gene_modules, methods | Comparison results | Compares different module scoring methods |

### Deprecated/Legacy Functions

The following functions have legacy versions in the legacy.R file:
- `pseudobulk_linear_fit_legacy`
- `post_hoc_slope_comparison_legacy`
- `run_pseudobulk_deg_legacy`
- `myhm_genesets2_legacy`

## Summary Statistics

- **Total Functions in new_master**: ~177 functions
- **Core Analysis Functions**: ~50 functions
- **Visualization Functions**: ~30 functions
- **Utility Functions**: ~20 functions
- **Legacy Functions**: ~25 functions (in legacy.R)

## Refactoring Summary

The refactoring branches (refactoring-v0.2 and refactor) introduced:
1. Reorganization of functions into core, deprecated, signatures, and utilities directories
2. New signature scoring functions
3. Improved documentation and parameter handling
4. Consolidation of DEG methods

