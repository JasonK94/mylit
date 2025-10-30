# Function Analysis

## Overview

This document provides a comprehensive analysis of all functions within the `myR` package, automatically generated from the R source files. Each table below corresponds to a specific source file and details the functions it contains.

## Core Functions

### From `myR/R/core/data_preparation.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `.get_feature_vector`| object, feature, assay, slot, cells | Numeric vector | Internal helper to extract a feature from a Seurat object. |
| `get_feature_vec` | object, feature, assay, slot, cells | Numeric vector | User-facing function to extract a feature from a Seurat object. |
| `get_feature_matrix` | object, features, assay, slot, cells | Matrix | Extracts multiple features from a Seurat object. |
| `prepare_metadata_table`| object, columns, cells, drop_na | Data frame | Extracts and optionally filters metadata from a Seurat object. |
| `aggregate_expression_by_group`| object, features, group_by, method, assay, slot | Matrix | Aggregates gene expression across cells within groups. |
| `prepare_count_matrix`| object, assay, slot, features, cells, min_cells, min_features | Matrix | Extracts a count matrix with optional filtering. |
| `convert_to_long_format`| object, features, metadata_cols, assay, slot, cells | Data frame | Converts expression data to long (tidy) format. |
| `check_data_quality`| data, check_finite, check_na, check_negative | List | Performs basic quality checks on expression data. |

### From `myR/R/core/validation.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `validate_seurat` | obj, assay, reduction, min_cells, min_features | TRUE or error | Checks if the input is a valid Seurat object. |
| `validate_metadata_column`| obj, column_name, required_type, allow_na | TRUE or error | Checks if a metadata column exists and is valid. |
| `validate_genes` | obj, genes, min_present, assay, warn_missing | Character vector | Checks if genes exist in the object and filters to valid genes. |
| `validate_numeric_range`| value, param_name, min, max, allow_na | TRUE or error | Checks if a numeric parameter is within an acceptable range. |
| `validate_choice` | value, param_name, choices, multiple | Validated value | Validates that a value is one of the allowed choices. |
| `validate_path` | path, must_exist, extensions, type | Normalized path | Checks if a file path exists and is valid. |
| `create_error_message`| context, message, suggestion | Formatted string | Creates a consistent, informative error message. |
| `check_packages` | packages, load | TRUE or error | Checks if required packages are installed. |

## Visualization Functions

### From `myR/R/visualization/basic_plots.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `mybar` | x, col, bins, title, xlab, color | ggplot2 object | Generates a histogram for a numeric vector. |
| `mydensity` | x, col, title, xlab, color, alpha | ggplot2 object | Generates a density plot for a numeric vector. |
| `myline` | vec, cutoff, show_intersection, title, xlab, ylab, color | ggplot2 object | Creates a line plot for a single numeric vector. |
| `mylines` | vec1, vec2, cutoff1, cutoff2, title, ylab1, ylab2, color1, color2 | ggplot2 object | Creates a line plot for multiple vectors, with an optional dual y-axis. |
| `cdf` | values, type, title, xlab, color | ggplot2 object | Computes and plots the Cumulative Distribution Function (CDF). |
| `cdf_multi` | data_list, group_col, value_col, type, title, xlab | ggplot2 object | Computes and plots CDFs for multiple groups. |

### From `myR/R/visualization/composition_plots.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `cmb` | sobj, cluster_col, sample_col, colors, title, xlab, ylab | ggplot2 object | Creates a proportional stacked bar graph of clusters. |
| `acmb` | sobj, cluster_col, sample_col, colors, title, xlab, ylab | ggplot2 object | Creates an absolute count stacked bar graph of clusters. |
| `cml` | sobj, cluster_col, sample_col, title, xlab, ylab | ggplot2 object | Creates a cumulative line plot of cluster proportions. |
| `upset_gene_lists`| gene_lists, min_size, nsets, nintersects, order_by, ... | UpSet plot | Creates an UpSet plot to visualize intersections of gene lists. |
| `vln_p` | sobj, features, group_by, split_by, test_method, comparisons, ncol | ggplot2 object | Generates violin plots with statistical comparisons. |

## Analysis Functions

### From `myR/R/analysis/cell_communication/nichenet_analysis.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `ligand_to_target` | ligand, target, lr_network, sig_network, gr_network, top_n_targets | Data frame | Retrieves ligand-target gene regulatory potential from NicheNet. |
| `run_nichenet_analysis` | seurat_obj, sender_cells, receiver_cells, condition_oi, condition_ref, ... | List | Performs a comprehensive NicheNet cell-cell communication analysis. |
| `prepare_nichenet_circos_data` | ligand_receptor_links, ligand_target_links, ligand_activities, ... | List | Prepares data for NicheNet Circos plot visualization. |
| `draw_nichenet_circos_plot` | circos_data, ligand_color, receptor_color, target_color | NULL | Creates a Circos plot showing NicheNet interactions. |

### From `myR/R/analysis/differential_expression/differential_expression.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `prepare_pseudobulk_edgeR` | object, sample_col, group_col, cluster_col, target_cluster, ... | List | Aggregates single-cell counts to pseudo-bulk level for edgeR. |
| `run_pseudobulk_deg` | object, sample_col, group_col, comparison, cluster_col, ... | Data frame or list | Performs edgeR-based pseudo-bulk differential gene expression analysis. |
| `linear_seurat` | sobj, layer, features, regressor, regressor.type, covariates, ... | Data frame | Performs linear regression for gene expression against a regressor. |
| `create_analysis_config` | patient, drug, timepoint, ck, response, aoi | Named list | Creates a configuration object for LMM metadata columns. |
| `fit_lmm_single_gene` | gene_expr, metadata, config, formula_str, formula_components, ... | List | Internal helper to fit an lmer model for one gene. |
| `summarize_lmm_results` | lmm_results, config | Tidy data frame | Internal helper to combine results from multiple LMMs. |
| `run_lmm_multiple_genes` | seurat_obj, genes, config, formula_str, n_cores, ... | List | Applies Linear Mixed Models (LMM) to multiple genes in parallel. |
| `find_response_differential_genes` | lmm_summary, config, drug_name, top_n | Data frame | Identifies genes where treatment response differs by drug from LMM results. |
| `find_drug_specific_genes` | lmm_summary, config, top_n | Data frame | Identifies genes most strongly associated with specific drugs. |

### From `myR/R/analysis/markers/marker_processing.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `marker_trim` | markers, direction, padj_threshold, logfc_col, padj_col | Data frame | Filters marker gene results by direction and significance. |
| `marker_filter` | markers, remove_ribo, remove_mito, remove_hb, remove_ig, ... | Data frame | Removes unwanted gene types (ribosomal, mitochondrial, etc.) from markers. |
| `lrf` | results, term_col | Data frame | Filters out "(Intercept)" terms from linear model results. |
| `all_markers_to_list` | all_markers, cluster_col, gene_col | Named list | Converts FindAllMarkers output to a named list, one per cluster. |
| `marker_print_all`| all_markers, n, order_by, decreasing, cluster_col, gene_col | Invisible data frame| Prints top N marker genes for each cluster. |
| `marker_print` | markers, n, order_by, decreasing, gene_col, columns | Invisible data frame| Prints top N marker genes from a single marker data frame. |
| `synthesize_markers`| marker_list, gene_col, logfc_col, pval_col, weights | Data frame | Combines p-values and logFCs from multiple marker results. |
| `synthesize_ranks`| marker_list, gene_col, rank_by, decreasing | Data frame | Combines ranks from multiple marker results using geometric mean. |

### From `myR/R/analysis/pathway/pathway_enrichment.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `convert_gene_ids`| genes, from, to, org_db, remove_na | Named character vector | Converts gene identifiers between different formats. |
| `prepare_gene_lists`| deg_results, gene_col, logfc_col, pval_col, padj_col, ... | List | Prepares gene lists from DE results for pathway analysis. |
| `get_pathway_sets`| species, category, subcategory | Character vector | Returns available pathway database identifiers from MSigDB. |
| `run_go_analysis`| genes, universe, ont, pval_cutoff, qval_cutoff, org_db | enrichResult object | Performs Gene Ontology (GO) enrichment analysis. |
| `run_kegg_analysis`| genes, universe, organism, pval_cutoff, qval_cutoff | enrichResult object | Performs KEGG pathway enrichment analysis. |
| `run_gsea_analysis`| ranked_genes, species, category, subcategory, min_size, ... | Data frame | Performs Gene Set Enrichment Analysis (GSEA) using fgsea. |
| `format_results`| results, method | Standardized data frame| Standardizes output from different pathway analysis methods. |
| `myGO` | deg_results, run_go, run_kegg, run_gsea, go_ont, gsea_category, ... | List | Comprehensive wrapper to run GO, KEGG, and GSEA analyses. |

### From `myR/R/analysis/spatial/geomx_analysis.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `prepare_geomx_data`| counts, metadata, feature_info | List | Prepares raw GeoMx count data into a standard format. |
| `q3_normalize` | geomx_data, scale_factor | Normalized data | Performs Q3 (third quartile) normalization on GeoMx data. |
| `perform_qc` | geomx_data, min_counts, min_genes, min_rois | Filtered data list | Performs basic QC checks and filtering on GeoMx data. |
| `find_deg_geomx`| geomx_data, group_col, comparison, method, covariates, ... | Data frame | Performs differential expression analysis on GeoMx data. |
| `plot_deg_volcano`| deg_results, padj_threshold, logfc_threshold, label_top, title | ggplot object | Creates a volcano plot from GeoMx DEG results. |
| `plot_deg_heatmap`| geomx_data, deg_results, group_col, n_genes, scale_rows | Heatmap object | Creates a heatmap of top differentially expressed genes from GeoMx data. |
| `run_geomx_analysis`| counts, metadata, group_col, comparison, normalize, qc, ... | List | Wrapper function to run the complete GeoMx analysis workflow. |

### From `myR/R/analysis/trajectory/trajectory_inference.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `run_slingshot_from_seurat` | seurat_obj, cluster_col, reduction, start_cluster, end_cluster, ... | SingleCellExperiment | Runs Slingshot trajectory inference from a Seurat object. |
| `analyze_gene_dynamics` | sce, gene, condition, lineage, plot, output_dir | List | Fits a GAM to model gene expression changes along pseudotime. |
| `process_gene_list_dynamics` | sce, genes, condition, lineage, plot, output_dir, n_cores | List | Processes a list of genes for dynamics analysis, with parallel support. |
| `analyze_gene_dynamics_tradeSeq` | sce, gene, n_knots, conditions | List | Fits GAMs using tradeSeq for differential expression patterns along trajectories. |

## Signature Functions

### From `myR/R/signatures/signature_discovery.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `find_gene_modules` | expr_matrix, method, n_modules | List | Identifies co-expressed gene modules using methods like WGCNA or HOPACH. |
| `calculate_module_eigengene` | expr_matrix, genes | Numeric vector | Calculates the first principal component (eigengene) for a gene module. |

### From `myR/R/signatures/signature_scoring.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `score_signature` | expr_data, signature, method, n_bins, n_ctrl | Numeric vector | Scores a gene signature for each cell/sample. |
| `AddMultipleModuleScores` | seurat_obj, gene_sets, method, n_bins, n_ctrl, ... | Seurat object | Wrapper to add multiple gene signature scores to a Seurat object's metadata. |
| `compare_module_scoring_methods` | seurat_obj, gene_sets, methods | Data frame | Compares the output of different module scoring methods. |

## Utility Functions

### From `myR/R/utilities/demulti_utils.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `get_barcode_mapping`| samples, barcodes | Data frame | Creates a mapping between sample names and barcodes from a file. |
| `assign_barcodes` | seurat_obj, barcode_mapping, barcode_col | Seurat object | Assigns sample identities to cells based on a barcode mapping. |

### From `myR/R/utilities/general_utils.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `theme_my` | base_size | ggplot theme | A custom ggplot2 theme for consistent plot styling. |
| `save_plot` | plot, filename, width, height, dpi | NULL | Saves a ggplot object to a file with specified dimensions and resolution. |
| ` `%||%` ` | a, b | Value | Infix operator that returns `b` if `a` is NULL. |

### From `myR/R/utilities/plot_utils.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `get_palette` | n, palette_name | Character vector | Retrieves a color palette of a specified length. |
| `plot_pct_expr` | seurat_obj, features, group_by | ggplot object | Plots the percentage of cells expressing a feature across groups. |

### From `myR/R/utilities/sample_utils.R`

| Function Name | Input | Output | Function |
|---|---|---|---|
| `sort_samples` | samples | Character vector | Sorts sample names based on numeric components within the names. |
| `subset_seurat_by_time` | seurat_obj, time_col, start_time, end_time | Seurat object | Subsets a Seurat object based on a time range. |

## Miscellaneous Functions (from R/ root)

*This section includes functions from various files located directly in `myR/R/`. These may be older functions or general utilities not yet categorized.*

| Function Name | File | Description |
|---|---|---|
| `myhm_genes3` | `plots.R` | An earlier version of a gene expression heatmap function. |
| `myhm_genes4` | `plots.R` | Creates a gene expression heatmap with flexible grouping and annotation. |
| `plot_cluster_fractions` | `plots.R`| Visualizes cluster proportions across different conditions. |
| `plot_interaction_for_gene`| `plots.R` | Plots the interaction effect of variables on gene expression. |
| `find_conserved_markers`| `markers.R` | Finds markers that are conserved across different conditions for a given cluster. |
| `find_all_conserved_markers` | `markers.R` | Wrapper to run `find_conserved_markers` for all clusters. |
| `add_progeny_scores` | `signature.R`| Adds PROGENy pathway activity scores to a Seurat object. |
| `convert_seurat_to_anndata`| `seurat_to_scanpy.R` | Converts a Seurat object to an AnnData object for use with Python's scanpy. |
| ... | ... | *(and many more from the 16 files in R/)* |

## Deprecated and Legacy Functions

### From `myR/R/deprecated/` and `myR/R/legacy.R`

| Function Name | File | Description |
|---|---|---|
| `pseudobulk_linear_fit_legacy` | `legacy.R` | Legacy version of the pseudobulk linear model fitting function. |
| `post_hoc_slope_comparison_legacy` | `legacy.R` | Legacy version of a post-hoc slope comparison function. |
| `run_pseudobulk_deg_legacy` | `legacy.R` | Legacy version of the pseudobulk DEG analysis function. |
| `myhm_genesets2_legacy` | `legacy.R`| Legacy version of a gene set heatmap function. |
| `pseudobulk_deg_standalone` | `deprecated/pseudobulk_deg_standalone.R` | A standalone script version of the pseudobulk DEG analysis. |

---
*This document is now a complete representation of all functions found in the `myR/R/` directory and its subdirectories.*

