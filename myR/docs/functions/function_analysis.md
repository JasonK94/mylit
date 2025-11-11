# myR Function Analysis

## Directory: `core`

### File: `myR/R/core/data_preparation.R`

#### Function: `.get_feature_vector`
- **Description:** Internal helper to extract a feature (gene or metadata column) from a Seurat object.
- **Parameters:**
  - `object`: Seurat object
  - `feature`: Feature name (gene or metadata column)
  - `assay`: Assay to use (default: DefaultAssay)
  - `slot`: Slot to use (default: "data")
  - `cells`: Subset of cells to extract (default: all cells)
- **Returns:** Numeric vector of feature values

#### Function: `get_feature_vec`
- **Description:** Extracts a feature (gene or metadata column) from a Seurat object.
- **Parameters:**
  - `object`: Seurat object
  - `feature`: Feature name (gene or metadata column)
  - `assay`: Assay to use (default: DefaultAssay)
  - `slot`: Slot to use (default: "data")
  - `cells`: Subset of cells to extract (default: all cells)
- **Returns:** Numeric vector of feature values

#### Function: `get_feature_matrix`
- **Description:** Extracts multiple features from a Seurat object and returns as a matrix.
- **Parameters:**
  - `object`: Seurat object
  - `features`: Character vector of feature names
  - `assay`: Assay to use (default: DefaultAssay)
  - `slot`: Slot to use (default: "data")
  - `cells`: Subset of cells to extract (default: all cells)
- **Returns:** Matrix with features as rows and cells as columns

#### Function: `prepare_metadata_table`
- **Description:** Extracts and optionally filters metadata from a Seurat object.
- **Parameters:**
  - `object`: Seurat object
  - `columns`: Specific columns to extract (NULL = all columns)
  - `cells`: Subset of cells (NULL = all cells)
  - `drop_na`: Whether to drop rows with NA values (default: FALSE)
- **Returns:** Data frame of metadata

#### Function: `aggregate_expression_by_group`
- **Description:** Aggregates gene expression across cells within groups (e.g., clusters, samples).
- **Parameters:**
  - `object`: Seurat object
  - `features`: Features (genes) to aggregate
  - `group_by`: Metadata column to group by
  - `method`: Aggregation method: "mean", "median", "sum" (default: "mean")
  - `assay`: Assay to use (default: DefaultAssay)
  - `slot`: Slot to use (default: "data")
- **Returns:** Matrix with features as rows and groups as columns

#### Function: `prepare_count_matrix`
- **Description:** Extracts a count matrix from a Seurat object with optional filtering.
- **Parameters:**
  - `object`: Seurat object
  - `assay`: Assay to use (default: "RNA")
  - `slot`: Slot to use (default: "counts")
  - `features`: Features to include (NULL = all)
  - `cells`: Cells to include (NULL = all)
  - `min_cells`: Minimum number of cells expressing a feature (default: 0)
  - `min_features`: Minimum number of features expressed in a cell (default: 0)
- **Returns:** Sparse or dense matrix of counts

#### Function: `convert_to_long_format`
- **Description:** Converts expression data to long (tidy) format for ggplot2.
- **Parameters:**
  - `object`: Seurat object
  - `features`: Features to include
  - `metadata_cols`: Metadata columns to include
  - `assay`: Assay to use (default: DefaultAssay)
  - `slot`: Slot to use (default: "data")
  - `cells`: Cells to include (NULL = all)
- **Returns:** Data frame in long format with columns: cell, feature, expression, and metadata

#### Function: `check_data_quality`
- **Description:** Performs basic quality checks on expression data.
- **Parameters:**
  - `data`: Matrix or data frame of expression data
  - `check_finite`: Check for infinite values (default: TRUE)
  - `check_na`: Check for NA values (default: TRUE)
  - `check_negative`: Check for negative values (default: TRUE)
- **Returns:** List with logical flags and counts of issues found

### File: `myR/R/core/validation.R`

#### Function: `validate_seurat`
- **Description:** Checks if the input is a valid Seurat object with expected components.
- **Parameters:**
  - `obj`: Object to validate
  - `assay`: Optional assay name to check for
  - `reduction`: Optional reduction name to check for
  - `min_cells`: Minimum number of cells required
  - `min_features`: Minimum number of features required
- **Returns:** TRUE if valid, stops with error message otherwise

#### Function: `validate_metadata_column`
- **Description:** Checks if a metadata column exists and optionally validates its type.
- **Parameters:**
  - `obj`: Seurat object or data frame
  - `column_name`: Column name to validate
  - `required_type`: Optional required type ("numeric", "factor", "character")
  - `allow_na`: Whether NA values are allowed
- **Returns:** TRUE if valid, stops with error message otherwise

#### Function: `validate_genes`
- **Description:** Checks if genes exist in the object and optionally filters to valid genes.
- **Parameters:**
  - `obj`: Seurat object or character vector of available genes
  - `genes`: Character vector of genes to validate
  - `min_present`: Minimum number of genes that must be present (default: all)
  - `assay`: Assay to check genes in (for Seurat objects)
  - `warn_missing`: Whether to warn about missing genes
- **Returns:** Character vector of valid genes present in the object

#### Function: `validate_numeric_range`
- **Description:** Checks if a numeric parameter is within acceptable range.
- **Parameters:**
  - `value`: Numeric value to validate
  - `param_name`: Name of parameter (for error messages)
  - `min`: Minimum allowed value (inclusive)
  - `max`: Maximum allowed value (inclusive)
  - `allow_na`: Whether NA is allowed
- **Returns:** TRUE if valid, stops with error message otherwise

#### Function: `validate_choice`
- **Description:** Validates that a value is one of allowed choices (like match.arg but more informative).
- **Parameters:**
  - `value`: Value to validate
  - `param_name`: Name of parameter (for error messages)
  - `choices`: Vector of allowed values
  - `multiple`: Whether multiple choices are allowed
- **Returns:** The validated value (or values if multiple=TRUE)

#### Function: `validate_path`
- **Description:** Checks if a file path exists and optionally validates extension.
- **Parameters:**
  - `path`: File path to validate
  - `must_exist`: Whether file must already exist
  - `extensions`: Optional vector of allowed extensions (e.g., c("csv", "txt"))
  - `type`: Type of path ("file" or "directory")
- **Returns:** Normalized path if valid, stops with error otherwise

#### Function: `create_error_message`
- **Description:** Helper to create consistent, informative error messages.
- **Parameters:**
  - `context`: Context where error occurred (e.g., function name)
  - `message`: Main error message
  - `suggestion`: Optional suggestion for fixing the error
- **Returns:** Formatted error message

#### Function: `check_packages`
- **Description:** Checks if required packages are installed and optionally loads them.
- **Parameters:**
  - `packages`: Character vector of package names
  - `load`: Whether to load the packages (default: FALSE)
- **Returns:** TRUE if all packages available, stops with error otherwise

## Directory: `analysis`

### File: `myR/R/analysis/cell_communication/nichenet_analysis.R`

#### Function: `ligand_to_target`
- **Description:** Retrieves ligand-target gene regulatory potential from NicheNet databases.
- **Parameters:**
  - `ligand`: Character vector of ligand genes
  - `target`: Character vector of target genes
  - `lr_network`: Ligand-receptor network (from NicheNet)
  - `sig_network`: Signaling network (from NicheNet)
  - `gr_network`: Gene regulatory network (from NicheNet)
  - `top_n_targets`: Number of top targets to return per ligand (default: 250)
- **Returns:** Data frame with ligand-target regulatory potential scores

#### Function: `run_nichenet_analysis`
- **Description:** Performs a comprehensive NicheNet cell-cell communication analysis workflow, including ligand activity prediction, receptor inference, and visualization.
- **Parameters:**
  - `seurat_obj`: Seurat object
  - `sender_cells`: Vector of sender cell identities or logical vector
  - `receiver_cells`: Vector of receiver cell identities or logical vector
  - `condition_oi`: Condition of interest for receiver DE analysis
  - `condition_ref`: Reference condition for receiver DE analysis
  - `ident_col`: Identity column in metadata (default: "seurat_clusters")
  - `condition_col`: Condition column in metadata (required for DE)
  - `lr_network`: Ligand-receptor network (from NicheNet)
  - `sig_network`: Signaling network (from NicheNet)
  - `gr_network`: Gene regulatory network (from NicheNet)
  - `ligand_target_matrix`: Pre-computed ligand-target matrix (optional)
  - `expressed_pct`: Expression percentage threshold (default: 0.10)
  - `top_n_ligands`: Number of top ligands to analyze (default: 20)
  - `top_n_targets`: Number of top targets per ligand (default: 200)
  - `plot_circos`: Generate Circos plot (default: FALSE)
  - `output_dir`: Directory for saving plots (default: NULL)
- **Returns:** List containing: ligand_activities, best_ligands, ligand_target_links, ligand_receptor_links, plots, warnings

#### Function: `prepare_nichenet_circos_data`
- **Description:** Prepares ligand-receptor and ligand-target data for Circos visualization.
- **Parameters:**
  - `ligand_receptor_links`: Data frame with ligand-receptor pairs
  - `ligand_target_links`: Data frame with ligand-target regulatory links
  - `ligand_activities`: Data frame with ligand activity scores
  - `top_n_ligands`: Number of top ligands to include (default: 10)
  - `top_n_targets`: Number of top targets per ligand (default: 20)
- **Returns:** List with formatted data for Circos plot

#### Function: `draw_nichenet_circos_plot`
- **Description:** Creates a Circos plot showing ligand-receptor and ligand-target interactions.
- **Parameters:**
  - `circos_data`: List from prepare_nichenet_circos_data
  - `ligand_color`: Color for ligands (default: "#E41A1C")
  - `receptor_color`: Color for receptors (default: "#377EB8")
  - `target_color`: Color for targets (default: "#4DAF4A")
- **Returns:** NULL (plot is drawn to current device)

### File: `myR/R/analysis/differential_expression/differential_expression.R`

#### Function: `prepare_pseudobulk_edgeR`
- **Description:** Aggregates single-cell counts to pseudo-bulk level for differential expression analysis.
- **Parameters:**
  - `object`: Seurat object
  - `sample_col`: Metadata column identifying samples/replicates
  - `group_col`: Metadata column for grouping (e.g., condition, treatment)
  - `cluster_col`: Optional cluster column for per-cluster analysis (default: NULL)
  - `target_cluster`: If cluster_col is specified, analyze only this cluster (default: NULL)
  - `assay`: Assay to use (default: "RNA")
  - `slot`: Slot to use (default: "counts")
  - `min_cells`: Minimum cells per sample to include (default: 10)
  - `min_counts`: Minimum total counts per gene (default: 10)
- **Returns:** List containing counts, metadata, and cluster info

#### Function: `run_pseudobulk_deg`
- **Description:** Performs edgeR-based pseudo-bulk differential gene expression analysis.
- **Parameters:**
  - `object`: Seurat object or prepared pseudobulk list
  - `sample_col`: Sample identifier column
  - `group_col`: Group comparison column
  - `comparison`: Two-element vector for comparison (e.g., c("Treatment", "Control"))
  - `cluster_col`: Optional cluster column
  - `target_cluster`: Specific cluster to analyze
  - `mode`: Analysis mode: "overall", "per_cluster", "specific_cluster"
  - `assay`: Assay to use (default: "RNA")
  - `slot`: Slot to use (default: "counts")
  - `min_cells`: Minimum cells per sample (default: 10)
  - `min_counts`: Minimum total counts per gene (default: 10)
  - `fdr_threshold`: FDR threshold (default: 0.05)
  - `logfc_threshold`: Log fold change threshold (default: 0)
- **Returns:** Data frame with DE results or list of data frames if per_cluster

#### Function: `.run_edger_analysis`
- **Description:** Internal edgeR Analysis.
- **Parameters:**
  - `counts`: counts
  - `metadata`: metadata
  - `group_col`: group_col
  - `comparison`: comparison
  - `fdr_threshold`: fdr_threshold
  - `logfc_threshold`: logfc_threshold
- **Returns:** Data frame with DE results

#### Function: `linear_seurat`
- **Description:** Performs linear regression analysis for gene expression against a regressor with support for continuous, categorical, and ordinal predictors.
- **Parameters:**
  - `sobj`: Seurat object
  - `layer`: Expression layer: "counts", "data", "scale.data"
  - `features`: Features to test (default: "all")
  - `regressor`: Regressor variable name in metadata
  - `regressor.type`: Type: "continuous", "categorical", "ordinal"
  - `reference.level`: Reference level for categorical
  - `ordinal.method`: Method for ordinal: "linear", "polynomial", "spline"
  - `link.function`: Link function: "linear", "poisson", "negative.binomial"
  - `effect`: Effect type: "fixed", "random"
  - `covariates`: Covariate column names
  - `min.cells`: Minimum cells expressing gene (default: 10)
  - `return.full`: Return full results including Seurat object
- **Returns:** Data frame with regression results

#### Function: `create_analysis_config`
- **Description:** Creates a configuration object to manage metadata column names consistently throughout complex experimental designs (e.g., GeoMx with patient, drug, timepoint).
- **Parameters:**
  - `patient`: Column name for patient ID
  - `drug`: Column name for drug/treatment
  - `timepoint`: Column name for timepoint (e.g., pre/post)
  - `ck`: Column name for stratification variable (e.g., CK status)
  - `response`: Column name for treatment response
  - `aoi`: Column name for AOI (Area of Interest) ID
- **Returns:** Named list containing column names

#### Function: `fit_lmm_single_gene`
- **Description:** Internal helper to fit an lmer model for one gene with complex experimental design.
- **Parameters:**
  - `gene_expr`: Numeric vector of expression values
  - `metadata`: Metadata data frame
  - `config`: Analysis configuration from create_analysis_config()
  - `formula_str`: Optional custom formula string
  - `formula_components`: List with fixed, interactions, and random components
  - `use_config_names`: Whether to map generic names to config names
- **Returns:** List with model, effects, anova, convergence status

#### Function: `summarize_lmm_results`
- **Description:** Internal helper to combine results from multiple LMMs and calculate adjusted p-values.
- **Parameters:**
  - `lmm_results`: List of results from fit_lmm_single_gene
  - `config`: Analysis configuration
- **Returns:** Tidy data frame summarizing all model effects

#### Function: `run_lmm_multiple_genes`
- **Description:** Applies LMM to multiple genes in parallel. Main workhorse for LMM analysis with complex experimental designs (e.g., patient, drug, timepoint, response).
- **Parameters:**
  - `seurat_obj`: Seurat object
  - `genes`: Character vector of gene names to analyze
  - `config`: Analysis configuration from create_analysis_config()
  - `formula_str`: Optional custom formula string
  - `formula_components`: List specifying model formula components
  - `use_config_names`: Whether to use config names in formula
  - `n_cores`: Number of CPU cores for parallel processing
  - `verbose`: Whether to print progress messages
- **Returns:** List containing: raw_results, summary, converged_genes, total_genes

#### Function: `find_response_differential_genes`
- **Description:** Identifies genes where treatment response differs by drug from LMM results.
- **Parameters:**
  - `lmm_summary`: Summary data frame from run_lmm_multiple_genes
  - `config`: Analysis configuration
  - `drug_name`: Optional specific drug name to focus on
  - `top_n`: Number of top genes to return
- **Returns:** Data frame of top genes ranked by effect size

#### Function: `find_drug_specific_genes`
- **Description:** Identifies genes most strongly associated with specific drugs.
- **Parameters:**
  - `lmm_summary`: Summary data frame from run_lmm_multiple_genes
  - `config`: Analysis configuration
  - `top_n`: Number of top genes to return
- **Returns:** Data frame of top genes ranked by effect size

### File: `myR/R/analysis/markers/marker_processing.R`

#### Function: `marker_trim`
- **Description:** Filters marker gene results based on log fold change direction and adjusted p-value threshold.
- **Parameters:**
  - `markers`: Data frame from FindMarkers/FindAllMarkers
  - `direction`: Filter direction: "up" (positive logFC), "down" (negative logFC), or "both" (default: "both")
  - `padj_threshold`: Adjusted p-value threshold (default: 0.05)
  - `logfc_col`: Name of log fold change column (default: "avg_log2FC")
  - `padj_col`: Name of adjusted p-value column (default: "p_val_adj")
- **Returns:** Filtered data frame

#### Function: `marker_filter`
- **Description:** Removes ribosomal, mitochondrial, hemoglobin, and other unwanted genes from marker results.
- **Parameters:**
  - `markers`: Data frame from FindMarkers/FindAllMarkers
  - `remove_ribo`: Remove ribosomal genes (^RP[SL]) (default: TRUE)
  - `remove_mito`: Remove mitochondrial genes (^MT-) (default: TRUE)
  - `remove_hb`: Remove hemoglobin genes (^HB[AB]) (default: TRUE)
  - `remove_ig`: Remove immunoglobulin genes (^IG[HKL]) (default: FALSE)
  - `remove_tr`: Remove T-cell receptor genes (^TR[ABGD]) (default: FALSE)
  - `custom_patterns`: Additional regex patterns to remove (default: NULL)
  - `gene_col`: Name of gene column (default: "gene" or rownames if not present)
- **Returns:** Filtered data frame

#### Function: `lrf`
- **Description:** Filters out "(Intercept)" terms from linear mixed model results.
- **Parameters:**
  - `results`: Data frame with model results
  - `term_col`: Name of column containing term names (default: "term")
- **Returns:** Filtered data frame

#### Function: `all_markers_to_list`
- **Description:** Converts FindAllMarkers output to a named list, with one element per cluster.
- **Parameters:**
  - `all_markers`: Data frame from FindAllMarkers
  - `cluster_col`: Name of cluster column (default: "cluster")
  - `gene_col`: Name of gene column (default: "gene" or rownames)
- **Returns:** Named list of marker data frames, one per cluster

#### Function: `marker_print_all`
- **Description:** Prints top N marker genes for each cluster from FindAllMarkers output.
- **Parameters:**
  - `all_markers`: Data frame from FindAllMarkers
  - `n`: Number of top markers to print per cluster (default: 10)
  - `order_by`: Column to order by (default: "avg_log2FC")
  - `decreasing`: Sort in decreasing order (default: TRUE)
  - `cluster_col`: Name of cluster column (default: "cluster")
  - `gene_col`: Name of gene column (default: "gene" or rownames)
- **Returns:** Invisibly returns the input data frame

#### Function: `marker_print`
- **Description:** Prints top N marker genes from a single marker data frame.
- **Parameters:**
  - `markers`: Data frame from FindMarkers
  - `n`: Number of top markers to print (default: 10)
  - `order_by`: Column to order by (default: "avg_log2FC")
  - `decreasing`: Sort in decreasing order (default: TRUE)
  - `gene_col`: Name of gene column (default: "gene" or rownames)
  - `columns`: Columns to display (default: c("avg_log2FC", "p_val_adj"))
- **Returns:** Invisibly returns the input data frame

#### Function: `synthesize_markers`
- **Description:** Combines p-values and log fold changes from multiple FindMarkers results using Fisher's method for p-values and weighted mean for log fold changes.
- **Parameters:**
  - `marker_list`: Named list of marker data frames
  - `gene_col`: Name of gene column (default: "gene" or rownames)
  - `logfc_col`: Name of log fold change column (default: "avg_log2FC")
  - `pval_col`: Name of p-value column (default: "p_val")
  - `weights`: Optional weights for each dataset (default: equal weights)
- **Returns:** Data frame with synthesized results

#### Function: `synthesize_ranks`
- **Description:** Combines ranks from multiple FindMarkers results using geometric mean of ranks.
- **Parameters:**
  - `marker_list`: Named list of marker data frames
  - `gene_col`: Name of gene column (default: "gene" or rownames)
  - `rank_by`: Column to rank by (default: "avg_log2FC")
  - `decreasing`: Rank in decreasing order (default: TRUE)
- **Returns:** Data frame with synthesized ranks

### File: `myR/R/analysis/pathway/pathway_enrichment.R`

#### Function: `convert_gene_ids`
- **Description:** Converts gene identifiers between different formats (e.g., SYMBOL, ENTREZID, ENSEMBL).
- **Parameters:**
  - `genes`: Character vector of gene identifiers
  - `from`: Source ID type (default: "SYMBOL")
  - `to`: Target ID type (default: "ENTREZID")
  - `org_db`: Organism database (default: org.Hs.eg.db for human)
  - `remove_na`: Remove genes that couldn't be converted (default: TRUE)
- **Returns:** Named character vector of converted IDs (names are original IDs)

#### Function: `prepare_gene_lists`
- **Description:** Prepares gene lists from differential expression results for pathway analysis.
- **Parameters:**
  - `deg_results`: Data frame with DE results (must contain logFC and p-value columns)
  - `gene_col`: Name of gene column (default: "gene" or rownames)
  - `logfc_col`: Name of log fold change column (default: "logFC" or "avg_log2FC")
  - `pval_col`: Name of p-value column (default: "pvalue" or "p_val")
  - `padj_col`: Name of adjusted p-value column (default: "padj" or "p_val_adj")
  - `padj_threshold`: Adjusted p-value threshold (default: 0.05)
  - `logfc_threshold`: Log fold change threshold (default: 0)
  - `convert_ids`: Convert gene IDs to ENTREZID (default: TRUE)
- **Returns:** List containing: ranked_genes, up_genes, down_genes, sig_genes, background_genes

#### Function: `get_pathway_sets`
- **Description:** Returns available pathway database identifiers for GSEA.
- **Parameters:**
  - `species`: Species code (default: "Homo sapiens")
  - `category`: MSigDB category (e.g., "H", "C2", "C5")
  - `subcategory`: MSigDB subcategory (e.g., "CP:KEGG", "GO:BP")
- **Returns:** Character vector of available gene set names

#### Function: `run_go_analysis`
- **Description:** Performs Gene Ontology (GO) enrichment analysis using clusterProfiler.
- **Parameters:**
  - `genes`: Character vector of gene IDs (ENTREZID format)
  - `universe`: Background gene universe (optional)
  - `ont`: GO ontology: "BP", "MF", "CC", or "ALL" (default: "BP")
  - `pval_cutoff`: P-value cutoff (default: 0.05)
  - `qval_cutoff`: Q-value cutoff (default: 0.05)
  - `org_db`: Organism database (default: org.Hs.eg.db)
- **Returns:** enrichResult object from clusterProfiler

#### Function: `run_kegg_analysis`
- **Description:** Performs KEGG pathway enrichment analysis using clusterProfiler.
- **Parameters:**
  - `genes`: Character vector of gene IDs (ENTREZID format)
  - `universe`: Background gene universe (optional)
  - `organism`: KEGG organism code (default: "hsa" for human)
  - `pval_cutoff`: P-value cutoff (default: 0.05)
  - `qval_cutoff`: Q-value cutoff (default: 0.05)
- **Returns:** enrichResult object from clusterProfiler

#### Function: `run_gsea_analysis`
- **Description:** Performs GSEA using fgsea with MSigDB gene sets.
- **Parameters:**
  - `ranked_genes`: Named numeric vector of genes ranked by score (e.g., logFC)
  - `species`: Species for MSigDB (default: "Homo sapiens")
  - `category`: MSigDB category (default: "H" for Hallmark)
  - `subcategory`: MSigDB subcategory (optional)
  - `min_size`: Minimum gene set size (default: 15)
  - `max_size`: Maximum gene set size (default: 500)
  - `nperm`: Number of permutations (default: 10000)
- **Returns:** Data frame with GSEA results

#### Function: `format_results`
- **Description:** Standardizes output from different pathway analysis methods.
- **Parameters:**
  - `results`: Results object from GO, KEGG, or GSEA analysis
  - `method`: Method used: "GO", "KEGG", or "GSEA"
- **Returns:** Standardized data frame

#### Function: `myGO`
- **Description:** Runs GO, KEGG, and GSEA analyses and returns combined results.
- **Parameters:**
  - `deg_results`: Data frame with differential expression results
  - `run_go`: Perform GO analysis (default: TRUE)
  - `run_kegg`: Perform KEGG analysis (default: TRUE)
  - `run_gsea`: Perform GSEA (default: TRUE)
  - `go_ont`: GO ontology (default: "BP")
  - `gsea_category`: MSigDB category for GSEA (default: "H")
  - `padj_threshold`: Adjusted p-value threshold for gene selection (default: 0.05)
  - `logfc_threshold`: Log fold change threshold (default: 0)
- **Returns:** List containing results from each analysis method

### File: `myR/R/analysis/spatial/geomx_analysis.R`

#### Function: `prepare_geomx_data`
- **Description:** Prepares raw GeoMx count data for analysis by organizing into a standard format.
- **Parameters:**
  - `counts`: Count matrix (genes x ROIs)
  - `metadata`: Data frame with ROI-level metadata
  - `feature_info`: Optional data frame with gene/probe information
- **Returns:** List containing: counts, metadata, features

#### Function: `q3_normalize`
- **Description:** Performs Q3 (third quartile) normalization on GeoMx count data.
- **Parameters:**
  - `geomx_data`: List from prepare_geomx_data or count matrix
  - `scale_factor`: Scaling factor after normalization (default: 1000)
- **Returns:** Normalized data in same format as input

#### Function: `perform_qc`
- **Description:** Performs basic QC checks and filtering on GeoMx data.
- **Parameters:**
  - `geomx_data`: List from prepare_geomx_data
  - `min_counts`: Minimum total counts per ROI (default: 1000)
  - `min_genes`: Minimum number of detected genes per ROI (default: 100)
  - `min_rois`: Minimum number of ROIs expressing a gene (default: 2)
- **Returns:** Filtered GeoMx data list with QC metrics added to metadata

#### Function: `find_deg_geomx`
- **Description:** Performs differential expression analysis on GeoMx data using limma or Wilcoxon test.
- **Parameters:**
  - `geomx_data`: List from prepare_geomx_data (normalized)
  - `group_col`: Metadata column for grouping
  - `comparison`: Comparison to make (e.g., c("Treatment", "Control"))
  - `method`: Method to use: "limma" or "wilcox" (default: "limma")
  - `covariates`: Optional covariate columns to include in limma model
  - `log_transform`: Apply log2(x+1) transformation (default: TRUE for limma)
- **Returns:** Data frame with differential expression results

#### Function: `.geomx_limma`
- **Description:** Internal limma Analysis.
- **Parameters:**
  - `expr_data`: Expression data
  - `metadata`: Metadata
  - `group_col`: Group column
  - `comparison`: Comparison
  - `covariates`: Covariates
- **Returns:** Data frame with limma results

#### Function: `.geomx_wilcox`
- **Description:** Internal Wilcoxon Test.
- **Parameters:**
  - `expr_data`: Expression data
  - `metadata`: Metadata
  - `group_col`: Group column
  - `comparison`: Comparison
- **Returns:** Data frame with Wilcoxon results

#### Function: `plot_deg_volcano`
- **Description:** Creates a volcano plot from differential expression results.
- **Parameters:**
  - `deg_results`: Data frame from find_deg_geomx
  - `padj_threshold`: Adjusted p-value threshold (default: 0.05)
  - `logfc_threshold`: Log fold change threshold (default: 1)
  - `label_top`: Number of top genes to label (default: 10)
  - `title`: Plot title
- **Returns:** ggplot object

#### Function: `plot_deg_heatmap`
- **Description:** Creates a heatmap of top differentially expressed genes.
- **Parameters:**
  - `geomx_data`: GeoMx data list (normalized)
  - `deg_results`: Data frame from find_deg_geomx
  - `group_col`: Metadata column for grouping
  - `n_genes`: Number of top genes to plot (default: 50)
  - `scale_rows`: Z-score scale rows (default: TRUE)
- **Returns:** Heatmap (ComplexHeatmap or pheatmap object)

#### Function: `run_geomx_analysis`
- **Description:** Wrapper function to run the complete GeoMx analysis workflow.
- **Parameters:**
  - `raw_data`: Raw data matrix.
  - `metadata_SP`: Metadata data frame.
  - `group_var`: Column for differential expression.
  - `group1`: First group for comparison.
  - `group2`: Second group for comparison (optional).
  - `qc_filter`: Perform QC filtering (default: TRUE).
  - `min_genes`: Minimum number of detected genes per ROI (default: 100).
  - `min_counts`: Minimum total counts per ROI (default: 1000).
- **Returns:** List containing processed data, DEG results, and plots.

### File: `myR/R/analysis/trajectory/trajectory_inference.R`

#### Function: `run_slingshot_from_seurat`
- **Description:** Runs Slingshot trajectory inference from a Seurat object and returns a SingleCellExperiment object with computed trajectories.
- **Parameters:**
  - `seurat_obj`: Seurat object
  - `cluster_col`: Cluster column in metadata (default: "seurat_clusters")
  - `reduction`: Dimensionality reduction to use (default: "umap")
  - `start_cluster`: Starting cluster for trajectory (optional)
  - `end_cluster`: Ending cluster(s) for trajectory (optional)
  - `assay`: Assay to use (default: "RNA")
  - `slot`: Slot to use (default: "counts")
- **Returns:** SingleCellExperiment object with Slingshot results

#### Function: `analyze_gene_dynamics`
- **Description:** Fits a Generalized Additive Model (GAM) to model gene expression changes along pseudotime, accounting for conditions.
- **Parameters:**
  - `sce`: SingleCellExperiment object from run_slingshot_from_seurat
  - `gene`: Gene name to analyze
  - `condition`: Optional condition variable for comparison
  - `lineage`: Lineage number to analyze (default: 1)
  - `plot`: Create plots (default: TRUE)
  - `output_dir`: Directory for saving plots (optional)
- **Returns:** List containing: model, summary, plot

#### Function: `process_gene_list_dynamics`
- **Description:** Processes a list of genes using analyze_gene_dynamics, with optional parallel processing.
- **Parameters:**
  - `sce`: SingleCellExperiment object
  - `genes`: Character vector of gene names
  - `condition`: Optional condition variable
  - `lineage`: Lineage number (default: 1)
  - `plot`: Create plots (default: FALSE)
  - `output_dir`: Directory for saving results
  - `n_cores`: Number of cores for parallel processing (default: 1)
- **Returns:** List of results for each gene

#### Function: `analyze_gene_dynamics_tradeSeq`
- **Description:** Fits GAMs using tradeSeq for differential expression patterns along trajectories.
- **Parameters:**
  - `sce`: SingleCellExperiment object from Slingshot
  - `gene`: Gene name to analyze
  - `n_knots`: Number of knots for spline (default: 6)
  - `conditions`: Optional condition variable
- **Returns:** List with test results

## Directory: `signatures`

### File: `myR/R/signatures/signature_discovery.R`

#### Function: `find_gene_signature`
- **Description:** Discovers gene signatures using multiple methods including Random Forest, LASSO, limma, NMF, Wilcoxon, GAM, and PCA.
- **Parameters:**
  - `data`: Seurat object, count matrix, or data.frame
  - `meta.data`: Optional metadata data.frame (required if data is not Seurat)
  - `target_var`: Column name in metadata representing the target variable
  - `target_group`: For numeric: quantile cutoff or list(low, high). For factor: levels to compare
  - `method`: One of: "tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings"
  - `n_features`: Number of top features to return (default: 50)
  - `preprocess`: Whether to normalize/scale data (default: TRUE)
  - `min_cells`: Minimum cells expressing gene (default: 10)
  - `min_pct`: Minimum percentage of cells expressing gene (default: 0.01)
  - `return_model`: Return full model object (default: FALSE)
  - `seed`: Random seed for reproducibility (default: 42)
  - `...`: Additional method-specific parameters
- **Returns:** List containing: genes, weights, scores, performance, method, model

#### Function: `print.gene_signature`
- **Description:** Print Method for gene_signature Objects.
- **Parameters:**
  - `x`: gene_signature object
  - `...`: Additional arguments (unused)
- **Returns:** None (prints to console)

#### Function: `.find_signature_rf`
- **Description:** Internal helper for Random Forest signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
  - `return_model`: Whether to return the model
- **Returns:** List of signature results

#### Function: `.find_signature_lasso`
- **Description:** Internal helper for LASSO signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
  - `return_model`: Whether to return the model
- **Returns:** List of signature results

#### Function: `.find_signature_limma`
- **Description:** Internal helper for limma signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
  - `return_model`: Whether to return the model
- **Returns:** List of signature results

#### Function: `.find_signature_wilcoxon`
- **Description:** Internal helper for Wilcoxon signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
- **Returns:** List of signature results

#### Function: `.find_signature_nmf`
- **Description:** Internal helper for NMF signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
  - `return_model`: Whether to return the model
- **Returns:** List of signature results

#### Function: `.find_signature_gam`
- **Description:** Internal helper for GAM signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
- **Returns:** List of signature results

#### Function: `.find_signature_pca`
- **Description:** Internal helper for PCA loadings signature discovery.
- **Parameters:**
  - `expr_mat`: Expression matrix
  - `target_binary`: Binary target vector
  - `n_features`: Number of features to return
  - `n_groups`: Number of groups
  - `return_model`: Whether to return the model
- **Returns:** List of signature results

### File: `myR/R/signatures/signature_scoring.R`

#### Function: `AddMultipleModuleScores`
- **Description:** Calculates module scores for multiple gene sets using Seurat's AddModuleScore. Automatically handles naming and removes the default "1" suffix.
- **Parameters:**
  - `seurat_object`: A Seurat object
  - `feature_sets`: Named or unnamed list of character vectors (gene sets)
  - `assay`: Name of the assay to use (default: DefaultAssay)
  - `slot`: Slot to pull expression data from (default: "data")
  - `nbin`: Number of bins for AddModuleScore (default: 24)
  - `ctrl`: Number of control features (default: 100)
  - `seed`: Random seed (default: 1)
  - `search`: Passed to Seurat::[.Assay (default: FALSE)
  - `...`: Additional arguments passed to Seurat::AddModuleScore
- **Returns:** Seurat object with added module scores in metadata

#### Function: `add_signature_enrichit`
- **Description:** Calculates a gene signature score using escape::enrichIt and adds it to the Seurat object's metadata. Flexibly accepts different gene ID types.
- **Parameters:**
  - `seurat_obj`: A Seurat object
  - `gene_source`: File path, data.frame, or vector containing gene IDs
  - `signature_name`: Name for the new metadata column
  - `input_keytype`: Type of input gene IDs (default: "ENSEMBL")
  - `gene_col`: Column index/name if gene_source is file or data.frame (default: 1)
  - `sheet_name`: Sheet name/index if xlsx file (default: 1)
  - `assay`: Assay to use (default: "RNA")
  - `layer`: Layer (slot) to use (default: "data")
  - `...`: Additional arguments passed to escape::enrichIt
- **Returns:** Seurat object with signature score added to metadata

#### Function: `add_progeny_scores`
- **Description:** Infers pathway activities using progeny and adds them to the Seurat object's metadata.
- **Parameters:**
  - `seurat_obj`: A Seurat object
  - `organism`: Organism: "Human" or "Mouse" (default: "Human")
  - `topn`: Number of top genes per pathway (default: 100)
  - `...`: Additional arguments passed to progeny::progeny
- **Returns:** Seurat object with pathway activity scores in metadata

#### Function: `score_signature`
- **Description:** Applies a gene signature (from find_gene_signature) to score new expression data.
- **Parameters:**
  - `expr_data`: Seurat object or expression matrix
  - `signature`: gene_signature object from find_gene_signature
  - `normalize`: Whether to z-score normalize scores (default: TRUE)
- **Returns:** Named numeric vector of signature scores

## Directory: `utilities`

### File: `myR/R/utilities/demulti_utils.R`

#### Function: `get_best_two`
- **Description:** Finds the two highest probabilities from a vector.
- **Parameters:**
  - `probs`: Numeric vector of probabilities
- **Returns:** List containing: best, best_idx, second, second_idx

#### Function: `get_barcode_mapping`
- **Description:** Assigns barcodes (cell identities) to samples based on probability matrix.
- **Parameters:**
  - `prob_matrix`: Matrix where rows are cells and columns are samples, containing assignment probabilities
  - `singlet_threshold`: Minimum probability for singlet assignment (default: 0.5)
  - `doublet_threshold`: Minimum probability for doublet assignment (default: 0.3)
  - `return_probs`: Whether to return probability values (default: FALSE)
- **Returns:** Data frame with columns: barcode, assignment, type, prob1, prob2

#### Function: `is_doublet`
- **Description:** Determines if a sample identifier represents a doublet (contains "+").
- **Parameters:**
  - `sample_name`: Character vector of sample names
- **Returns:** Logical vector indicating doublet status

#### Function: `generate_sample_values`
- **Description:** Creates a list of sample identifiers including both singlets and doublets.
- **Parameters:**
  - `n_samples`: Number of samples
  - `include_doublets`: Whether to include doublet combinations (default: TRUE)
  - `prefix`: Prefix for sample names (default: "Sample")
- **Returns:** Character vector of sample identifiers

#### Function: `generate_sample_names`
- **Description:** Creates formatted sample name lists for display/plotting.
- **Parameters:**
  - `n_samples`: Number of samples
  - `format`: Format string with %d placeholder (default: "Sample %d")
  - `include_doublets`: Whether to include doublet combinations (default: TRUE)
- **Returns:** Character vector of formatted sample names

#### Function: `parse_doublet_name`
- **Description:** Extracts the two sample identifiers from a doublet name.
- **Parameters:**
  - `doublet_name`: Character string representing a doublet (e.g., "Sample1+Sample2")
  - `separator`: Character used to separate samples (default: "+")
- **Returns:** Character vector of length 2 with the two sample names, or NULL if not a doublet

#### Function: `filter_demulti_results`
- **Description:** Filters barcode assignments based on type and quality criteria.
- **Parameters:**
  - `assignments`: Data frame from get_barcode_mapping
  - `keep_singlets`: Keep singlet assignments (default: TRUE)
  - `keep_doublets`: Keep doublet assignments (default: FALSE)
  - `keep_negative`: Keep negative/unassigned cells (default: FALSE)
  - `min_prob`: Minimum probability threshold (if prob columns exist)
- **Returns:** Filtered data frame

#### Function: `summarize_demulti_results`
- **Description:** Generates summary statistics for demultiplexing results.
- **Parameters:**
  - `assignments`: Data frame from get_barcode_mapping
- **Returns:** List containing: n_total, n_singlets, n_doublets, n_negative, doublet_rate, assignment_rate

### File: `myR/R/utilities/general_utils.R`

#### Function: `downsample_sobj`
- **Description:** This function randomly samples cells from a Seurat object to create a smaller subset. The sampling is done without replacement and maintains the original data structure.
- **Parameters:**
  - `sobj`: A Seurat object to be downsampled
  - `ratio`: Integer indicating the downsampling ratio (1:ratio). Default is 10, meaning the output will contain 1/10th of the original cells
  - `seed`: Integer specifying the random seed for reproducibility. Default is 1234
- **Returns:** A downsampled Seurat object containing a subset of cells from the input object

#### Function: `printmy`
- **Description:** Convenience function to print marker genes in a format ready for copy-paste.
- **Parameters:**
  - `markers`: Data frame with marker gene results
  - `sign`: Sign of log2FC to print ("+" for positive, "-" for negative)
  - `num`: Number of genes to print (default: 100)
  - `pseudobulk`: Whether markers are from pseudobulk analysis (uses "logFC" instead of "avg_log2FC")
- **Returns:** Prints genes to console

#### Function: `printMy`
- **Description:** Print marker genes from a named list of marker data frames.
- **Parameters:**
  - `markers_list`: Named list of marker data frames
  - `...`: Additional arguments passed to printmy()
- **Returns:** Prints genes to console with list names as headers

#### Function: `print_gene_combinations`
- **Description:** For a list of gene sets, print genes that are unique to each combination (e.g., genes in A only, genes in A & B only, etc.)
- **Parameters:**
  - `gene_list`: Named list of character vectors. Each element is a gene vector.
  - `num_print`: Integer. Maximum number of genes to print per combination (default: 100).
- **Returns:** None (prints to console)

#### Function: `%||%`
- **Description:** Returns the first non-NULL value.
- **Parameters:**
  - `a`: First value
  - `b`: Second value (returned if a is NULL)
- **Returns:** a if not NULL, otherwise b

### File: `myR/R/utilities/plot_utils.R`

#### Function: `save_plot_with_conflict_resolution`
- **Description:** Saves a plot to a file, automatically handling filename conflicts by appending a number to the filename if it already exists.
- **Parameters:**
  - `plot_object`: ggplot2 object or other plot object
  - `base_filename`: Base filename (with extension)
  - `output_dir`: Directory to save plot to
  - `width`: Plot width in inches
  - `height`: Plot height in inches
  - `dpi`: Resolution in dots per inch
  - `device`: Device to use for saving (auto-detected from extension if NULL)
- **Returns:** The full path to the saved file

#### Function: `get_default_plot_dims`
- **Description:** Returns sensible default plot dimensions based on plot type.
- **Parameters:**
  - `plot_type`: Type of plot ("heatmap", "scatter", "boxplot", "default")
  - `n_items`: Number of items (e.g., genes, clusters) - used for sizing
- **Returns:** List with width and height in inches

#### Function: `sanitize_filename`
- **Description:** Removes or replaces characters that are problematic in filenames.
- **Parameters:**
  - `filename`: Filename to sanitize
  - `replacement`: Character to replace problematic characters with (default: "_")
- **Returns:** Sanitized filename

#### Function: `create_plot_grid`
- **Description:** Helper to calculate optimal grid layout dimensions.
- **Parameters:**
  - `n_plots`: Number of plots to arrange
  - `ncol`: Desired number of columns (if NULL, calculated automatically)
  - `nrow`: Desired number of rows (if NULL, calculated automatically)
- **Returns:** List with ncol and nrow

#### Function: `add_significance_stars`
- **Description:** Converts p-values to significance stars.
- **Parameters:**
  - `p_values`: Numeric vector of p-values
  - `thresholds`: Named vector of thresholds (default: c("***"=0.001, "**"=0.01, "*"=0.05))
  - `ns_symbol`: Symbol for non-significant (default: "ns")
- **Returns:** Character vector of significance symbols

#### Function: `get_color_palette`
- **Description:** Returns a color palette for consistent styling across plots.
- **Parameters:**
  - `name`: Palette name ("default", "viridis", "Set1", "custom")
  - `n`: Number of colors needed
- **Returns:** Vector of color codes

### File: `myR/R/utilities/sample_utils.R`

#### Function: `sort_samples`
- **Description:** Sorts character sample identifiers, handling both single-number and dual-number strings (e.g., "1", "2+3", "10"). Numeric samples are sorted numerically, with doublets (containing "+") sorted after singlets.
- **Parameters:**
  - `samples`: Character vector of sample identifiers
- **Returns:** Sorted character vector of sample IDs

#### Function: `generate_sample_names`
- **Description:** Helper function to generate standardized sample names.
- **Parameters:**
  - `n_samples`: Number of samples to generate names for
  - `prefix`: Prefix for sample names (default: "Sample")
  - `start_index`: Starting index (default: 1)
- **Returns:** Character vector of sample names

#### Function: `parse_sample_ids`
- **Description:** Extracts numeric components from sample identifiers.
- **Parameters:**
  - `samples`: Character vector of sample identifiers
  - `return_doublets`: Whether to return doublet information
- **Returns:** Data frame with columns: sample, numeric1, numeric2 (if doublet), is_doublet

## Directory: `visualization`

### File: `myR/R/visualization/basic_plots.R`

#### Function: `mybar`
- **Description:** Generates a histogram for a numeric vector or a column in a data frame.
- **Parameters:**
  - `x`: Numeric vector or data frame
  - `col`: If x is data frame, column name or index
  - `bins`: Number of bins (default: 30)
  - `title`: Plot title
  - `xlab`: X-axis label
  - `color`: Fill color (default: "steelblue")
- **Returns:** ggplot2 object

#### Function: `mydensity`
- **Description:** Generates a density plot for a numeric vector or a column in a data frame.
- **Parameters:**
  - `x`: Numeric vector or data frame
  - `col`: If x is data frame, column name or index
  - `title`: Plot title
  - `xlab`: X-axis label
  - `color`: Fill color (default: "steelblue")
  - `alpha`: Transparency (default: 0.5)
- **Returns:** ggplot2 object

#### Function: `myline`
- **Description:** Creates a line plot for a numeric vector with optional cutoff line.
- **Parameters:**
  - `vec`: Numeric vector
  - `cutoff`: Optional horizontal cutoff line value
  - `show_intersection`: Show intersection point with cutoff (default: TRUE)
  - `title`: Plot title
  - `xlab`: X-axis label (default: "Index")
  - `ylab`: Y-axis label (default: "Value")
  - `color`: Line color (default: "steelblue")
- **Returns:** ggplot2 object

#### Function: `mylines`
- **Description:** Creates a line plot for multiple numeric vectors with optional dual y-axis.
- **Parameters:**
  - `vec1`: First numeric vector
  - `vec2`: Second numeric vector (optional)
  - `cutoff1`: Cutoff for vec1
  - `cutoff2`: Cutoff for vec2
  - `title`: Plot title
  - `ylab1`: Y-axis label for vec1
  - `ylab2`: Y-axis label for vec2
  - `color1`: Color for vec1 (default: "steelblue")
  - `color2`: Color for vec2 (default: "darkred")
- **Returns:** ggplot2 object

#### Function: `cdf`
- **Description:** Computes and plots the CDF for probability, logit, or ratio values.
- **Parameters:**
  - `values`: Numeric vector of values
  - `type`: Type of values: "prob", "logit", or "ratio" (default: "prob")
  - `title`: Plot title
  - `xlab`: X-axis label
  - `color`: Line color (default: "steelblue")
- **Returns:** ggplot2 object

#### Function: `cdf_multi`
- **Description:** Computes and plots CDFs across multiple or grouped datasets.
- **Parameters:**
  - `data_list`: List of numeric vectors or data frame with grouping column
  - `group_col`: If data is data.frame, column name for grouping
  - `value_col`: If data is data.frame, column name for values
  - `type`: Type of values: "prob", "logit", or "ratio" (default: "prob")
  - `title`: Plot title
  - `xlab`: X-axis label
- **Returns:** ggplot2 object

### File: `myR/R/visualization/composition_plots.R`

#### Function: `cmb`
- **Description:** Creates a proportional stacked bar graph of clusters across samples.
- **Parameters:**
  - `sobj`: Seurat object
  - `cluster_col`: Cluster identity column (default: "seurat_clusters")
  - `sample_col`: Sample grouping column
  - `colors`: Color palette (optional)
  - `title`: Plot title
  - `xlab`: X-axis label
  - `ylab`: Y-axis label (default: "Proportion")
- **Returns:** ggplot2 object

#### Function: `acmb`
- **Description:** Creates an absolute count stacked bar graph of clusters across samples.
- **Parameters:**
  - `sobj`: Seurat object
  - `cluster_col`: Cluster identity column (default: "seurat_clusters")
  - `sample_col`: Sample grouping column
  - `colors`: Color palette (optional)
  - `title`: Plot title
  - `xlab`: X-axis label
  - `ylab`: Y-axis label (default: "Cell Count")
- **Returns:** ggplot2 object

#### Function: `cml`
- **Description:** Produces a cumulative line plot of per-cluster cell proportions.
- **Parameters:**
  - `sobj`: Seurat object
  - `cluster_col`: Cluster identity column (default: "seurat_clusters")
  - `sample_col`: Sample grouping column
  - `title`: Plot title
  - `xlab`: X-axis label (default: "Sample")
  - `ylab`: Y-axis label (default: "Cumulative Proportion")
- **Returns:** ggplot2 object

#### Function: `upset_gene_lists`
- **Description:** Creates an UpSet plot for visualizing intersections of multiple gene lists.
- **Parameters:**
  - `gene_lists`: Named list of character vectors (gene lists)
  - `min_size`: Minimum intersection size to display (default: 0)
  - `nsets`: Number of sets to display (default: all)
  - `nintersects`: Number of intersections to display (default: 40)
  - `order_by`: Order by: "freq" or "degree" (default: "freq")
  - `...`: Additional arguments passed to UpSetR::upset
- **Returns:** UpSet plot

#### Function: `vln_p`
- **Description:** Generates violin plots with statistical testing using ggpubr.
- **Parameters:**
  - `sobj`: Seurat object
  - `features`: Features to plot
  - `group_by`: Grouping variable
  - `split_by`: Optional splitting variable
  - `test_method`: Statistical test: "wilcox", "t.test", "anova", "kruskal" (default: "wilcox")
  - `comparisons`: List of comparisons (default: all pairwise)
  - `ncol`: Number of columns for faceting (default: NULL)
- **Returns:** ggplot2 object

## Directory: `deprecated`

### File: `myR/R/deprecated/pseudobulk_deg_standalone.R`

#### Function: `prepare_pseudobulk_edgeR`
- **Description:** Aggregates single-cell counts to pseudo-bulk level for differential expression analysis.
- **Parameters:**
  - `object`: Seurat object
  - `sample_col`: Metadata column identifying samples/replicates
  - `group_col`: Metadata column for grouping (e.g., condition, treatment)
  - `cluster_col`: Optional cluster column for per-cluster analysis (default: NULL)
  - `target_cluster`: If cluster_col is specified, analyze only this cluster (default: NULL)
  - `assay`: Assay to use (default: "RNA")
  - `slot`: Slot to use (default: "counts")
  - `min_cells`: Minimum cells per sample to include (default: 10)
  - `min_counts`: Minimum total counts per gene (default: 10)
- **Returns:** List containing: counts, metadata, cluster

#### Function: `run_pseudobulk_deg`
- **Description:** Performs edgeR-based pseudo-bulk differential gene expression analysis.
- **Parameters:**
  - `object`: Seurat object or prepared pseudobulk list from prepare_pseudobulk_edgeR
  - `sample_col`: Metadata column identifying samples (required if object is Seurat)
  - `group_col`: Metadata column for comparison (required if object is Seurat)
  - `comparison`: Comparison to make (e.g., c("Treatment", "Control"))
  - `cluster_col`: Optional cluster column (default: NULL)
  - `target_cluster`: Specific cluster to analyze (default: NULL for overall)
  - `mode`: Analysis mode: "overall", "per_cluster", "specific_cluster" (default: "overall")
  - `assay`: Assay to use (default: "RNA")
  - `slot`: Slot to use (default: "counts")
  - `min_cells`: Minimum cells per sample (default: 10)
  - `min_counts`: Minimum total counts per gene (default: 10)
  - `fdr_threshold`: FDR threshold for significance (default: 0.05)
  - `logfc_threshold`: Log fold change threshold (default: 0)
- **Returns:** Data frame with differential expression results, or list of data frames if mode="per_cluster"

#### Function: `.run_edger_analysis`
- **Description:** Internal edgeR Analysis Function.
- **Parameters:**
  - `counts`: counts
  - `metadata`: metadata
  - `group_col`: group_col
  - `comparison`: comparison
  - `fdr_threshold`: fdr_threshold
  - `logfc_threshold`: logfc_threshold
- **Returns:** Data frame with edgeR results

## Directory: `myR/R/`

### File: `myR/R/cluster_frequency.R`

#### Function: `seurat_group_stats`
- **Description:** This function performs comprehensive statistical analysis on scRNA-seq data comparing groups across different categorical variables.
- **Parameters:**
  - `data`: A Seurat object or metadata dataframe
  - `grouping_var`: Character. Variable to group by (e.g., "patient")
  - `categorizing_var`: Character. Dependent variable (e.g., "seurat_clusters")
  - `comparative_var`: Character. Independent variable (e.g., "condition", "smoking_status")
  - `test_use`: Character vector. Statistical tests to use, options: "t", "u", or both c("t", "u")
  - `...`: Additional arguments (not used)
- **Returns:** A data frame with statistical test results

#### Function: `seurat_posthoc_analysis`
- **Description:** This function performs post-hoc analysis based on results from seurat_group_stats.
- **Parameters:**
  - `data`: A Seurat object or metadata dataframe
  - `results`: Results dataframe from seurat_group_stats function
  - `alpha`: Significance level (default: 0.05)
- **Returns:** A list of post-hoc test results

#### Function: `extract_sample_metadata`
- **Description:** Extract sample-level metadata from cell-level metadata
- **Parameters:**
  - `cell_metadata`: Cell-level metadata (Seurat object or data.frame)
  - `sample_key`: Column name for sample identification
  - `target_key`: Column name for the target variable (e.g., cluster)
  - `group_key`: Column name for grouping variable (e.g., prognosis)
  - `fun`: Aggregation function (default: proportion calculation)
- **Returns:** Data frame with sample-level aggregated data

#### Function: `perform_prior_tests`
- **Description:** Perform prior tests for statistical assumptions
- **Parameters:**
  - `data`: Data frame with sample-level data
  - `value_col`: Column name for the values to test
  - `group_col`: Column name for grouping variable
- **Returns:** List with test results

#### Function: `plot_cluster_fractions`
- **Description:** Main function to create boxplot with statistical comparisons
- **Parameters:**
  - `sobj_metadata`: Seurat object or metadata data.frame
  - `sample_key`: Column name for sample identification
  - `cluster_key`: Column name for clusters
  - `group_key`: Column name for grouping (e.g., prognosis)
  - `clusters_to_plot`: Vector of clusters to include (NULL for all)
  - `groups_to_plot`: Vector of groups to include (NULL for all)
  - `method`: Statistical test method
  - `prior_test`: Whether to perform prior tests
  - `palette`: Color palette for groups
  - `show_brackets`: Whether to show brackets for comparisons
  - `p_label_size`: Size of p-value labels
  - `...`: Additional arguments for ggplot
- **Returns:** A list containing the ggplot object and a dataframe with statistics.

#### Function: `example_usage`
- **Description:** Example usage function for `plot_cluster_fractions`
- **Parameters:** None
- **Returns:** Prints a plot and statistics table.

#### Function: `calculate_sample_signatures`
- **Description:** Calculate custom signatures at sample level
- **Parameters:**
  - `seurat_obj`: Seurat object
  - `gene_sets`: Named list of gene sets
  - `sample_key`: Column name for sample identification
  - `method`: Scoring method (mean, ssgsea, etc.)
- **Returns:** Data frame with sample-level signatures

#### Function: `adjust_confounders`
- **Description:** Adjust for confounding variables
- **Parameters:**
  - `data`: Data frame with sample-level data
  - `target_col`: Column name for target variable
  - `group_col`: Column name for grouping variable
  - `adjust_cols`: Vector of columns to adjust for
  - `method`: Adjustment method (residuals, stratification)
- **Returns:** Adjusted data frame

#### Function: `add_simple_pvalues`
- **Description:** Add p-values without brackets (similar to first image style)
- **Parameters:**
  - `plot_obj`: ggplot object
  - `stats_df`: Statistics dataframe
  - `plot_data`: Plot data with x_labels
  - `p_label_size`: Size of p-value text
- **Returns:** Modified ggplot object

### File: `myR/R/pseudotime.R`

#### Function: `run_slingshot_from_seurat`
- **Description:** This function takes a Seurat object (typically after PCA, UMAP, and clustering), runs Slingshot to infer trajectories, and returns a SingleCellExperiment object containing the Slingshot results (pseudotime and cell weights).
- **Parameters:**
  - `seurat_obj`: A Seurat object.
  - `cluster_col`: A character string specifying the column name in Seurat object's metadata that contains the cluster labels.
  - `reduced_dim_name`: A character string specifying the name of the dimensionality reduction to use (e.g., "UMAP", "PCA"). Default is "UMAP".
  - `start_cluster`: A character string or numeric value specifying the identity of the starting cluster.
  - `end_clusters`: Optional. A character vector or numeric vector specifying identities of known end clusters.
  - `counts_assay_name`: Character string, the name of the assay in the Seurat object containing counts. Default is "RNA".
  - `main_trajectory_only`: Logical. Placeholder for future implementation. Default is FALSE.
  - `...`: Additional arguments to pass to `slingshot::slingshot()`.
- **Returns:** A SingleCellExperiment object with Slingshot's results.

#### Function: `save_plot_with_conflict_resolution`
- **Description:** Saves a plot to a file, automatically handling filename conflicts by appending a number to the filename if it already exists.
- **Parameters:**
  - `plot_object`: ggplot2 object or other plot object.
  - `base_filename`: Base filename (with extension).
  - `output_dir`: Directory to save plot to.
  - `width`: Plot width in inches.
  - `height`: Plot height in inches.
  - `dpi`: Resolution in dots per inch.
- **Returns:** The full path to the saved file.

#### Function: `analyze_gene_dynamics`
- **Description:** This function fits a Generalized Additive Model (GAM) to model gene expression changes along pseudotime, potentially accounting for different conditions.
- **Parameters:**
  - `gene_id`: A character string specifying the gene ID.
  - `cds_obj`: A cell_data_set object from Monocle3.
  - `condition_col_name`: A character string naming the column in `colData(cds_obj)` that contains the condition information.
  - `pseudotime_method`: A function or character string to specify how pseudotime is extracted. Defaults to `monocle3::pseudotime`.
  - `sample_col_name`: Optional. A character string naming the column for sample/patient IDs.
  - `output_dir`: A character string specifying the directory to save plots.
  - `k_val`: An integer, the number of basis functions for GAM splines.
  - `min_cells_for_fit`: An integer, the minimum number of cells required to fit the GAM.
  - `plot_split`: Logical, if TRUE, plots for different conditions will be faceted.
  - `plot_width`: Width of the saved plot in inches.
  - `plot_height`: Height of the saved plot in inches.
  - `plot_dpi`: DPI of the saved plot.
  - `scale_DR`: Logical, if TRUE, Dynamic Range (DR) will be scaled by the mean expression.
- **Returns:** A list containing status, metrics, plot path, and the plot object.

#### Function: `process_gene_list_dynamics`
- **Description:** Helper to process a list of genes using `analyze_gene_dynamics`.
- **Parameters:**
  - `gene_list`: A character vector of gene IDs.
  - `cds_obj`: A cell_data_set object from Monocle3.
  - `condition_col_name`: Column name for conditions in `colData(cds_obj)`.
  - `output_dir`: Directory to save plots.
  - `num_cores`: Number of cores for parallel processing.
  - `...`: Other arguments to pass to `analyze_gene_dynamics`.
- **Returns:** A data.frame summarizing metrics for all successfully processed genes.

#### Function: `analyze_gene_dynamics_tradeSeq`
- **Description:** Fits GAMs using the tradeSeq framework for a specific gene and performs tests for differential progression or expression patterns along pseudotime trajectories.
- **Parameters:**
  - `gene_id`: A character string, the ID of the gene to analyze.
  - `sce_obj`: A SingleCellExperiment object, typically output from `run_slingshot_from_seurat`.
  - `condition_col`: A character string specifying the column name in `colData(sce_obj)` that contains the condition factor.
  - `lineage_names`: A character vector specifying which lineages to test.
  - `nknots`: An integer, the number of knots for GAM splines in `fitGAM`.
  - `test_to_perform`: A character string. Currently supports "patternTest".
  - `pseudotime_assay_name`: Character. Name of the reducedDim in SCE containing pseudotime. Default "slingshot".
  - `weights_col_prefix`: Character. Prefix for colData columns containing cell weights. Default "slingWeight".
  - `output_dir`: A character string specifying the directory to save plots.
  - `plot_split`: Logical, if TRUE, plots for different conditions will be faceted. Default FALSE.
  - `scale_DR`: Logical, if TRUE, Dynamic Range (DR) will be scaled by the mean expression.
  - `fitGAM_args`: A list of additional arguments to pass to `tradeSeq::fitGAM`.
  - `test_args`: A list of additional arguments to pass to the chosen test function.
- **Returns:** A list containing status, metrics, plot path, and the plot object.

### File: `myR/R/GeoMx.R`

#### Function: `prepare_geomx_data`
- **Description:** Prepares raw GeoMx count data for analysis by organizing into a standard format.
- **Parameters:**
  - `raw_data`: Raw data matrix.
  - `metadata_SP`: Metadata data frame.
- **Returns:** List containing: expression, metadata, and gene_info.

#### Function: `q3_normalize`
- **Description:** Performs Q3 (third quartile) normalization on GeoMx count data.
- **Parameters:**
  - `expr_matrix`: Count matrix.
  - `scaling_factor`: Scaling factor after normalization (default: 1000).
- **Returns:** List containing normalized and log-normalized data.

#### Function: `perform_qc`
- **Description:** Performs basic QC checks and filtering on GeoMx data.
- **Parameters:**
  - `expr_matrix`: Count matrix.
  - `metadata`: Metadata data frame.
- **Returns:** List with QC metrics and plots.

#### Function: `find_deg_geomx`
- **Description:** Performs differential expression analysis on GeoMx data using limma or Wilcoxon test.
- **Parameters:**
  - `norm_expr`: Normalized expression matrix.
  - `metadata`: Metadata data frame.
  - `group_var`: Metadata column for grouping.
  - `group1`: First group for comparison.
  - `group2`: Second group for comparison (optional).
  - `method`: "limma" or "wilcox" (default: "limma").
  - `logFC_threshold`: Log fold change threshold (default: 0.5).
  - `pval_threshold`: P-value threshold (default: 0.05).
- **Returns:** Data frame with differential expression results.

#### Function: `plot_deg_volcano`
- **Description:** Creates a volcano plot from differential expression results.
- **Parameters:**
  - `deg_results`: Data frame from `find_deg_geomx`.
  - `title`: Plot title.
- **Returns:** ggplot object.

#### Function: `plot_deg_heatmap`
- **Description:** Creates a heatmap of top differentially expressed genes.
- **Parameters:**
  - `norm_expr`: Normalized expression matrix.
  - `deg_results`: Data frame from `find_deg_geomx`.
  - `metadata`: Metadata data frame.
  - `group_var`: Metadata column for grouping.
  - `top_n`: Number of top genes to plot (default: 50).
  - `scale_rows`: Z-score scale rows (default: TRUE).
- **Returns:** Heatmap plot.

#### Function: `run_geomx_analysis`
- **Description:** Wrapper function to run the complete GeoMx analysis workflow.
- **Parameters:**
  - `raw_data`: Raw data matrix.
  - `metadata_SP`: Metadata data frame.
  - `group_var`: Column for differential expression.
  - `group1`: First group for comparison.
  - `group2`: Second group for comparison (optional).
  - `qc_filter`: Perform QC filtering (default: TRUE).
  - `min_genes`: Minimum number of detected genes per ROI (default: 100).
  - `min_counts`: Minimum total counts per ROI (default: 1000).
- **Returns:** List containing processed data, DEG results, and plots.

### File: `myR/R/markers.R`

#### Function: `marker_trim`
- **Description:** Filter and process marker genes from Seurat's FindMarkers or FindAllMarkers results.
- **Parameters:**
  - `markers`: A data frame from Seurat's FindMarkers or FindAllMarkers results.
  - `sign`: Character string indicating the direction of log fold change to keep. Options: NULL (keep all), "+" (positive only), "-" (negative only).
  - `p_cutoff`: Numeric value for adjusted p-value cutoff. Default is NULL (no filtering).
  - `filter`: Additional filtering criteria (not used in this function).
- **Returns:** A processed data frame with filtered markers and added pct.diff column.

#### Function: `marker_filter`
- **Description:** Filter out unwanted genes from marker results.
- **Parameters:**
  - `markers`: A data frame from Seurat's FindMarkers or FindAllMarkers results.
  - `filter`: Character vector specifying which gene types to filter out. Options: "rb", "mt", "hb", "AC", "ENSG", "LINC".
- **Returns:** A filtered data frame with unwanted genes removed.

#### Function: `lrf`
- **Description:** filterout (Intercept) terms.
- **Parameters:**
  - `lmm_result`: The result from an LMM analysis.
- **Returns:** The LMM result with "(Intercept)" terms removed from the summary.

#### Function: `all_markers_to_list`
- **Description:** Convert FindAllMarkers results to a list organized by cluster.
- **Parameters:**
  - `markers`: A data frame from Seurat's FindAllMarkers results.
- **Returns:** A list where each element is a data frame containing markers for one cluster.

#### Function: `marker_print_all`
- **Description:** Print top marker genes for each cluster.
- **Parameters:**
  - `markers`: Either a list of marker data frames or a FindAllMarkers data frame.
  - `n`: Integer specifying the number of top markers to print per cluster.
  - `cluster_to_print`: Character vector specifying which clusters to print. If NULL, prints all clusters.
- **Returns:** NULL (prints results to console).

#### Function: `marker_print`
- **Description:** Print top n marker genes.
- **Parameters:**
  - `marker`: A marker data frame.
  - `n`: Number of genes to print (default: 100).
  - `sign`: "+" for positive logFC, "-" for negative.
- **Returns:** NULL (prints results to console).

#### Function: `synthesize_markers`
- **Description:** Combines p-values and log fold changes from multiple FindMarkers results using Fisher's method for p-values and weighted mean for log fold changes.
- **Parameters:**
  - `marker_list`: Named list of marker data frames
  - `p_sig`: The p-value column to use for significance check (default: "p_val_adj").
- **Returns:** Data frame with synthesized results

#### Function: `synthesize_ranks`
- **Description:** Combines ranks from multiple FindMarkers results using geometric mean of ranks.
- **Parameters:**
  - `marker_list`: Named list of marker data frames
  - `rank_by`: Column to rank by (default: "p_val").
- **Returns:** Data frame with synthesized ranks.

### File: `myR/R/plots.R`

#### Function: `mybar`
- **Description:** Generates a histogram for a numeric vector or a column in a data frame.
- **Parameters:**
  - `data`: Numeric vector or data frame.
  - `column`: If data is data frame, column name or index.
  - `bins`: Number of bins (default: NULL).
  - `x_unit`: Unit for x-axis label.
  - `y_unit`: Unit for y-axis label.
  - `xlab`: X-axis label.
  - `main`: Plot title.
- **Returns:** Prints a ggplot object.

#### Function: `myline`
- **Description:** Creates a line plot for a single vector with cutoff intersection.
- **Parameters:**
  - `vector`: Numeric vector.
  - `main`: Plot title.
  - `ylab`: Y-axis label.
  - `cutoff`: Optional horizontal cutoff line value.
  - `cutoff_col`: Color for the cutoff line.
  - `show_intersection`: Show intersection point with cutoff (default: TRUE).
  - `...`: Additional arguments to `plot()`.
- **Returns:** A base R plot.

#### Function: `mylines`
- **Description:** Creates a line plot for multiple numeric vectors with dual y-axis support.
- **Parameters:**
  - `...`: Numeric vectors for the left y-axis.
  - `vectors_right`: A single vector or a list of vectors for the right y-axis.
  - `main`: Plot title.
  - `ylab`: Left y-axis label.
  - `ylab_right`: Right y-axis label.
  - `cutoff`: Cutoff for left y-axis vectors.
  - `cutoff_right`: Cutoff for right y-axis vectors.
  - `cutoff_col`: Color for the left cutoff line.
  - `cutoff_col_right`: Color for the right cutoff line.
  - `legend_pos`: Position of the legend.
  - `show_intersection`: Show intersection points (default: TRUE).
- **Returns:** A base R plot with multiple lines.

#### Function: `find_and_mark_intersection`
- **Description:** Helper function to find and mark the intersection point on a line plot.
- **Parameters:**
  - `x_vals`: X values of the line.
  - `y_vals`: Y values of the line.
  - `cutoff`: The cutoff value.
  - `col`: Color for the intersection point and label.
  - `label_prefix`: Prefix for the intersection label.
- **Returns:** Adds points and text to an existing plot.

#### Function: `mydensity`
- **Description:** Generates a density plot for a numeric vector or a column in a data frame.
- **Parameters:**
  - `data`: Numeric vector or data frame.
  - `column`: If data is data frame, column name or index.
  - `adjust`: Bandwidth adjustment.
  - `x_unit`: Unit for x-axis label.
  - `y_unit`: Unit for y-axis label.
  - `xlab`: X-axis label.
  - `main`: Plot title.
- **Returns:** Prints a ggplot object.

#### Function: `mybox`
- **Description:** Generate boxplots of average feature expression from a Seurat object.
- **Parameters:**
  - `sobj`: A Seurat object.
  - `features`: A character vector of features to plot.
  - `sample_col`: The name of the metadata column identifying samples/patients.
  - `group.by`: The name of the metadata column for primary grouping on the x-axis.
  - `split.by`: The name of the metadata column for secondary grouping.
  - `idents`: A character vector of values from the `group.by` column to include.
  - `assay`: The assay to use for gene features.
  - `layer`: The layer to use for gene features.
  - `ncol`: Number of columns for arranging plots.
  - `pt.size`: Size of points for individual sample averages.
  - `violin`: Whether to overlay violin plots.
- **Returns:** A ggplot object.

#### Function: `mybox_df`
- **Description:** Generate boxplots of average feature values from a data frame.
- **Parameters:**
  - `df`: A data frame.
  - `features`: A character vector of feature column names to plot.
  - `sample_col`: The name of the column identifying samples/patients.
  - `group.by`: The name of the column for primary grouping on the x-axis.
  - `split.by`: The name of the column for secondary grouping.
  - `idents`: A character vector of values from the `group.by` column to include.
  - `ncol`: Number of columns for arranging plots.
  - `pt.size`: Size of points for individual sample averages.
  - `violin`: Whether to overlay violin plots.
- **Returns:** A ggplot object.

#### Function: `upset_gene_lists`
- **Description:** Creates an UpSet plot for visualizing intersections of multiple gene lists.
- **Parameters:**
  - `gene_lists`: Named list of character vectors (gene lists).
  - `sets`: Optional character vector to specify the order of sets.
  - `min_size`: Minimum intersection size to display.
  - `width_ratio`: Ratio of main plot and bar plot widths.
  - `keep_sets`: Character vector of specific sets to include in the plot.
- **Returns:** A ggplot object.

#### Function: `vln_p`
- **Description:** Creates violin plots with statistical comparisons.
- **Parameters:**
  - `sobj`: Seurat object.
  - `feature`: Features to plot.
  - `group.by`: Grouping variable.
  - `split.by`: Splitting variable.
  - `pt.size`: Point size.
  - `ncol`: Number of columns for faceting.
  - `...`: Additional arguments passed to `VlnPlot`.
- **Returns:** A ggplot object.

#### Function: `cmb`
- **Description:** Creates a proportional stacked bar graph of clusters.
- **Parameters:**
  - `sobj`: A Seurat object.
  - `identity`: The identity to use for clustering (default: "seurat_clusters").
  - `group.by`: The metadata column to group by (default: "sample").
  - `idents`: Which identities to include.
  - `df`: If TRUE, returns the data frame used for plotting.
  - `vlines`: X-axis positions for vertical lines.
  - `vline_color`: Color for vertical lines.
- **Returns:** A ggplot object or a data frame.

#### Function: `acmb`
- **Description:** Creates an absolute count stacked bar graph of clusters.
- **Parameters:**
  - `sobj`: A Seurat object.
  - `identity`: The identity to use for clustering (default: "seurat_clusters").
  - `group.by`: The metadata column to group by (default: "sample").
  - `idents`: Which identities to include.
  - `df`: If TRUE, returns the data frame used for plotting.
  - `vlines`: X-axis positions for vertical lines.
  - `vline_color`: Color for vertical lines.
- **Returns:** A ggplot object or a data frame.

#### Function: `myhm_genesets2`
- **Description:** Creates a heatmap of gene set expression per cluster.
- **Parameters:**
  - `sobj`: A Seurat object.
  - `group`: The identity to use for clustering.
  - `value`: How to aggregate expression values ("average" or "sum").
  - `assay`: Which assay to use.
  - `gene_sets`: A list of character vectors (gene sets).
  - `title`: Plot title.
  - `x_label`: X-axis label.
  - `y_label`: Y-axis label.
- **Returns:** A data frame with normalized expression values.

#### Function: `myhm_genes2`
- **Description:** Creates a heatmap of individual gene expression per cluster.
- **Parameters:**
  - `sobj`: A Seurat object.
  - `group`: The identity to use for clustering.
  - `value`: How to aggregate expression values ("average" or "sum").
  - `assay`: Which assay to use.
  - `genes`: Character vector of gene names to plot.
  - `title`: Plot title.
  - `x_label`: X-axis label.
  - `y_label`: Y-axis label.
- **Returns:** A data frame with normalized expression values.

#### Function: `cml`
- **Description:** Produces a cumulative line plot of per-cluster cell proportions.
- **Parameters:**
  - `sobj`: A Seurat object.
  - `cluster_col`: Metadata column with cluster IDs.
  - `group.by`: Metadata column for grouping.
  - `sort.by`: Cluster ordering strategy ("name" or "frequency").
  - `df`: If TRUE, returns the data frame used for plotting.
  - `color_palette`: Optional vector of colors.
  - `n_patterns`, `n_shapes`: Number of linetypes/shapes to cycle.
- **Returns:** A ggplot object or a data frame.

#### Function: `cdf`
- **Description:** Computes and plots the Cumulative Distribution Function (CDF).
- **Parameters:**
  - `data`: A data frame.
  - `probability_col`: Column with probability values.
  - `ratio_col`: Column with probability ratio values.
  - `plot_type`: "probability", "logit", or "ratio".
  - `output_file`: Path to save the plot.
- **Returns:** A list with the ggplot object and a summary data frame.

#### Function: `cdf_multi`
- **Description:** Computes and plots CDFs across multiple or grouped datasets.
- **Parameters:**
  - `data_list`: A data frame or a list of data frames.
  - `probability_col`: Column with probability values.
  - `ratio_col`: Column with probability ratio values.
  - `plot_type`: "probability", "logit", or "ratio".
  - `group_by_col`: Column to group by if `data_list` is a single data frame.
  - `sample_col`: Column to name datasets.
  - `output_file`: Path to save the plot.
- **Returns:** A list with the ggplot object and a summary data frame.

#### Function: `scatter_smooth_colored`
- **Description:** Averaged Expression vs Numeric Covariate (with optional colour).
- **Parameters:**
  - `object`: A Seurat object or a data.frame.
  - `feature`: Gene symbol or expression column name.
  - `group.by`: Column used to aggregate cells into samples.
  - `x_var`: Column providing the numeric predictor per cell.
  - `transpose`: If TRUE, swap X and Y axes.
  - `color_by`: Column name to color points.
  - `palette`: Optional palette vector/name.
  - `transparency`: Map `color_by` (numeric only) to point alpha.
  - `transparency_desc`: If TRUE, higher values become more transparent.
  - `fitted_line`: "linear", "loess", "lasso", or NULL.
- **Returns:** A ggplot scatter-smooth plot.

#### Function: `mybox_geomx`
- **Description:** Generate boxplots of average feature values from GeoMx data.
- **Parameters:**
  - `data_norm`: Normalized data matrix.
  - `metadata`: Metadata data frame.
  - `features`: A character vector of feature column names to plot.
  - `group.by`: The name of the column for primary grouping.
  - `split.by`: The name of the column for secondary grouping.
  - `test_method`: "t.test" or "wilcox.test".
  - `p_adjust_method`: P-value adjustment method.
  - `hide_ns`: Hide non-significant p-values.
  - `show_points`: Show individual data points.
  - `comparisons`: "all" or a list of specific comparisons.
- **Returns:** Prints ggplot objects and invisibly returns a list of statistical results.

#### Function: `myhm_genesets4`
- **Description:** Creates a heatmap of gene set expression per cluster (version 4).
- **Parameters:**
  - `sobj`: A Seurat object.
  - `group`: The identity to use for clustering.
  - `value`: How to aggregate expression values ("average" or "sum").
  - `assay`: Which assay to use.
  - `gene_sets`: A list of character vectors (gene sets).
  - `title`: Plot title.
  - `x_label`: X-axis label.
  - `y_label`: Y-axis label.
- **Returns:** A data frame with normalized expression values.

#### Function: `myhm_genes4`
- **Description:** Creates a heatmap of individual gene expression per cluster (version 4).
- **Parameters:**
  - `sobj`: A Seurat object.
  - `group`: The identity to use for clustering.
  - `value`: How to aggregate expression values ("average" or "sum").
  - `assay`: Which assay to use.
  - `genes`: Character vector of gene names to plot.
  - `title`: Plot title.
  - `x_label`: X-axis label.
  - `y_label`: Y-axis label.
- **Returns:** A data frame with normalized expression values.

#### Function: `%||%`
- **Description:** NULL coalescing operator.
- **Parameters:**
  - `a`: First value.
  - `b`: Second value.
- **Returns:** `a` if not NULL, otherwise `b`.

### File: `myR/R/pathway_analysis.R`

#### Function: `convert_gene_ids`
- **Description:** This function converts gene identifiers between different formats using the org.Hs.eg.db annotation package.
- **Parameters:**
  - `genes`: Character vector of gene identifiers to convert.
  - `from`: Character string specifying the input ID type (default: "SYMBOL").
  - `to`: Character string specifying the output ID type (default: "ENTREZID").
- **Returns:** Data frame with mapping between input and output gene IDs.

#### Function: `prepare_gene_lists`
- **Description:** This function processes differential expression results to create gene lists suitable for various pathway enrichment analyses.
- **Parameters:**
  - `deg_df`: Data frame containing differential expression results.
  - `fc_threshold`: Numeric threshold for log2 fold change (default: 0.25).
  - `p_use`: Character string specifying which p-value to use (default: "p_val_adj").
  - `pval_threshold`: Numeric threshold for p-value significance (default: 0.05).
- **Returns:** List containing ranked genes, up/down-regulated genes, and background.

#### Function: `get_pathway_sets`
- **Description:** Returns available pathway database identifiers for enrichment analysis.
- **Parameters:**
  - `pathway_set`: Character vector of specific pathway sets to return, or NULL to return all.
  - `species`: Character string specifying species (default: "Homo sapiens").
- **Returns:** Character vector of pathway set identifiers.

#### Function: `run_go_analysis`
- **Description:** Performs Gene Ontology enrichment analysis using clusterProfiler.
- **Parameters:**
  - `genes`: Character vector of gene symbols for analysis.
  - `gene_type`: Character string describing gene set type (for logging).
  - `ont`: Character string specifying GO ontology: "BP", "MF", or "CC".
  - `background`: Character vector of background genes (optional).
  - `pval_cutoff`: Numeric p-value cutoff for significance (default: 0.05).
- **Returns:** enrichResult object from clusterProfiler, or NULL if no results.

#### Function: `run_kegg_analysis`
- **Description:** Performs KEGG pathway enrichment analysis using clusterProfiler.
- **Parameters:**
  - `genes`: Character vector of gene symbols for analysis.
  - `background`: Character vector of background genes (optional).
  - `pval_cutoff`: Numeric p-value cutoff for significance (default: 0.05).
- **Returns:** enrichResult object from clusterProfiler, or NULL if no results.

#### Function: `run_gsea_analysis`
- **Description:** Performs GSEA using fgsea with MSigDB gene sets.
- **Parameters:**
  - `ranked_genes`: Named numeric vector of genes ranked by fold change.
  - `pathway_set`: Character string specifying pathway database (default: "H").
  - `min_size`: Minimum gene set size (default: 15).
  - `max_size`: Maximum gene set size (default: 500).
- **Returns:** Data frame with GSEA results.

#### Function: `format_results`
- **Description:** Standardizes the output format from different pathway analysis methods.
- **Parameters:**
  - `result`: Results object from pathway analysis.
  - `analysis_type`: Character string specifying analysis type: "GO", "KEGG", or "GSEA".
- **Returns:** Data frame with standardized columns.

#### Function: `myGO`
- **Description:** Main function to perform comprehensive pathway enrichment analysis on differential expression results.
- **Parameters:**
  - `DEG`: Data frame with differential expression results.
  - `seurat_obj`: Seurat object (currently not used).
  - `pathway_set`: Character vector specifying which pathway databases to use.
  - `pathway`: Character string to filter results for specific pathway names (optional).
  - `analysis_type`: "GO", "KEGG", "GSEA", or "ALL".
  - `go_ontology`: GO ontologies: "BP", "MF", "CC".
  - `fc_threshold`: Log2 fold change threshold.
  - `p_use`: Which p-value column to use.
  - `pval_threshold`: P-value significance threshold.
  - `gsea_min_size`: Minimum gene set size for GSEA.
  - `gsea_max_size`: Maximum gene set size for GSEA.
  - `return_plots`: Logical indicating whether to return plots (not implemented).
- **Returns:** Named list containing analysis results for each method.

### File: `myR/R/signature.R`

#### Function: `AddMultipleModuleScores`
- **Description:** Calculates module scores for multiple gene sets using Seurat's AddModuleScore.
- **Parameters:**
  - `seurat_object`: A Seurat object.
  - `feature_sets`: A named or unnamed list of character vectors (gene sets).
  - `assay`: Name of the assay to use.
  - `slot`: Slot to pull expression data from.
  - `nbin`: Number of bins for AddModuleScore.
  - `ctrl`: Number of control features.
  - `seed`: Random seed.
  - `search`: Passed to Seurat::`[.Assay`.
  - `...`: Additional arguments passed to Seurat::AddModuleScore.
- **Returns:** A Seurat object with added module scores in the metadata.

#### Function: `PlotModuleScoreHeatmap`
- **Description:** Visualizes module scores as a heatmap with z-score normalization and statistical testing.
- **Parameters:**
  - `sobj`: Seurat object.
  - `gene_sets`: Named list of gene sets.
  - `group`: Grouping variable.
  - `assay`: Assay to use.
  - `test_method`: "wilcox" or "t".
  - `p_adjust`: P-value adjustment method.
  - `show_pval`: Whether to show p-values on the heatmap.
  - `scale_method`: "feature" or "group".
  - `color_limits`: Manual color scale limits.
  - `title`: Plot title.
  - `x_label`: X-axis label.
  - `y_label`: Y-axis label.
  - `...`: Additional arguments passed to AddModuleScore.
- **Returns:** A list containing the heatmap plot and statistical results.

#### Function: `CompareModuleScoringMethods`
- **Description:** Compares the simple averaging method with AddModuleScore.
- **Parameters:**
  - `sobj`: Seurat object.
  - `gene_sets`: Named list of gene sets.
  - `group`: Grouping variable.
  - `assay`: Assay to use.
- **Returns:** A list containing results from both methods and a comparison plot.

#### Function: `add_signature_enrichit`
- **Description:** Calculates a gene signature score using escape::enrichIt and adds it to metadata.
- **Parameters:**
  - `seurat_obj`: A Seurat object.
  - `gene_source`: Path to a file or an R object with gene IDs.
  - `signature_name`: Name for the new metadata column.
  - `input_keytype`: Type of input gene IDs (e.g., "ENSEMBL", "SYMBOL").
  - `gene_col`: Column index/name for genes.
  - `sheet_name`: Sheet name/index for xlsx files.
  - `assay`: Assay to use.
  - `layer`: Layer (slot) to use.
  - `...`: Additional arguments for `escape::enrichIt`.
- **Returns:** A Seurat object with the new signature score.

#### Function: `add_progeny_scores`
- **Description:** Infers pathway activities using progeny and adds them to metadata.
- **Parameters:**
  - `seurat_obj`: A Seurat object.
  - `organism`: "Human" or "Mouse".
  - `topn`: Number of top genes per pathway.
  - `...`: Additional arguments for `progeny::progeny`.
- **Returns:** A Seurat object with pathway activity scores.

#### Function: `linear_seurat`
- **Description:** Performs linear regression for gene expression.
- **Parameters:**
  - `sobj`: Seurat object.
  - `layer`: "counts", "data", or "scale.data".
  - `features`: Features to test.
  - `regressor`: Regressor variable in metadata.
  - `regressor.type`: "continuous", "categorical", or "ordinal".
  - `reference.level`: Reference level for categorical regressor.
  - `ordinal.method`: "linear", "polynomial", or "spline".
  - `link.function`: "linear", "poisson", or "negative.binomial".
  - `effect`: "fixed" or "random".
  - `covariates`: Covariate column names.
  - `min.cells`: Minimum cells expressing a gene.
  - `return.full`: Return full results including Seurat object.
  - `...`: Additional arguments for model fitting.
- **Returns:** Data frame with regression results.

#### Function: `plot_top_genes`
- **Description:** Helper function to visualize top results from `linear_seurat`.
- **Parameters:**
  - `results`: Results from `linear_seurat`.
  - `sobj`: Seurat object.
  - `layer`: Layer to get expression data from.
  - `top_n`: Number of top genes to plot.
- **Returns:** A grid of ggplot objects.

#### Function: `example_usage`
- **Description:** Example usage for `linear_seurat`.
- **Parameters:** None.
- **Returns:** None.

#### Function: `find_gene_signature`
- **Description:** Discovers gene signatures that best separate a target variable.
- **Parameters:**
  - `data`: Seurat object, count matrix, or data.frame.
  - `meta.data`: Optional metadata.
  - `target_var`: Target variable column name.
  - `target_group`: Groups to compare.
  - `method`: "tree_based", "lasso", "limma", etc.
  - `n_features`: Number of top features to return.
  - `preprocess`: Whether to normalize/scale data.
  - `min_cells`: Minimum cells expressing a gene.
  - `min_pct`: Minimum percentage of cells expressing a gene.
  - `return_model`: Whether to return the full model object.
  - `seed`: Random seed.
  - `...`: Additional method-specific parameters.
- **Returns:** A list containing the gene signature and performance metrics.

#### Function: `score_signature`
- **Description:** Applies a gene signature to score new expression data.
- **Parameters:**
  - `expr_data`: Seurat object or expression matrix.
  - `signature`: A `gene_signature` object from `find_gene_signature`.
  - `normalize`: Whether to z-score normalize scores.
- **Returns:** Named numeric vector of signature scores.

#### Function: `print.gene_signature`
- **Description:** Print method for `gene_signature` objects.
- **Parameters:**
  - `x`: A `gene_signature` object.
  - `...`: Additional arguments.
- **Returns:** Prints a summary to the console.

### File: `myR/R/test_claude.R`

#### Function: `create_analysis_config`
- **Description:** Centralizes metadata column names for consistent reference throughout analysis.
- **Parameters:**
  - `patient`: String. Column name for patient/subject identifier.
  - `drug`: String. Column name for drug/treatment type.
  - `timepoint`: String. Column name for timepoint (e.g., "pre"/"post").
  - `ck`: String. Column name for stratification variable (e.g., CK status).
  - `response`: String. Column name for treatment response classification.
  - `aoi`: String. Column name for AOI (Area of Interest) unique identifier.
- **Returns:** A list containing standardized column name mappings.

#### Function: `validate_config`
- **Description:** Checks that all required columns exist in metadata and have valid values.
- **Parameters:**
  - `metadata`: Data frame. Seurat metadata or standalone metadata.
  - `config`: List. Configuration object from create_analysis_config().
  - `required_cols`: Character vector. Which config columns are required.
- **Returns:** Invisibly returns TRUE if valid, otherwise throws informative error.

#### Function: `prepare_geomx_data`
- **Description:** Reads count matrix and metadata from files (Excel or CSV), creates a Seurat object, and adds derived metadata columns for convenience in downstream analysis.
- **Parameters:**
  - `count_file`: String. Path to count data file (.xlsx or .csv).
  - `metadata_file`: String. Path to metadata file (.xlsx or .csv).
  - `count_matrix`: Matrix. Raw count matrix (genes x samples).
  - `metadata`: Data.frame. Sample metadata.
  - `config`: List. Configuration from create_analysis_config().
  - `normalize_method`: String. One of "none", "log", or "quantile".
  - `min_cells`: Integer. Minimum cells expressing a gene to keep it.
  - `min_features`: Integer. Minimum features per cell to keep it.
- **Returns:** A Seurat object with added metadata columns.

#### Function: `diagnose_sample_parity`
- **Description:** For paired/longitudinal analyses, checks whether each patient (or group) has complete and balanced representation of required conditions.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `config`: List. Configuration from create_analysis_config().
  - `grouping_vars`: Character vector. Metadata columns defining groups.
  - `check_var`: String. The metadata column to check for parity.
  - `required_values`: Character vector. Values that MUST be present in check_var for each group.
  - `allow_extra`: Logical. If FALSE, groups with additional values will fail validation.
- **Returns:** A list with summary, passed samples, failed groups, and a message.

#### Function: `screen_genes`
- **Description:** Performs rapid statistical screening to identify candidate genes for computationally intensive LMM analysis.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `config`: List. Configuration from create_analysis_config().
  - `grouping_var`: String. Which variable to compare.
  - `subset_var`: String. Optional. Metadata column to subset by.
  - `subset_value`: Value to filter subset_var by.
  - `min_pct`: Numeric. Minimum fraction of samples expressing gene in either group.
  - `logfc_threshold`: Numeric. Minimum log fold-change threshold.
  - `top_n`: Integer. Number of top genes to return.
  - `test_method`: String. Statistical test method (default "wilcox").
- **Returns:** A list containing all markers, top genes, and the subsetted Seurat object.

#### Function: `fit_single_gene_lmm`
- **Description:** Internal helper function that fits an lmer model for one gene.
- **Parameters:**
  - `gene`: String. Gene name.
  - `expr_vector`: Numeric vector. Expression values for the gene.
  - `metadata`: Data frame. Sample metadata.
  - `formula_str`: String. Complete model formula.
  - `factor_cols`: Character vector. Columns to convert to factors.
  - `reference_levels`: Named list. Factor reference levels.
- **Returns:** List with model results or error information.

#### Function: `run_lmm_analysis`
- **Description:** Main function for LMM analysis. Fits mixed-effects models to account for patient-level random effects.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `genes`: Character vector. Gene names to analyze.
  - `config`: List. Configuration from create_analysis_config().
  - `fixed_effects`: Character vector. Fixed effect terms.
  - `interactions`: Character vector. Interaction terms.
  - `random_effects`: String. Random effects formula.
  - `reference_levels`: Named list. Reference levels for factors.
  - `n_cores`: Integer. Number of CPU cores for parallel processing.
  - `verbose`: Logical. Print progress messages.
- **Returns:** A list with LMM results.

#### Function: `extract_significant_genes`
- **Description:** Filters LMM summary table to genes with significant effects for specific terms.
- **Parameters:**
  - `lmm_summary`: Data frame. The $summary component from run_lmm_analysis().
  - `term_pattern`: String. Regex pattern to match model terms.
  - `p_threshold`: Numeric. Adjusted p-value threshold.
  - `effect_threshold`: Numeric. Minimum absolute effect size.
  - `top_n`: Integer. Maximum number of genes to return.
  - `rank_by`: String. How to rank genes: "effect" or "pvalue".
- **Returns:** Data frame of significant genes with their statistics.

#### Function: `find_drug_response_genes`
- **Description:** Identifies genes where treatment response differs by drug type.
- **Parameters:**
  - `lmm_summary`: Data frame from run_lmm_analysis()$summary.
  - `config`: List. Configuration object.
  - `focus_drug`: String. Optional. Specific drug to focus on.
  - `...`: Additional arguments passed to extract_significant_genes().
- **Returns:** Data frame of genes with drug-specific response patterns.

#### Function: `find_temporal_response_genes`
- **Description:** Identifies genes where expression changes over time depend on treatment response.
- **Parameters:**
  - `lmm_summary`: Data frame from run_lmm_analysis()$summary.
  - `config`: List. Configuration object.
  - `...`: Additional arguments passed to extract_significant_genes().
- **Returns:** Data frame of genes with time-by-response interactions.

#### Function: `plot_volcano`
- **Description:** Visualizes effect sizes and p-values for a specific model term.
- **Parameters:**
  - `lmm_summary`: Data frame from run_lmm_analysis()$summary.
  - `term_pattern`: String. Regex to filter model terms.
  - `effect_threshold`: Numeric. Threshold for "large effect".
  - `p_threshold`: Numeric. P-value significance threshold.
  - `label_top`: Integer. Number of top genes to label.
  - `title`: String. Plot title.
- **Returns:** ggplot object.

#### Function: `plot_lmm_results`
- **Description:** Visualizes individual patient trajectories and estimated marginal means from a fitted LMM.
- **Parameters:**
  - `lmm_result`: List. Result for a single gene from run_lmm_analysis()$results.
  - `gene`: String. Gene name (for titles).
  - `seurat_obj`: Seurat object (for raw data).
  - `config`: List. Configuration object.
  - `facet_by`: String. Optional metadata column to facet by.
- **Returns:** List with three ggplot objects: emmeans, trajectories, and contrasts.

#### Function: `plot_gene_boxplot`
- **Description:** Simple visualization for a single gene's expression across timepoints.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `gene`: String. Gene name.
  - `config`: List. Configuration object.
  - `split_by`: String. Optional metadata column to facet by.
  - `add_stats`: Logical. Add statistical comparison.
- **Returns:** ggplot object.

#### Function: `create_delta_matrix`
- **Description:** Alternative analysis approach that calculates per-patient change scores (post - pre) for each gene.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `config`: List. Configuration object.
  - `assay`: String. Assay name.
  - `slot`: String. Data slot.
  - `subset_var`: String. Optional variable to subset by.
  - `subset_value`: Value to filter subset_var by.
  - `aggregate_fun`: Function to aggregate multiple AOIs per patient-timepoint.
- **Returns:** List with delta_matrix and patient_metadata.

#### Function: `export_lmm_results`
- **Description:** Saves LMM analysis results to a multi-sheet Excel file.
- **Parameters:**
  - `lmm_results`: List from run_lmm_analysis().
  - `output_file`: String. Output Excel file path.
  - `include_model_objects`: Logical. Save model objects to RDS file.
  - `top_genes_n`: Integer. Number of top genes to highlight.
- **Returns:** Invisibly returns the output file path.

#### Function: `generate_lmm_report`
- **Description:** Creates an HTML summary report with key visualizations and tables.
- **Parameters:**
  - `lmm_results`: List from run_lmm_analysis().
  - `seurat_obj`: Seurat object (needed for plotting).
  - `config`: List. Configuration object.
  - `output_file`: String. Output HTML file path.
  - `top_genes`: Integer. Number of top genes to visualize.
- **Returns:** Invisibly returns the output file path.

#### Function: `print.geomx_config`
- **Description:** Print method for geomx_config objects.
- **Parameters:**
  - `config`: List from create_analysis_config().
- **Returns:** Invisibly returns the config object.

#### Function: `summarize_geomx_object`
- **Description:** Quick summary of Seurat object for GeoMx data.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `config`: List. Configuration object.
- **Returns:** Invisibly returns the Seurat object.

#### Function: `extract_gene_coefficients`
- **Description:** Convenience function to get detailed coefficient info for genes of interest.
- **Parameters:**
  - `lmm_results`: List from run_lmm_analysis().
  - `genes`: Character vector of gene names.
- **Returns:** Data frame of coefficients for specified genes.

#### Function: `compare_screening_vs_lmm`
- **Description:** Diagnostic function to see how well initial screening predicts LMM significance.
- **Parameters:**
  - `screening_results`: List from screen_genes().
  - `lmm_results`: List from run_lmm_analysis().
  - `term_pattern`: String. Which LMM terms to focus on.
- **Returns:** Data frame comparing both methods.

#### Function: `run_complete_analysis`
- **Description:** End-to-end analysis pipeline combining screening and LMM fitting.
- **Parameters:**
  - `seurat_obj`: Seurat object.
  - `config`: List from create_analysis_config().
  - `screening_params`: List of parameters for screen_genes().
  - `lmm_params`: List of parameters for run_lmm_analysis().
  - `check_parity`: Logical. Run parity diagnostics first.
  - `export_results`: Logical. Export results to Excel.
  - `output_prefix`: String. Prefix for output files.
- **Returns:** List with screening results, LMM results, and file paths.

### File: `myR/R/demulti_utils.R`

#### Function: `get_best_two`
- **Description:** This function takes a row of probability values and returns the indices of the two highest values.
- **Parameters:**
  - `row`: A numeric vector of probability values.
- **Returns:** A numeric vector of length 2 containing the indices of the two highest probability values.

#### Function: `get_barcode_mapping`
- **Description:** This function processes demultiplexing data to create a mapping between barcodes and their most likely sample assignments.
- **Parameters:**
  - `demux_data`: A data frame containing demultiplexing results.
- **Returns:** A data frame containing barcode mapping information.

#### Function: `is_doublet`
- **Description:** This function checks if a sample name contains a '+' character, indicating it is a doublet.
- **Parameters:**
  - `sample`: Character string containing the sample name.
- **Returns:** Logical value: TRUE if the sample is a doublet, FALSE otherwise.

#### Function: `generate_sample_values`
- **Description:** This function generates a vector of sample values including both singlets and doublets.
- **Parameters:**
  - `start_num`: Integer specifying the starting sample number.
  - `end_num`: Integer specifying the ending sample number.
- **Returns:** A character vector containing all possible sample combinations.

#### Function: `generate_sample_names`
- **Description:** This function generates a vector of sample names including both singlets and doublets.
- **Parameters:**
  - `vector`: Character vector containing sample names.
- **Returns:** A character vector containing all possible sample combinations.

### File: `myR/R/utils.R`

#### Function: `downsample_sobj`
- **Description:** This function randomly samples cells from a Seurat object to create a smaller subset.
- **Parameters:**
  - `sobj`: A Seurat object to be downsampled.
  - `ratio`: Integer indicating the downsampling ratio (1:ratio). Default is 10.
  - `seed`: Integer specifying the random seed for reproducibility. Default is 1234.
- **Returns:** A downsampled Seurat object.

#### Function: `printmy`
- **Description:** Prints marker genes in a comma-separated format.
- **Parameters:**
  - `markers`: Data frame with marker gene results.
  - `sign`: Sign of logFC to print ("+" or "-").
  - `num`: Number of genes to print (default: 100).
  - `pseudobulk`: Whether markers are from pseudobulk analysis.
- **Returns:** Prints genes to the console.

#### Function: `printMy`
- **Description:** Prints marker genes from a named list of marker data frames.
- **Parameters:**
  - `markers_list`: Named list of marker data frames.
  - `...`: Additional arguments passed to `printmy()`.
- **Returns:** Prints genes to the console.

#### Function: `sort_samples`
- **Description:** Sorts character sample identifiers, handling single-number and dual-number strings.
- **Parameters:**
  - `samples`: Character vector of sample identifiers.
- **Returns:** Sorted character vector of sample IDs.

#### Function: `print_gene_combinations`
- **Description:** For a list of gene sets, print genes that are unique to each combination.
- **Parameters:**
  - `gene_list`: Named list of character vectors.
  - `num_print`: Maximum number of genes to print per combination (default: 100).
- **Returns:** Prints gene combinations to the console.

#### Function: `sort_samples` (Second Definition)
- **Description:** Helper: Sort Sample Identifiers. Sorts character sample identifiers, handling single-number and dual-number strings (e.g., "1+2").
- **Parameters:**
  - `samples`: Character vector of sample identifiers.
- **Returns:** Sorted character vector of sample IDs.
