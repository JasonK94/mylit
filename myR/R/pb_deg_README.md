# Pseudobulk Analysis Functions (pb_deg.R)

## Overview

This module provides a comprehensive suite of pseudobulk analysis functions for single-cell RNA-seq data. All functions follow a clear naming convention with the `pb_` prefix for easy identification and organization.

## Refactoring Summary

**Original file:** `pseudobulk_deg.R` (3,245 lines, 17 functions with 6 duplicates)
**Refactored file:** `pb_deg.R` (2,346 lines, 16 unique functions)
**Reduction:** ~28% fewer lines, 0 duplicates, improved modularity

---

## Naming Convention

All functions follow a structured naming pattern:

- `pb_utils_*`: Utility/helper functions
- `pb_prepare_*`: Data preparation functions
- `pb_run_*`: Statistical backend functions
- `pb_deg_*`: Differential expression analysis functions
- `pb_linear_*`: Linear modeling functions

---

## Function Categories

### 1. Utility Functions (`pb_utils_*`)

#### `pb_utils_validate_metadata(obj, required_cols, metadata_name = "metadata")`
Validates that required columns exist in Seurat object metadata.

**Parameters:**
- `obj`: Seurat object or data frame
- `required_cols`: Character vector of required column names
- `metadata_name`: Optional name for error messages

**Returns:** NULL (stops with error if validation fails)

---

#### `pb_utils_check_samples(meta, group_col, min_samples = 2, verbose = TRUE)`
Checks if each group has sufficient samples.

**Parameters:**
- `meta`: Metadata data frame
- `group_col`: Column name for grouping
- `min_samples`: Minimum required samples per group
- `verbose`: Logical, whether to print messages

**Returns:** Named logical vector indicating which groups have sufficient samples

---

#### `pb_utils_build_design(condition_col, covariates = NULL, custom_formula = NULL)`
Builds design formula for statistical models.

**Parameters:**
- `condition_col`: Main condition column name
- `covariates`: Optional vector of covariate column names
- `custom_formula`: Optional custom formula (overrides automatic building)

**Returns:** Formula object

---

#### `pb_utils_format_results(results, method = "edgeR", add_gene_col = TRUE)`
Standardizes DEG results to common format.

**Parameters:**
- `results`: Results data frame from DEG analysis
- `method`: Method used ("edgeR", "DESeq2", "wilcox", "t-test")
- `add_gene_col`: Logical, whether to add gene column from rownames

**Returns:** Standardized results data frame with columns: `gene`, `avg_log2FC`, `p_val`, `p_val_adj`

---

#### `pb_utils_create_binary_groups(meta, col, target_val, labels = c("GROUP1_TARGET", "GROUP2_REFERENCE"))`
Creates binary group assignments.

**Parameters:**
- `meta`: Metadata data frame
- `col`: Column name to check
- `target_val`: Target value for group 1
- `labels`: Character vector of length 2 for group labels

**Returns:** Character vector of group assignments

---

### 2. Preparation Functions (`pb_prepare_*`)

#### `pb_prepare_edgeR(seurat_obj, assay = "SCT", slot = "counts", sample_col, cluster_col, group_col, min_count = 10, norm_method = "TMM", design_formula = NULL, verbose = TRUE)`

**Purpose:** Converts single-cell Seurat object to pseudobulk data and prepares for edgeR analysis.

**Parameters:**
- `seurat_obj`: Seurat object
- `assay`: Assay name (default: "SCT")
- `slot`: Slot name (default: "counts")
- `sample_col`: Sample identifier column in metadata
- `cluster_col`: Cluster identifier column in metadata
- `group_col`: Group comparison column in metadata
- `min_count`: Minimum count for filterByExpr (default: 10)
- `norm_method`: Normalization method (default: "TMM")
- `design_formula`: Optional custom design formula
- `verbose`: Print progress messages

**Returns:** List containing:
- `pb`: Pseudobulk expression matrix (genes × sample_cluster)
- `meta`: Pseudobulk sample metadata
- `dge`: Filtered and normalized DGEList object
- `design`: Design matrix
- `contrast_levels`: Group levels for contrast setting

**Example:**
```r
prep_data <- pb_prepare_edgeR(
  seurat_obj = sobj,
  assay = "SCT",
  slot = "counts",
  sample_col = "exp_sample_no",
  cluster_col = "annotation2_big",
  group_col = "group3"
)
```

**Renamed from:** `prepare_pseudobulk_edgeR`

---

#### `pb_prepare_matrix(seurat_obj, ident.1, ident.2 = NULL, group.by = NULL, sample.by, aggregate.by = NULL, assay = NULL, slot = "counts", covariates = NULL, min.cells.per.sample = 10, verbose = TRUE)`

**Purpose:** Creates pseudobulk count matrix with flexible aggregation options.

**Parameters:**
- `seurat_obj`: Seurat object
- `ident.1`: First group identifier
- `ident.2`: Second group identifier (optional, for pairwise comparison)
- `group.by`: Grouping variable for subsetting
- `sample.by`: Column for sample aggregation (REQUIRED)
- `aggregate.by`: Additional aggregation variable (e.g., cell type)
- `assay`: Assay to use (NULL = default assay)
- `slot`: Data slot (default: "counts")
- `covariates`: Additional metadata columns to include
- `min.cells.per.sample`: Minimum cells per pseudobulk sample
- `verbose`: Print progress messages

**Returns:** List containing:
- `matrix`: Pseudobulk count matrix
- `metadata`: Associated sample metadata

**Example:**
```r
pb_data <- pb_prepare_matrix(
  seurat_obj = sobj,
  ident.1 = "TypeA",
  ident.2 = "TypeB",
  group.by = "cell_type",
  sample.by = "patient_id",
  min.cells.per.sample = 20
)
```

**Renamed from:** `create_pseudobulk_matrix_advanced`

---

### 3. Statistical Backend Functions (`pb_run_*`)

#### `pb_run_edgeR(pb_matrix, metadata, design.formula = NULL, covariates = NULL, ident.1, ident.2, contrast.type = "coef", logfc.threshold = 0, verbose = TRUE)`

**Purpose:** Executes edgeR statistical analysis on pseudobulk data.

**Parameters:**
- `pb_matrix`: Pseudobulk count matrix
- `metadata`: Sample metadata
- `design.formula`: Custom design formula (optional)
- `covariates`: Covariate column names
- `ident.1`: First group identifier
- `ident.2`: Second group identifier
- `contrast.type`: Type of contrast ("coef" or "contrast")
- `logfc.threshold`: Minimum |log2FC| to report
- `verbose`: Print messages

**Returns:** Data frame with DEG results (gene, avg_log2FC, p_val, p_val_adj, etc.)

**Example:**
```r
edgeR_results <- pb_run_edgeR(
  pb_matrix = pb_data$matrix,
  metadata = pb_data$metadata,
  ident.1 = "Treatment",
  ident.2 = "Control",
  logfc.threshold = 0.5
)
```

**Renamed from:** `run_edgeR_pseudobulk_advanced`

---

#### `pb_run_DESeq2(pb_matrix, metadata, design.formula = NULL, covariates = NULL, ident.1, ident.2, logfc.threshold = 0, verbose = TRUE)`

**Purpose:** Executes DESeq2 statistical analysis on pseudobulk data.

**Parameters:**
- `pb_matrix`: Pseudobulk count matrix
- `metadata`: Sample metadata
- `design.formula`: Custom design formula (optional)
- `covariates`: Covariate column names
- `ident.1`: First group identifier
- `ident.2`: Second group identifier
- `logfc.threshold`: Minimum |log2FC| to report
- `verbose`: Print messages

**Returns:** Data frame with DEG results

**Example:**
```r
deseq_results <- pb_run_DESeq2(
  pb_matrix = pb_data$matrix,
  metadata = pb_data$metadata,
  ident.1 = "Disease",
  ident.2 = "Healthy"
)
```

**Renamed from:** `run_DESeq2_pseudobulk_advanced`

---

### 4. DEG Analysis Functions (`pb_deg_*`)

#### `pb_deg_edgeR(analysis_level = c("overall", "per_cluster", "specific_cluster"), contrast, target_cluster = NULL, min_samples_per_group_cluster = 2, min_count = 10, ..., verbose = TRUE)`

**Purpose:** Main orchestrator for edgeR-based DEG analysis at multiple levels.

**Parameters:**
- `analysis_level`: Analysis scope
  - `"overall"`: All clusters combined
  - `"per_cluster"`: Each cluster separately
  - `"specific_cluster"`: Single specified cluster
- `contrast`: Contrast vector for glmQLFTest
- `target_cluster`: Target cluster (for specific_cluster mode)
- `min_samples_per_group_cluster`: Minimum samples per group
- `min_count`: Minimum count for filtering
- `...`: Arguments for `pb_prepare_edgeR` (seurat_obj, sample_col, etc.)
- `verbose`: Print progress

**Returns:** Tibble with DEG results (includes `cluster` column for per_cluster/specific_cluster)

**Example:**
```r
# Overall analysis
overall_deg <- pb_deg_edgeR(
  analysis_level = "overall",
  contrast = c(-1, 1),  # Treatment vs Control
  seurat_obj = sobj,
  sample_col = "patient_id",
  cluster_col = "cell_type",
  group_col = "condition"
)

# Per-cluster analysis
per_cluster_deg <- pb_deg_edgeR(
  analysis_level = "per_cluster",
  contrast = c(-1, 1),
  seurat_obj = sobj,
  sample_col = "patient_id",
  cluster_col = "cell_type",
  group_col = "condition"
)
```

**Renamed from:** `run_pseudobulk_deg`

---

#### `pb_deg_cluster(sobj, BULK_COL, SAMPLE_COL, GROUP_COL, COMPARISON_MODE = c("FindAllMarkers", "OneVsRest", "GroupVsGroup"), TEST_METHOD = c("wilcox", "t-test"), TARGET_CELLTYPE_OR_CLUSTER = NULL, TARGET_GROUP1 = NULL, TARGET_GROUP2 = NULL, MIN_SAMPLES_THRESHOLD = 2, MIN_CELL_COUNT = 10, MIN_PCT = 0.1, LOGFC_THRESHOLD = 0.25, edgeR.filterByExpr.min.count = 10, edgeR.calcNormFactors.method = "TMM", VERBOSE = TRUE)`

**Purpose:** Identifies cluster marker genes using pseudobulk approach.

**Comparison Modes:**
1. `"FindAllMarkers"`: Find markers for each cluster vs all other clusters
2. `"OneVsRest"`: Specific cluster vs all others
3. `"GroupVsGroup"`: Compare two specific groups within a cluster

**Parameters:**
- `sobj`: Seurat object
- `BULK_COL`: Column defining clusters/cell types
- `SAMPLE_COL`: Sample identifier column
- `GROUP_COL`: Group comparison column
- `COMPARISON_MODE`: Type of comparison
- `TEST_METHOD`: Statistical test ("wilcox" uses edgeR, "t-test" uses t-test)
- `TARGET_CELLTYPE_OR_CLUSTER`: Target cluster (for OneVsRest)
- `TARGET_GROUP1`, `TARGET_GROUP2`: Groups to compare (for GroupVsGroup)
- `MIN_SAMPLES_THRESHOLD`: Minimum samples required
- `MIN_CELL_COUNT`: Minimum cells per sample
- Other parameters for filtering and edgeR settings

**Returns:** Data frame with marker genes, log2FC, p-values, and expression percentages

**Example:**
```r
# Find all markers
all_markers <- pb_deg_cluster(
  sobj = sobj,
  BULK_COL = "cell_type",
  SAMPLE_COL = "patient_id",
  GROUP_COL = "condition",
  COMPARISON_MODE = "FindAllMarkers",
  TEST_METHOD = "wilcox"
)

# One vs rest
tcell_markers <- pb_deg_cluster(
  sobj = sobj,
  BULK_COL = "cell_type",
  SAMPLE_COL = "patient_id",
  GROUP_COL = "condition",
  COMPARISON_MODE = "OneVsRest",
  TARGET_CELLTYPE_OR_CLUSTER = "T_cells"
)
```

**Renamed from:** `cluster_pseudobulk_deg`

---

#### `pb_deg_flexible(sobj, bulk.category, sample.col, group.col = NULL, ident.1 = NULL, ident.2 = NULL, test.use = "wilcox", min.samples = 2, min.cells.per.sample = 10, edgeR.filterByExpr.min.count = 10, edgeR.calcNormFactors.method = "TMM", DESeq2.test = "Wald", logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE, verbose = TRUE)`

**Purpose:** Flexible DEG analysis supporting multiple statistical backends.

**Supported Methods:**
- `"wilcox"`: Wilcoxon rank-sum test (edgeR-based)
- `"t-test"`: t-test
- `"roc"`: ROC analysis
- `"edgeR"`: edgeR quasi-likelihood
- `"DESeq2"`: DESeq2 analysis

**Parameters:**
- `sobj`: Seurat object
- `bulk.category`: Column for aggregation (e.g., cell type)
- `sample.col`: Sample identifier column
- `group.col`: Group comparison column
- `ident.1`, `ident.2`: Groups to compare
- `test.use`: Statistical test method
- `min.samples`: Minimum samples per group
- `min.cells.per.sample`: Minimum cells per pseudobulk sample
- Other filtering and method-specific parameters

**Returns:** Data frame with DEG results for each bulk category

**Example:**
```r
# Using edgeR
edgeR_deg <- pb_deg_flexible(
  sobj = sobj,
  bulk.category = "cell_type",
  sample.col = "patient_id",
  group.col = "disease_status",
  ident.1 = "Disease",
  ident.2 = "Healthy",
  test.use = "edgeR"
)

# Using DESeq2
deseq2_deg <- pb_deg_flexible(
  sobj = sobj,
  bulk.category = "cell_type",
  sample.col = "patient_id",
  group.col = "treatment",
  ident.1 = "Treated",
  ident.2 = "Control",
  test.use = "DESeq2"
)
```

**Renamed from:** `pseudobulk_deg`

---

#### `pb_deg_FindMarkers(seurat_obj, ident.1, ident.2 = NULL, group.by = NULL, sample.by, aggregate.by = NULL, assay = NULL, slot = "counts", test.use = "edgeR", design.formula = NULL, covariates = NULL, min.cells.per.sample = 10, logfc.threshold = 0, contrast.type = "coef", verbose = TRUE)`

**Purpose:** Seurat-style FindMarkers using pseudobulk aggregation.

**Parameters:**
- Similar to Seurat's `FindMarkers` but with pseudobulk aggregation
- `test.use`: "edgeR" or "DESeq2"
- `sample.by`: REQUIRED for pseudobulk aggregation
- `aggregate.by`: Optional additional aggregation (e.g., by cell type)
- `design.formula`: Custom model formula
- `covariates`: Covariates to include in model

**Returns:** Data frame with DEG results

**Example:**
```r
markers <- pb_deg_FindMarkers(
  seurat_obj = sobj,
  ident.1 = "Activated",
  ident.2 = "Resting",
  group.by = "activation_status",
  sample.by = "patient_id",
  test.use = "edgeR",
  logfc.threshold = 0.5
)
```

**Renamed from:** `FindMarkers_pseudobulk`

---

#### `pb_deg_FindAllMarkers(seurat_obj, sample.by, aggregate.by = NULL, assay = NULL, slot = "counts", test.use = "edgeR", design.formula = NULL, covariates = NULL, min.cells.per.sample = 10, logfc.threshold = 0, only.pos = TRUE, verbose = TRUE)`

**Purpose:** Seurat-style FindAllMarkers using pseudobulk (finds markers for all identities).

**Parameters:**
- Similar to `pb_deg_FindMarkers` but iterates through all identities
- `only.pos`: Return only positive markers (default: TRUE)

**Returns:** Data frame with markers for all clusters/identities

**Example:**
```r
all_markers <- pb_deg_FindAllMarkers(
  seurat_obj = sobj,
  sample.by = "patient_id",
  aggregate.by = "cell_type",
  test.use = "edgeR",
  only.pos = TRUE,
  logfc.threshold = 0.5
)
```

**Renamed from:** `FindAllMarkers_pseudobulk`

---

### 5. Linear Modeling Functions (`pb_linear_*`)

#### `pb_linear_fit(sobj, genes, sample_col = "sample", numeric_predictor = "severity_score", group_col = NULL, p_adjust_method = "BH", min_samples_per_group = 3, min_distinct_predictor_values = 2)`

**Purpose:** Performs linear regression of gene expression vs continuous predictor.

**Use Cases:**
- Correlate expression with disease severity scores
- Analyze expression trends across continuous variables
- Group-stratified correlation analysis

**Parameters:**
- `sobj`: Seurat object
- `genes`: Vector of gene names to analyze
- `sample_col`: Sample identifier column
- `numeric_predictor`: Continuous predictor variable column
- `group_col`: Optional grouping column (for stratified analysis)
- `p_adjust_method`: Multiple testing correction method
- `min_samples_per_group`: Minimum samples required
- `min_distinct_predictor_values`: Minimum unique predictor values required

**Returns:** Data frame with:
- `gene`: Gene name
- `group`: Group name (if group_col specified)
- `intercept`: Regression intercept
- `slope`: Regression slope
- `slope_se`: Standard error of slope
- `n_samples_in_model`: Sample count
- `p_value`: P-value for slope
- `adj_p_value`: Adjusted p-value
- `r_squared`: R² value

**Example:**
```r
# Overall analysis
severity_fit <- pb_linear_fit(
  sobj = sobj,
  genes = c("IL6", "TNF", "IFNG"),
  sample_col = "patient_id",
  numeric_predictor = "disease_severity_score"
)

# Group-stratified analysis
stratified_fit <- pb_linear_fit(
  sobj = sobj,
  genes = top_genes,
  sample_col = "patient_id",
  numeric_predictor = "severity_score",
  group_col = "treatment_group"
)
```

**Renamed from:** `pseudobulk_linear_fit`

---

#### `pb_linear_compare_slopes(results_df, gene_col = "gene", group_col = "group", slope_col = "slope", se_col = "slope_se", n_samples_col = "n_samples_in_model", p_adjust_method = "BH", adjustment_scope = "global")`

**Purpose:** Post-hoc pairwise comparison of regression slopes between groups.

**Parameters:**
- `results_df`: Output from `pb_linear_fit`
- `gene_col`: Gene column name
- `group_col`: Group column name
- `slope_col`: Slope column name
- `se_col`: Standard error column name
- `n_samples_col`: Sample count column name
- `p_adjust_method`: Multiple testing correction
- `adjustment_scope`: "global" (all comparisons) or "per_gene" (within each gene)

**Returns:** Data frame with pairwise slope comparisons:
- `gene`: Gene name
- `group1`, `group2`: Groups being compared
- `slope_diff`: Difference in slopes
- `t_statistic`: t-statistic
- `df`: Degrees of freedom
- `p_value`: P-value
- `adj_p_value`: Adjusted p-value

**Example:**
```r
# Run linear fit on multiple groups
fit_results <- pb_linear_fit(
  sobj = sobj,
  genes = marker_genes,
  sample_col = "patient_id",
  numeric_predictor = "time_point",
  group_col = "treatment"
)

# Compare slopes between groups
slope_comparisons <- pb_linear_compare_slopes(
  results_df = fit_results,
  adjustment_scope = "per_gene"
)

# Find genes with significantly different slopes
sig_diff <- slope_comparisons %>%
  filter(adj_p_value < 0.05)
```

**Renamed from:** `post_hoc_slope_comparison`

---

## Migration Guide

### Old Function → New Function Mapping

| Old Function Name | New Function Name | Category |
|-------------------|-------------------|----------|
| `prepare_pseudobulk_edgeR` | `pb_prepare_edgeR` | Preparation |
| `create_pseudobulk_matrix_advanced` | `pb_prepare_matrix` | Preparation |
| `run_edgeR_pseudobulk_advanced` | `pb_run_edgeR` | Statistical Backend |
| `run_DESeq2_pseudobulk_advanced` | `pb_run_DESeq2` | Statistical Backend |
| `run_pseudobulk_deg` | `pb_deg_edgeR` | DEG Analysis |
| `cluster_pseudobulk_deg` | `pb_deg_cluster` | DEG Analysis |
| `pseudobulk_deg` | `pb_deg_flexible` | DEG Analysis |
| `FindMarkers_pseudobulk` | `pb_deg_FindMarkers` | DEG Analysis |
| `FindAllMarkers_pseudobulk` | `pb_deg_FindAllMarkers` | DEG Analysis |
| `pseudobulk_linear_fit` | `pb_linear_fit` | Linear Modeling |
| `post_hoc_slope_comparison` | `pb_linear_compare_slopes` | Linear Modeling |

### Example Migration

**Before:**
```r
prep <- prepare_pseudobulk_edgeR(sobj, sample_col = "patient", ...)
results <- run_pseudobulk_deg(analysis_level = "overall", ...)
```

**After:**
```r
prep <- pb_prepare_edgeR(sobj, sample_col = "patient", ...)
results <- pb_deg_edgeR(analysis_level = "overall", ...)
```

---

## Workflow Examples

### Example 1: Basic DEG Analysis

```r
library(Seurat)
library(dplyr)

# Prepare data
prep_data <- pb_prepare_edgeR(
  seurat_obj = sobj,
  assay = "RNA",
  sample_col = "patient_id",
  cluster_col = "cell_type",
  group_col = "condition"
)

# Run overall DEG
deg_results <- pb_deg_edgeR(
  analysis_level = "overall",
  contrast = c(-1, 1),  # Treatment vs Control
  prepared_data = prep_data
)

# Filter significant genes
sig_genes <- deg_results %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
```

### Example 2: Per-Cluster Analysis with FindMarkers Style

```r
# Find markers for specific cell type
tcell_markers <- pb_deg_FindMarkers(
  seurat_obj = sobj,
  ident.1 = "CD4_T",
  ident.2 = "CD8_T",
  group.by = "cell_type",
  sample.by = "patient_id",
  test.use = "edgeR",
  logfc.threshold = 0.5
)

# Find markers for all cell types
all_cell_markers <- pb_deg_FindAllMarkers(
  seurat_obj = sobj,
  sample.by = "patient_id",
  aggregate.by = "cell_type",
  test.use = "edgeR",
  only.pos = TRUE
)
```

### Example 3: Linear Modeling with Continuous Variables

```r
# Define genes of interest
genes_of_interest <- c("IL6", "TNF", "IFNG", "IL1B")

# Fit linear models
linear_results <- pb_linear_fit(
  sobj = sobj,
  genes = genes_of_interest,
  sample_col = "patient_id",
  numeric_predictor = "disease_severity",
  group_col = "treatment_group"
)

# Compare slopes between treatment groups
slope_comparison <- pb_linear_compare_slopes(
  results_df = linear_results,
  adjustment_scope = "per_gene"
)

# Visualize results
library(ggplot2)
sig_slopes <- linear_results %>% filter(adj_p_value < 0.05)
ggplot(sig_slopes, aes(x = gene, y = slope, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal()
```

### Example 4: Flexible Multi-Method Comparison

```r
# Compare results from different methods
results_edgeR <- pb_deg_flexible(
  sobj = sobj,
  bulk.category = "cell_type",
  sample.col = "patient_id",
  group.col = "condition",
  ident.1 = "Disease",
  ident.2 = "Healthy",
  test.use = "edgeR"
)

results_deseq2 <- pb_deg_flexible(
  sobj = sobj,
  bulk.category = "cell_type",
  sample.col = "patient_id",
  group.col = "condition",
  ident.1 = "Disease",
  ident.2 = "Healthy",
  test.use = "DESeq2"
)

# Compare concordance
library(dplyr)
comparison <- inner_join(
  results_edgeR %>% select(bulk_cat, gene, edgeR_logFC = avg_log2FC, edgeR_padj = p_val_adj),
  results_deseq2 %>% select(bulk_cat, gene, DESeq2_logFC = avg_log2FC, DESeq2_padj = p_val_adj),
  by = c("bulk_cat", "gene")
)
```

---

## Best Practices

### 1. Sample Size Considerations
- Ensure minimum 2-3 biological replicates per group
- Use `min_samples_per_group_cluster` parameter to enforce thresholds
- Check sample sufficiency with `pb_utils_check_samples()`

### 2. Choosing Statistical Methods
- **edgeR**: General-purpose, handles low counts well, fast
- **DESeq2**: More conservative, good for small sample sizes
- **Linear models**: For continuous predictors (severity scores, time points)

### 3. Metadata Validation
- Always validate required columns with `pb_utils_validate_metadata()`
- Ensure factor levels are correctly set for contrasts
- Check for missing values in grouping variables

### 4. Quality Control
- Filter low-quality cells before pseudobulk aggregation
- Use `min.cells.per.sample` to ensure robust aggregation
- Check filtering statistics in verbose output

### 5. Multiple Testing Correction
- Always use adjusted p-values (`p_val_adj`)
- Consider per-gene adjustment for post-hoc tests
- Be aware of test multiplicity in per-cluster analyses

---

## Dependencies

Required R packages:
- `Seurat` (≥ 4.0)
- `edgeR`
- `DESeq2`
- `dplyr`
- `tidyr`
- `tibble`
- `rlang`

---

## Performance Notes

- **Pseudobulk aggregation** is much faster than cell-level analysis
- **edgeR** is generally faster than DESeq2 for large datasets
- Use `verbose = FALSE` to suppress progress messages in batch processing
- Consider parallel processing for per-cluster analyses with many clusters

---

## Troubleshooting

### Common Issues

**Error: "다음 컬럼이 메타데이터에 없습니다"**
- Check column names in `sobj@meta.data`
- Ensure no typos in column name arguments

**Warning: "그룹별 최소 샘플 수를 만족하지 못합니다"**
- Increase sample size or reduce `min_samples_per_group_cluster`
- Check for groups with insufficient replicates

**Error: "filterByExpr 결과 남은 유전자가 없습니다"**
- Reduce `min_count` parameter
- Check if expression data is properly normalized
- Ensure count data (not log-normalized) is used

**NA values in group mapping**
- Check for numeric-like sample IDs (may need 'g' prefix handling)
- Verify sample IDs match between metadata and expression data
- Check for samples with multiple group assignments

---

## Citation

If you use these functions in your research, please cite:

- **edgeR**: Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics, 26(1), 139-140.

- **DESeq2**: Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550.

- **Seurat**: Hao et al., (2021). "Integrated analysis of multimodal single-cell data." Cell, 184, 3573-3587.

---

## Version History

- **v1.0** (2025-10-31): Initial refactored release
  - Removed 6 duplicate functions
  - Standardized naming convention (pb_* prefix)
  - Added utility functions for common patterns
  - Improved documentation and examples
  - 28% code reduction while maintaining all functionality

---

## Contact & Support

For issues, questions, or contributions, please refer to the main package documentation.
