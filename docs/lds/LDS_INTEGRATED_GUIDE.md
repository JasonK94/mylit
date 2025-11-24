# Limma-Dream-SVA (LDS) Pipeline

## Overview

LDS (Limma-Dream-SVA) is an integrated pipeline for fitting Linear Mixed Models (LMM) with multiple random effects (e.g., patient, batch) in single-cell RNA-seq or GeoMx Digital Spatial Profiling data. It detects hidden covariates through SVA (Surrogate Variable Analysis) and performs differential expression analysis using limma's `dream` function.

**Key Features**:
- Suitable for cases with relatively few samples, such as GeoMx data
- Multiple random effects modeling (e.g., `(1|patient) + (1|batch)`)
- Automatic detection and correction of hidden variability through SVA
- Automatic determination of SV count based on residual variance

## Core Concepts

### SVA (Surrogate Variable Analysis)

SVA is a method for detecting unknown covariates (e.g., technical batch effects, hidden biological factors).

#### Mathematical Principles of SVA

SVA finds principal components that explain **residual variance**:

1. **Fixed Effects Model Fitting**: Fit model with user-specified fixed effects (e.g., `g3_clean + (1|hos_no) + (1|GEM)`)
   ```
   Y (genes × cells) = X (cells × covariates) × β + R (genes × cells)
   ```
   where R is the residual matrix.

2. **Residual Calculation**: Calculate residuals after removing fixed effects
   ```
   R = Y - X × β^T
   ```
   Calculate residuals after fitting regression model for each gene

3. **SVD of Residual Matrix**: Analyze the structure of the entire residual matrix
   ```
   R = U × D × V^T
   ```
   - U: Principal components in gene space (genes × SV)
   - D: Singular values (importance of SV)
   - V: Principal components in cell space (cells × SV) ← This is SV!

4. **SV Selection**: Calculate the proportion of residual variance explained by each SV
   ```
   Variance explained ratio = (D[i]^2) / sum(D^2)
   ```
   - Select minimum number of SVs where cumulative variance exceeds `sv_var_cutoff`
   - Or if `sv_var_cutoff = NULL`, use all significant SVs

5. **Model Correction**: Add extracted SVs as covariates to correct the model

**Important Points**:
- SVA does not explain the **total residual variance of all 291 genes**
- It explains **residual variance not explained by fixed effects**
- This residual variance may include batch effects, hidden covariates, etc.

#### Criteria for Significant SVs

The `sva()` function in the SVA package automatically determines this when `n.sv = NULL`:

- **Permutation Test Based**: Detects statistically significant SVs by comparing residuals of `mod` (full model) and `mod0` (null model)
- **Residual Variance Based**: Selects SVs in order of large singular values from SVD of residual matrix
- Generally determines significant SVs based on **p-value < 0.05**
- Additional filtering with `sv_var_cutoff` determines the actual number to use

### Limma-Dream

`dream` is an extension of the limma package that can fit LMMs with multiple random effects:

- **voomWithDreamWeights**: voom transformation considering weights
- **dream**: LMM fitting (supports lme4 formula syntax)
- **eBayes**: Empirical Bayes adjustment

### Pipeline Workflow

```
Input Data (Seurat/DGEList/Matrix)
    ↓
[1] Data Extraction and Validation
    ↓
[2] Formula Parsing (Extract Fixed Effects)
    ↓
[3] DGEList Creation, Filtering, Normalization
    ↓
[4] SVA Execution (after voom transformation)
    - Detect maximum SV count
    - Automatically determine SV count based on residual variance (or user-specified)
    ↓
[5] Final Formula Generation (original + SV)
    ↓
[6] limma-dream Pipeline
    - voomWithDreamWeights
    - dream (LMM fitting)
    - eBayes (Empirical Bayes)
    ↓
Final Result (MArrayLM object)
```

## Main Functions

### `LDS()`

Main function that executes the entire LDS pipeline.

#### Parameters

| Parameter | Description | Default |
|-----------|------------|---------|
| `sobj` | Seurat object, DGEList object, or Raw Count Matrix | Required |
| `formula` | Model formula (lme4 syntax) | Required |
| `meta.data` | Metadata (required for Matrix input) | `NULL` |
| `layer` | Count layer of Seurat object | `"counts"` |
| `n_sv` | Number of SVs to use | `NULL` (automatic) |
| `sv_var_cutoff` | Proportion of residual variance SVs should explain. If `NULL`, use all significant SVs | `0.5` |
| `n_cores` | Number of CPU cores for parallel processing | `max(1, detectCores()-2)` |
| `remove_na` | Whether to filter NA values | `TRUE` |
| `min.count` | `filterByExpr`'s `min.count` | `10` |
| `min.total.count` | `filterByExpr`'s `min.total.count` | `15` |
| `min.prop` | `filterByExpr`'s `min.prop` | `0.1` |
| `large.n` | `filterByExpr`'s `large.n` | `10` |
| `plot_sva_correlation` | Whether to perform SVA correlation analysis and generate Heatmaps | `TRUE` |

#### Formula Examples

```r
# Basic example: treatment effect, patient and batch as random effects
~ treatment + (1|patient) + (1|batch)

# Complex example: multiple fixed and random effects
~ response + age + sex + (1|emrid) + (1|set)

# Random effects only (no fixed effects)
~ (1|patient) + (1|batch)
```

#### Return Value

List object:
- `fit`: MArrayLM object (can use limma's `topTable()`, includes p-value)
- `voom`: EList object with weights
- `sva_obj`: Original SVA object
- `svs_used`: SV matrix actually used in the model
- `final_formula`: Final formula used (string)
- `dge`: Filtered/normalized DGEList
- `fixed_vars`: List of fixed effect variables (used as coef in topTable)
- `contrast_applied`: Whether contrast was applied
- `sv_correlation`: Correlation matrix between SVA and metadata (optional)
- `heatmaps`: Generated Heatmap file paths and correlation matrix (optional)

#### Usage Example

```r
# Run on Seurat object
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient) + (1|batch),
  n_sv = NULL,  # automatic determination
  sv_var_cutoff = 0.5
)

# Check results
top_genes <- limma::topTable(result$fit, number = 100)
head(top_genes)

# Check used SVs
head(result$svs_used)
result$final_formula
```

## SV Count Determination Methods

### 1. User-Specified (provide `n_sv`)

```r
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient),
  n_sv = 5  # use exactly 5 SVs
)
```

### 2. Automatic Determination (`n_sv = NULL`, default)

#### 2-1. Specify `sv_var_cutoff` (default 0.5)

Automatically determined based on cumulative proportion of residual variance:

1. Detect maximum SV count with SVA
2. Calculate proportion of residual variance explained by each SV
3. Select minimum number of SVs where cumulative variance exceeds `sv_var_cutoff` (default 50%)

```r
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient),
  n_sv = NULL,  # automatic determination
  sv_var_cutoff = 0.5  # explain 50% of residual variance
)
```

**Example**: If SV 1 explains 30% and SV 2 explains 55% → select 2

#### 2-2. `sv_var_cutoff = NULL` (use all SVs)

Use all significant SVs:

```r
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient),
  n_sv = NULL,
  sv_var_cutoff = NULL  # use all significant SVs
)
```

## Result Interpretation

### Using `topTable()`

```r
# Check fixed effect variables
result$fixed_vars

# Calculate p-value by specifying coef
all_results <- limma::topTable(
  result$fit,
  coef = result$fixed_vars[1],  # first fixed effect variable
  number = Inf,
  sort.by = "P"
)

# Or without coef (if p-value already calculated)
all_results <- limma::topTable(
  result$fit,
  number = Inf,
  sort.by = "P"
)
```

**Note**: For `dream` results, the `LDS()` function automatically calculates p-values internally using t-statistics and df.

### Key Columns

- `logFC`: log fold change
- `AveExpr`: average expression
- `t`: t-statistic
- `P.Value`: p-value
- `adj.P.Val`: FDR-adjusted p-value
- `B`: log-odds

## Warnings

### 1. Data Format

- **Seurat Object**: Extract raw counts with `layer = "counts"`
- **DGEList**: Use `$counts` and `$samples`
- **Matrix**: Must provide `meta.data`

### 2. Formula Writing

- Use lme4 formula syntax: `(1|group)` format
- Distinguish between fixed and random effects
- Cases with only random effects (`~ (1|patient)`) are possible, but SVA's `mod0` is set to `~ 1`

### 3. Sample Size

- Suitable for cases with relatively few samples, such as GeoMx data
- SV detection may be difficult with too few samples

### 4. Package Dependencies

Required packages:
- `limma` (includes dream, voomWithDreamWeights)
- `edgeR` (DGEList, filterByExpr, calcNormFactors, voom)
- `lme4` (formula parsing)
- `BiocParallel` (parallel processing)
- `sva` (SVA execution)

### 5. Memory and Execution Time

- Large datasets may use significant memory
- SVA and dream are computationally intensive
- Parallel processing speed can be improved by adjusting `n_cores`

## Troubleshooting

### SVA Fails to Find SVs

```
... SVA failed to find significant surrogate variables (SVs).
```

- Sample size may be too small
- Fixed effects model may already explain most variability
- May be normal (proceed without SVs)

### All Genes Filtered

```
All genes were filtered.
```

- `filterByExpr` conditions may be too strict
- Data quality issues
- Need to check `design_for_filter`

### dream Fitting Failure

- Check formula syntax
- Verify all required variables are in metadata
- Check if sample size is too small

## SVA Correlation Analysis and Visualization

### Automatic Heatmap Generation

When `plot_sva_correlation = TRUE` (default), the LDS function automatically generates 3 types of Heatmaps:

1. **Full Correlation Matrix** (`{prefix}_heatmap_full.png`)
   - Full correlation matrix including metadata and SVs
   - Size: (metadata + SV) × (metadata + SV)

2. **Metadata × SV** (`{prefix}_heatmap_metadata_x_sv.png`)
   - Correlations between all metadata variables and SVs
   - Rows: metadata variables, Columns: SVs

3. **Top n Metadata × SV** (`{prefix}_heatmap_top10_metadata_x_sv.png`)
   - Top 10 metadata variables with p-value < 0.05
   - Displays p-values with asterisks (* p<0.05, ** p<0.01, *** p<0.001)

### Helper Functions

#### `lds_corrplot()`

Independent function that generates correlation plots between SVA and metadata.

**Difference**: `lds_08_create_heatmaps()` runs automatically as part of the LDS pipeline, but `lds_corrplot()` is a helper function that can be used independently.

```r
# Use independently
lds_corrplot(
  sva = result$svs_used,
  meta.data = meta.data,
  method = "spearman",
  display.value = TRUE,
  output_file = "correlation_plot.png"
)
```

#### `lds_corrplot_summary()`

Returns a correlation summary table:

```r
summary_table <- lds_corrplot_summary(
  sva = result$svs_used,
  meta.data = meta.data,
  top_n = 20
)
```

### File Locations

Generated Heatmap files are saved in `output_dir`:
- Default: `tempdir()` (R session's temporary directory)
- When `save_intermediate = TRUE`: specified `output_dir`

```r
# Check file paths from results
result$heatmaps$heatmap_files
```

## References

- limma package documentation: https://bioconductor.org/packages/limma
- dream function: `?limma::dream`
- SVA paper: Leek et al. (2007) Nature Reviews Genetics
- Modularized implementation: `myR/R/lds.R`
- Helper functions: `myR/R/lds_corrplot.R`, `myR/R/lds_08_heatmaps.R`

