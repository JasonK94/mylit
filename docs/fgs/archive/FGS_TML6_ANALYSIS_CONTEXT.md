# FGS and TML6 Analysis Context

## Overview

This document provides context for running FGS (Find Gene Signature) and TML6 (Train Meta-Learner v6) analyses on the IS6 dataset.

## Data Information

### Dataset
- **File**: `/data/user3/sobj/IS6_sex_added_251110.qs`
- **Type**: Seurat object
- **Description**: Single-cell RNA-seq data with sex information added

### Key Variables

#### Main Fixed Effect
- **`g3`**: Target group variable
  - Values: `"1"`, `"2"`, `NA`
  - **IMPORTANT**: This is a **character/factor** variable, NOT numeric
  - `NA` values must be removed before analysis
  - Represents patient groups (e.g., good vs. bad prognosis)

#### Covariates
- **`sex`**: Biological cofactor for sex differences
  - Used to account for sex-specific effects
  
- **`anno3.scvi`**: Cluster annotation (scvi integrated)
  - Used to account for cluster-specific differences
  - Each cell belongs to a specific cluster

#### Batch/Random Effects
- **`GEM`**: Batch variable (maximum 8 levels)
  - Patients are nested within specific GEMs
  - Can be used as fixed effect (since FGS doesn't support random effects)
  
- **`hos_no`**: Patient ID (23 patients)
  - Patients are nested within specific groups (g3)
  - **Perfect confounding**: g3 and hos_no are perfectly confounded (each patient has one g3 value)
  - **Perfect nesting**: GEM completely contains hos_no (each patient belongs to one GEM)

#### Confounding Relationships
- **g3 ↔ hos_no**: Perfect confounding (each patient has one g3 value)
- **GEM ↔ hos_no**: Perfect nesting (each patient belongs to one GEM)
- **g3 ↔ GEM**: Some GEMs are perfectly separated by g3 values (some GEMs contain only g3=1 or only g3=2)

### Variable Keys
```r
sample_key <- "hos_no"      # Sample/patient ID
batch_key <- "GEM"          # Batch variable
group_key <- "g3"           # Target group variable
covar_key <- "sex"          # Covariate for sex
cluster_key <- "anno3.scvi" # Cluster annotation
```

## Model Specifications

### Ideal Model (Most Sophisticated)
```
formula1 <- ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)
```

**Components:**
- `g3`: Main fixed effect (target)
- `sex`: Sex covariate
- `anno3.scvi`: Cluster covariate
- `GEM`: Batch as fixed effect
- `g3:anno3.scvi`: Interaction between group and cluster
- `sex:anno3.scvi`: Interaction between sex and cluster
- `(1|GEM/hos_no)`: Nested random effects (GEM and patient within GEM)

### Fallback Model (FGS Compatible)
Since FGS and many other functions don't support:
- Random effects `(1|...)`
- Complex interactions in all methods

**Recommended fallback:**
```
Fixed effects only: ~ g3 + sex + anno3.scvi + GEM
```

**Note**: 
- GEM is used as a fixed effect (not random) because FGS doesn't support random effects
- Interactions may not be supported by all methods
- Some methods may only support simple additive models

## FGS Function

### Function Name
- **Primary**: `FGS()` (alias for `find_gene_signature_v5.3`)
- **Deprecated**: `FGS_v5.2()`, `FGS_v5.3()` (moved to test_to_delete.R)

### Key Features
- Supports multiple methods:
  - Tree-based: `random_forest`, `random_forest_ranger`, `xgboost`
  - Regularization: `lasso`, `ridge`, `elastic_net`
  - Dimensionality reduction: `pca_loadings`, `nmf_loadings`
  - Statistical modeling: `gam`, `limma`, `wilcoxon`
- Enhanced GAM with dynamic k adjustment (v5.3)
- Supports both single-cell and pseudobulk data

### Important Parameters
- `target_var`: "g3" (main fixed effect)
- `control_vars`: c("sex", "anno3.scvi", "GEM") (covariates)
- `method`: Vector of methods to use
- `n_features`: Number of top genes to return (default: 50)
- `gam.k`: NULL for dynamic k adjustment (v5.3 feature)
- `gam.k_dynamic_factor`: Factor for dynamic k calculation (default: 5)

### Data Cleaning Requirements
1. **Remove NA from g3**: Essential before analysis
   ```r
   is6_clean <- is6[, !is.na(is6@meta.data$g3)]
   ```

2. **Convert g3 to factor**: Prevents numeric confusion
   ```r
   is6_clean@meta.data$g3 <- factor(is6_clean@meta.data$g3, levels = c("1", "2"))
   ```

## TML6 Function

### Function Name
- `TML6()`: Train Meta-Learner v6

### Purpose
- Combines multiple L1 (Level-1) gene signatures into a unified prediction model
- Uses stacked ensemble approach with L2 (Level-2) models
- Selects best model via cross-validation

### Input Requirements
- `l1_signatures`: Named list of L1 signatures (from FGS results)
  - Each signature can be:
    - Character vector of gene names (uniform weights = 1)
    - Named numeric vector (gene names → weights)
    - List with `up` and/or `down` components
    - Data frame with gene and weight columns
- `holdout_data`: Expression data (Seurat object or matrix)
- `target_var`: Target variable name (for Seurat) or vector (for matrix)

### Key Parameters
- `l2_methods`: Vector of caret model methods (default: c("glm", "ranger", "xgbTree"))
- `k_folds`: Number of CV folds (default: 5)
- `metric`: Performance metric ("AUC", "ROC", "Accuracy", "Kappa")
- `allow_parallel`: Enable parallel execution (default: FALSE)

## Analysis Workflow

### Step 1: Data Loading and Cleaning
```r
# Load data
is6 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# Remove NA from g3
is6_clean <- is6[, !is.na(is6@meta.data$g3)]

# Convert g3 to factor
is6_clean@meta.data$g3 <- factor(is6_clean@meta.data$g3, levels = c("1", "2"))
```

### Step 2: FGS Analysis
```r
fgs_results <- FGS(
  data = is6_clean,
  target_var = "g3",
  control_vars = c("sex", "anno3.scvi", "GEM"),
  method = c("lasso", "random_forest", "limma", "gam"),
  n_features = 100,
  fgs_seed = 42,
  gam.k = NULL  # Use dynamic k
)
```

### Step 3: Extract Signatures for TML6
```r
l1_signatures <- list()
for (method_name in names(fgs_results)) {
  if (!is.null(fgs_results[[method_name]]$genes)) {
    genes <- fgs_results[[method_name]]$genes
    weights <- fgs_results[[method_name]]$weights
    sig_vec <- weights
    names(sig_vec) <- genes
    l1_signatures[[method_name]] <- sig_vec
  }
}
```

### Step 4: TML6 Meta-Learner
```r
meta_model <- TML6(
  l1_signatures = l1_signatures,
  holdout_data = is6_clean,
  target_var = "g3",
  l2_methods = c("glm", "ranger", "xgbTree"),
  k_folds = 5,
  metric = "AUC",
  fgs_seed = 42
)
```

### Step 5: Gene-Level Importance
```r
gene_importance <- compute_meta_gene_importance(meta_model, normalize = TRUE)
```

## Output Files

### FGS Results
- **File**: `/data/user3/sobj/FGS_results_IS6_251110.qs`
- **Content**: List of method-specific results
- **Structure**: Each method contains `genes`, `weights`, `scores`, `performance`, `method`, etc.

### TML6 Results
- **File**: `/data/user3/sobj/TML6_results_IS6_251110.qs`
- **Content**: Meta-learner model and training data
- **Structure**: 
  - `best_model`: Best L2 model (caret train object)
  - `best_model_name`: Name of best model
  - `best_metric_name`: Metric used for selection
  - `model_comparison`: Comparison of all models
  - `l2_train`: Training data with signature scores
  - `l1_signatures`: Standardized signatures

### Gene Importance
- **File**: `/data/user3/sobj/gene_importance_IS6_251110.qs`
- **Content**: Gene-level importance scores
- **Structure**:
  - `signature_importance`: L1 signature importances
  - `gene_importance`: Per-gene contributions
  - `gene_summary`: Aggregated contributions per gene
  - `positive_class`: Positive class name
  - `model_type`: Selected L2 model type

## Important Notes

### Data Type Warnings
1. **g3 is NOT numeric**: Always treat as character/factor
2. **NA handling**: Must remove NA from g3 before analysis
3. **Confounding**: Be aware of perfect confounding between g3 and hos_no

### Model Limitations
1. **No random effects**: FGS doesn't support `(1|...)` syntax
2. **Limited interactions**: Not all methods support interactions
3. **Fixed effects only**: GEM used as fixed effect, not random

### Method Selection
- Start with subset of methods for faster execution
- Full method list available but computationally expensive
- GAM with dynamic k (v5.3) recommended for better performance

### Performance Considerations
- FGS: Minutes to hours depending on methods and data size
- TML6: Additional time for cross-validation
- Parallel execution available but opt-in (`allow_parallel=TRUE`)

## Example Usage

See `scripts/run_FGS_TML6_analysis.R` for complete working example.

## References

- FGS documentation: `myR/R/signature.R` (FGS function)
- TML6 documentation: `myR/R/signature.R` (TML6 function)
- Test instructions: `docs/deg_consensus/TEST_INSTRUCTIONS_example.md`

