# Differential Expression Analysis Module

## Overview

This module provides a unified interface for differential expression analysis across multiple methods. Instead of having separate files for each method, all DE methods are consolidated here for easier maintenance and discovery.

## Available Methods

### 1. Pseudobulk Differential Expression (edgeR)
Best for comparing groups with biological replicates at the sample level.

```r
# Overall pseudobulk DE
results <- run_pseudobulk_deg(
  seurat_obj,
  sample_col = "sample_id",
  group_col = "condition",
  comparison = c("Treatment", "Control"),
  mode = "overall"
)

# Per-cluster analysis
results_list <- run_pseudobulk_deg(
  seurat_obj,
  sample_col = "sample_id",
  group_col = "condition",
  comparison = c("Treatment", "Control"),
  cluster_col = "seurat_clusters",
  mode = "per_cluster"
)
```

### 2. Linear Regression
Best for modeling gene expression against continuous, categorical, or ordinal predictors.

```r
# Continuous predictor (e.g., age)
results <- linear_seurat(
  seurat_obj,
  regressor = "age",
  regressor.type = "continuous",
  layer = "data"
)

# Categorical predictor with covariates
results <- linear_seurat(
  seurat_obj,
  regressor = "cell_type",
  regressor.type = "categorical",
  reference.level = "Control",
  covariates = c("batch", "sex")
)

# Ordinal predictor (e.g., disease stage)
results <- linear_seurat(
  seurat_obj,
  regressor = "disease_stage",
  regressor.type = "ordinal",
  ordinal.method = "linear"
)
```

### 3. Linear Mixed Models (LMM)
(To be added from test.R)
Best for handling random effects like patient ID or batch effects.

```r
# Coming soon:
# results <- run_lmm_deg(...)
```

## Method Comparison

| Method | Use Case | Pros | Cons |
|--------|----------|------|------|
| **Pseudobulk (edgeR)** | Group comparisons with replicates | Gold standard for bulk-like DE, handles overdispersion well | Requires sample-level replicates |
| **Linear Regression** | Continuous/ordinal predictors, multi-factor designs | Flexible, handles various predictor types | Doesn't account for overdispersion |
| **LMM** | Repeated measures, hierarchical data | Handles random effects, paired/matched designs | Computationally intensive |

## When to Use Which Method

### Use Pseudobulk when:
- You have biological replicates (multiple samples per condition)
- You want to compare two or more groups
- Your data has typical scRNA-seq characteristics (count data, overdispersion)
- You want robust, publication-quality DE results

### Use Linear Regression when:
- You want to model gene expression against continuous variables (age, pseudotime, etc.)
- You have ordinal predictors (disease stage: mild/moderate/severe)
- You need to control for multiple covariates
- You're exploring relationships rather than strict comparisons

### Use LMM when:
- You have repeated measures (same patient at multiple timepoints)
- You have hierarchical/nested data structure
- You need to account for random effects (patient-to-patient variation)
- You have paired or matched samples

## File Organization

```
differential_expression/
├── differential_expression.R   # Main unified module
├── README.md                   # This file
└── (future additions)
    └── advanced_methods.R      # MAST, DESeq2, etc. (if needed)
```

## Migration Notes

### From standalone pseudobulk_deg.R:
```r
# Old way (still works via legacy file):
source("myR/R/pseudobulk_deg.R")
results <- run_pseudobulk_deg(...)

# New way (recommended):
source("myR/R/analysis/differential_expression/differential_expression.R")
results <- run_pseudobulk_deg(...)

# Or with package:
library(myR)
results <- run_pseudobulk_deg(...)
```

### From signature.R linear_seurat:
```r
# Old way:
source("myR/R/signature.R")
results <- linear_seurat(...)

# New way:
source("myR/R/analysis/differential_expression/differential_expression.R")
results <- linear_seurat(...)
```

## Contributing

When adding new DE methods:
1. Add the function to `differential_expression.R`
2. Follow the existing naming convention
3. Document with roxygen2
4. Add example to this README
5. Update the comparison table

## See Also

- **Marker Analysis**: `analysis/markers/` - for FindMarkers-based analysis
- **Pathway Analysis**: `analysis/pathway/` - for functional enrichment of DE results
- **GeoMx DE**: `analysis/spatial/geomx_analysis.R` - for spatial transcriptomics DE

