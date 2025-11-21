# Progress Report: 2025-11-21

## 1. CMGI Earth Package Fix
- **Issue**: `earth` package was missing in `renv` environment, causing CMGI to fail for `earth` models.
- **Action**:
    - Activated `renv` (`/home/user3/GJC_KDW_250721`).
    - Installed `earth` package using `renv::install("earth")`.
    - Verified fix with `test_cmgi_earth_fix.R`.
- **Result**: CMGI successfully extracts importance from `earth` models.
- **Merge**: Merged `cmgi` branch into `fgs` (manual conflict resolution to keep FGS improvements and apply CMGI updates).

## 2. FGS Failed Clusters Resolution
- **Issue**: Some clusters failed during FGS/TML7 pipeline.
- **Action**:
    - Created `run_failed_clusters_v2.R` to identify and process missing clusters.
    - Used `min_cells=3` to handle small clusters.
    - Directly sourced `myR/R/signature.R` to ensure latest code is used (overriding `start.R`'s loaded package).
- **Status**: Script is running. Processing 7 missing clusters out of 24 total.
    - Currently processing: `Monocytes / Macrophages`.

## 3. TML7 Development (Phase 1)
- **Objective**: Implement model-specific normalization for gene importance.
- **Action**:
    - Created `myR/R/signature_dev.R` with `compute_meta_gene_importance_v2`.
    - Implemented normalization methods: `max_abs`, `min_max`, `softmax`, `rank`, `z_score`.
    - Created `test_cmgi_dev.R` to verify implementation.
- **Result**:
    - Successfully tested on real TML data (`CD4+ T-cells`).
    - `softmax` and `z_score` provide distinct rankings compared to `max_abs`.
    - Ready for integration into main codebase.

## 4. Ranger vs Random Forest Comparison
- **Objective**: Verify consistency between `randomForest` and `ranger` implementations.
- **Action**:
    - Created `compare_ranger_rf_test.R`.
    - Running on `Naive B-cells` (1600 cells).
- **Status**: Running. `random_forest` step is taking time due to dataset size.

## Next Steps
- Monitor `run_failed_clusters_v2.R` progress.
- Analyze `compare_ranger_rf_test.R` results once complete.
- Integrate `compute_meta_gene_importance_v2` into `signature.R` (after current runs finish).
