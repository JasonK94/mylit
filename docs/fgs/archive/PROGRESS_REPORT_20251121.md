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
    - Directly sourced `myR/R/signature.R` to ensure latest code is used.
    - **Update**: Excluded `random_forest` method (too slow) and focused on `ranger` and `nmf`.
- **Status**: Script is running (`run_failed_clusters_v4.log`).

## 3. TML7 Development (Phase 1)
- **Objective**: Implement model-specific normalization for gene importance.
- **Action**:
    - Created `myR/R/signature_dev.R` with `compute_meta_gene_importance_v2`.
    - Implemented normalization methods: `max_abs`, `min_max`, `softmax`, `rank`, `z_score`.
    - Verified with `test_cmgi_dev.R`.
    - **Integration**: Integrated the enhanced `compute_meta_gene_importance` into `myR/R/signature.R`.
    - **Documentation**: Created `docs/tml/MODEL_SPECIFIC_NORMALIZATION.md`.
- **Result**:
    - Successfully tested on real TML data (`CD4+ T-cells`).
    - Codebase now supports advanced normalization for gene importance.

## 4. Ranger vs Random Forest Comparison
- **Objective**: Verify consistency between `randomForest` and `ranger` implementations.
- **Action**:
    - Created `compare_ranger_rf_test.R`.
    - Downsampled to 300 cells for speed.
    - Fixed `intersect` conflict issue.
- **Result**:
    - **Jaccard Index: 1.0000**
    - `random_forest` and `random_forest_ranger` produced identical gene lists (Top 50 features).
    - Confirms `random_forest_ranger` implementation is correct and can replace the slower `random_forest`.

## Next Steps
- Monitor `run_failed_clusters_v2.R` completion.
- Finalize and push changes.
