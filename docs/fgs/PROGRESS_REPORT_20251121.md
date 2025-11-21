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
    - Initial run confirmed FGS execution success (verifying `random_forest_ranger` works).
    - Encountered `conflicted` package error during comparison step.
    - Fixed script and re-running (`compare_ranger_rf_v2.log`).
- **Status**: Running. Expected completion in ~6 minutes.

## Next Steps
- Monitor `run_failed_clusters_v2.R` progress.
- Analyze `compare_ranger_rf_v2.log` results.
- Once verified, push changes to remote (if applicable) or finalize session.
