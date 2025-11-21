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

## 3. Next Steps
- Monitor `run_failed_clusters_v2.R` progress.
- Verify NMF and Ranger L1 methods on full dataset (Step 3 of Plan).
- Compare Ranger vs Random Forest results (Step 4 of Plan).
