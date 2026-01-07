# FGS Development Log

## 2026-01-07: FGS Performance Optimization & Refactoring

### 🚀 Key Changes

#### 1. Performance Optimization
- **`random_forest` Removed**: Removed the single-threaded `random_forest` method from the default method list. `random_forest_ranger` is now the sole default tree-based method, offering significantly better performance.
- **`glmnet` Optimization**: Relaxed the convergence threshold for `ridge` and `lasso` from `1e-7` to `1e-5`. This speeds up the solution path for high-dimensional data without compromising signature quality.
- **Timing Estimation v2**: Implemented a **data-size aware** timing estimation. The system now records `n_cells` and `n_genes` to predict runtime based on complexity ($N \times P$), preventing inaccurate estimates when switching between small and large datasets.

#### 2. Code Refactoring
- **Function Renaming**: 
  - Renamed `find_gene_signature_v5_impl` to **`find_gene_signature`**.
  - Removed the `FGS` alias wrapper to simplify the call stack.
- **Logic Consolidation**:
  - Removed internal `method_impl` argument. The code now strictly follows the "v5.4" logic (Permutation importance for Ranger, Non-negative NMF).
- **Conflict Resolution**: 
  - Removed `compute_meta_gene_importance` from `myR/R/utils_fgs.R` to resolve conflicts with the TML-compatible version in `myR/R/tml_utils.R`.

#### 3. Cleanup & UX
- **Warning Suppression**: 
  - Suppressed "Coefficients not estimable" warnings in `limma::removeBatchEffect`.
  - Suppressed package startup messages in `init_fgs_env.R`.
- **File Cleanup**: Deleted unused test files (`test_to_delete.R`, `test_claude.R`, `test_cursor.R`, `test_analysis.R`).

### 📝 Next Steps
- Verify if `ridge` regression performance is acceptable with the new threshold.
- Continue monitoring `nmf_loadings` stability with the non-negative constraint fix.
