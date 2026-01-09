# Development Summary - Phase 2 Refactoring (Final)

## Objective
Refactor DEG analysis pipeline to remove `muscat` dependency (standalone mode), improve robustness (rank deficiency handling), and restore consolidated output structure.

## Key Changes

### 1. Robust Design Matrix Handling (`handle_design_v1`)
- **Problem**: Analysis failed due to rank deficient design matrices (e.g., collinearity between Batch/Set and Group) and numeric covariate scaling issues.
- **Solution**: Implemented `handle_design_v1` utility in `deg_methods_base.R`.
  - Automatically scales numeric covariates.
  - Detects rank deficiency and greedily drops problematic terms (e.g., `set`) to ensure full rank.
  - Returns cleaned metadata and formula.
- **Integration**: Applied to `edgeR` (LRT, QLF), `DESeq2` (Wald, LRT), `limma` (voom, trend), and `dream`.

### 2. Standalone Methods Refactoring
- All methods now use `create_pseudobulk_v1` natively, removing `muscat` wrapper dependency for core logic.
- **runNEBULA**: Updated to support `max_genes = NULL` (no limit) and use `doFuture` for parallelism.
- **runDREAM**: Refactored to use `handle_design_v1` for Fixed Effects while preserving Random Effects (`(1|batch)`), and to generally track parameters.

### 3. Metadata & Traceability
- Added `run_info` (captured via `match.call()`) and `formula` attributes to results objects of all methods.
- Ensures exact parameters used for each run are preserved.

### 4. Consolidated Output Restoration
- Updated `run_production_is6_refined.R` to perform aggregation steps explicitly:
  - `standardize_deg_results`
  - `build_deg_matrices`
  - `compute_agreement_scores`
  - `compute_consensus_scores`
  - `generate_consensus_deg_list`
- Resulting `.qs` files now contain the full suite of `consensus_scores` and `consensus_deg_list` as requested.

### 5. Verification
- **Test Run**: Executed `run_production_is6_refined.R` on IS6 dataset (Global level).
- **Result**: `edgeR`, `limma` methods successfully identified and resolved rank deficiency ("Dropping candidate term: set"), converting failures into successful completions.

## Next Steps (Phase 3)
- **Method Comparison**: Execute the comparison plan (`DEV_CONS_PLAN.md`) using the generated results.
- **Optimization**: Monitor `dream` performance on larger subsets.
