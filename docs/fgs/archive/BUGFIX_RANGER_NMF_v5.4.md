# FGS v5.4: random_forest_ranger and nmf_loadings Bug Fixes

## Date
2025-11-21

## Issues Fixed

### 1. `random_forest_ranger` - Name Mapping Issue

**Problem**: 
- When using `control_vars`, `ranger` internally applies `make.names()` to column names when creating the model via formula interface
- This caused importance names to not match original gene names, leading to "subscript out of bounds" errors

**Solution**:
- Added explicit name mapping before and after `ranger` model training
- Store original column names before creating `ranger_data`
- Create mapping: `make.names(orig) -> orig`
- After extracting importance, map names back to original if they don't match
- Only keep valid mappings to avoid NA issues

**Code Changes** (`myR/R/signature.R`, lines 1287-1329):
```r
# Store original names and create mapping for restoration
orig_colnames <- colnames(X_ranger)
# Mapping: make.names(orig) -> orig
name_mapping <- setNames(orig_colnames, make.names(orig_colnames))

# ... ranger training ...

weights_magnitude_all <- ranger::importance(rf_model)

# Map names back to original before filtering/sorting
if (!all(names(weights_magnitude_all) %in% orig_colnames)) {
  mapped_names <- name_mapping[names(weights_magnitude_all)]
  valid_idx <- !is.na(mapped_names)
  weights_magnitude_all <- weights_magnitude_all[valid_idx]
  names(weights_magnitude_all) <- mapped_names[valid_idx]
}
```

**Test Results**:
- Tested with 2,500 cells, 51,795 genes
- Successfully identified 50 gene signatures
- No errors during importance extraction or gene subsetting

---

### 2. `nmf_loadings` - do.call() Function Reference Issue

**Problem**:
- `do.call(NMF::nmf, nmf_args)` resulted in error: `'what' must be a function or character string`
- This occurs when `NMF::nmf` is not properly resolved as a function object in certain R environments

**Solution**:
- Store function reference explicitly before calling `do.call()`
- Use direct function reference instead of namespace-qualified call

**Code Changes** (`myR/R/signature.R`, lines 1580-1582):
```r
# Use direct function reference instead of string
nmf_func <- NMF::nmf
nmf_res <- do.call(nmf_func, nmf_args)
```

**Test Results**:
- Tested with 153 cells (Naive B-cells cluster)
- Successfully completed NMF factorization in 11.6 seconds
- Identified 30 gene signatures
- No errors during `do.call()` execution

---

## Verification

Both methods were tested on real data:
- **Dataset**: IS6_sex_added_251110_ds2500.qs (2,500 cells, 51,795 genes)
- **Target variable**: g3 (2 groups)
- **Success**: Both methods completed successfully and returned valid gene signatures

## Related Issues
- Closes issue with `random_forest_ranger` failing when using control variables
- Closes issue with `nmf_loadings` immediate failure on do.call()

## Breaking Changes
None - these are bug fixes that restore intended functionality.

## Notes
- `ranger` warning about formula interface for high-dimensional data is expected and can be ignored
- For large datasets, consider using the x/y interface directly instead of formula interface to avoid potential "protection stack overflow" errors
