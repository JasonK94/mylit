# Deprecated Files

This directory contains files that have been deprecated during the refactoring process.

## Files in this directory:

### 1. `pseudobulk_deg_standalone.R`
**Status:** Deprecated - Functions consolidated  
**Reason:** All pseudobulk differential expression functions have been consolidated into:
- `R/analysis/differential_expression/differential_expression.R`

**Key functions that were moved:**
- `prepare_pseudobulk_edgeR()`
- `run_pseudobulk_deg()`

**What to do:** Use the new consolidated `differential_expression.R` module instead.

---

### 2. `test_with_parse_error.R`
**Status:** Deprecated - Parse error + Functions consolidated  
**Reason:** 
1. This file had a parse error that prevented `devtools::document()` from working
   - Error: "attempt to use zero-length variable name"
   - Parse error around line 35 with NSE (non-standard evaluation) code
   
2. The main LMM functions from this file have already been consolidated into:
   - `R/analysis/differential_expression/differential_expression.R`

**Key functions that were moved:**
- `create_analysis_config()` → Now in differential_expression.R
- `fit_lmm_single_gene()` → Now in differential_expression.R  
- `run_lmm_multiple_genes()` → Now in differential_expression.R
- `find_response_differential_genes()` → Now in differential_expression.R
- `find_drug_specific_genes()` → Now in differential_expression.R

**Additional functions in test.R (not yet consolidated):**
- `diagnosis_parity_legacy()` - Metadata parity checking
- `diagnosis_parity()` - Improved version
- `prepare_geomx_data()` - GeoMx data preparation
- `quick_screen_genes()` - Gene screening
- And several other GeoMx-specific helper functions

**What to do:**
- **For LMM analysis:** Use `differential_expression.R` - all main functions are there
- **For GeoMx-specific utilities:** These are also in `test_claude.R` which is still active
- **If you need the legacy functions:** You can manually source this file, but note it has a parse error that prevents package building

**Note:** The file was kept in deprecated because:
1. It contains some utility functions that might be useful for reference
2. The parse error makes it unsuitable for inclusion in the package
3. All critical functionality has been moved to properly working modules

---

## Package Building

With these files moved to deprecated/, the package now:
- ✅ Builds successfully with `devtools::document()`
- ✅ Loads successfully with `devtools::load_all()`  
- ✅ All 100+ functions properly documented
- ✅ No parse errors

## Migration Guide

If you were using functions from these deprecated files:

### From `pseudobulk_deg_standalone.R`:
```r
# OLD (deprecated)
source("myR/R/pseudobulk_deg.R")
results <- run_pseudobulk_deg(...)

# NEW (recommended)
library(myR)  # or devtools::load_all()
results <- run_pseudobulk_deg(...)  # Same function, now in differential_expression.R
```

### From `test_with_parse_error.R`:
```r
# OLD (had parse errors)
source("myR/R/test.R")
config <- create_analysis_config(...)
results <- run_lmm_multiple_genes(...)

# NEW (recommended)
library(myR)  # or devtools::load_all()
config <- create_analysis_config(...)  # Now in differential_expression.R
results <- run_lmm_multiple_genes(...)  # Now in differential_expression.R
```

---

**Date deprecated:** October 10, 2025  
**Refactoring branch:** refactoring-v0.2

