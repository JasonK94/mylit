# Function Analysis Error Log

This document records issues encountered during the automated analysis of R scripts in the `myR/R/` directory.

## 1. Analysis Process Adjustments

### 1.1. Handling of Duplicate Functions

**Issue:** The analysis process initially documented only the last occurrence of functions that were defined multiple times within the same file or across different files.

**Resolution:** The process has been updated to document **all** occurrences of every function, noting the source file for each. This is crucial for the refactoring process to have a complete overview of all definitions. A specific instance of this was the `sort_samples` function, which is defined twice in `myR/R/utils.R`.

## 2. File-Specific Issues

### 2.1. `myR/R/deprecated/test_with_parse_error.R`

- **Problem:** **Parsing Error**.
- **Details:** This file contains a mix of R script code and R Markdown chunks (` ```{r} ... ``` `). This non-standard format prevents successful parsing as a pure R script, so it has been excluded from the function analysis.

### 2.2. `myR/R/markers.R`

- **Problem:** **File Write Error (Previous attempts)**.
- **Details:** Previous attempts to write the analysis results for this file to `myR/function_analysis.md` failed. The analysis will be re-attempted.

### 2.3. `myR/R/pseudobulk_deg.R`

- **Problem:** **File Read/Grep/Write Error (Previous attempts)**.
- **Details:** This file is exceptionally large and has a complex, non-standard structure, causing failures in reading, searching (`grep`), and writing its analysis. The analysis will be re-attempted methodically.

### 2.4. Duplicate Filenames (`demulti_utils.R`)

- **Files:**
    1. `myR/R/demulti_utils.R`
    2. `myR/R/utilities/demulti_utils.R`
- **Details:** Two files with the same name exist. They contain functions with identical names (e.g., `get_best_two`, `get_barcode_mapping`) but **different implementations and function signatures**. Both sets of functions have been documented under their respective file paths in the main analysis document.
