# Differential Expression Function Guide

## Overview
- All DE helpers that mention **MAST**, **MUSCAT**, or **NEBULA** now live in `myR/R/analysis.R`.
- Older snippets from `test_analysis.R`, `test.R`, `test_cursor.R`, and `test_working.R` were suffixed with `_legacy` and archived in `myR/R/test_to_delete.R`.
- This consolidation ensures a single source of truth for NEBULA, MUSCAT, limma/edgeR/DESeq2 pipelines while keeping the older experiments available for reference.

## Function Inventory
| Function | Purpose | Key Inputs | Previous Location |
| --- | --- | --- | --- |
| `runMAST_v1` | Single-cell hurdle model via MAST | `formula`, `lrt_variable`, `min_cells_expr` | `myR/R/test_analysis.R` (pre 7fe20d3) |
| `runMUSCAT_v5` | Pseudobulk DE via muscat (edgeR/DESeq2/limma) | `cluster_id`, `sample_id`, `group_id`, `contrast`, `method` | `myR/R/test_analysis.R` |
| `runMUSCAT2_v1` | Hardened muscat pipeline with NA scrubbing and stricter filtering | Adds `remove_na_groups`, `pb_min_cells`, `filter_genes` | `myR/R/test_analysis.R` |
| `runNEBULA_v1` | Original NEBULA GLMM helper | `fixed_effects`, `covar_effects`, `patient_col` | `myR/R/test_analysis.R` (pre 7fe20d3) |
| `runNEBULA2_v1` | Improved NEBULA with NA checks & separation diagnostics | `fixed_effects`, `covar_effects`, `patient_col`, `min_count` | `myR/R/test_analysis.R` |
| `runNEBULA2_v1_with_pseudobulk` | Two-step muscat + NEBULA bridge | `cluster_id`, `sample_id`, `group_id`, `nebula_fixed` | `myR/R/test_analysis.R` |
| `runNEBULA2_v1_with_formula` | Formula-driven NEBULA interface with nested random-effect parser | `formula`, `patient_col`, `offset`, `min_count` | `myR/R/test_analysis.R` |

Legacy wrappers that previously lived in `test.R` or `test_cursor.R` are now preserved as `run<NAME>_legacy()` functions inside `myR/R/test_to_delete.R` for historical debugging.

## Complex Formula Examples
The engines accept full lme4-style formulas. Tested templates:

- **runMAST_v1** (hurdle model, no random effects inside the formula — random intercepts handled via `random = ~1|patient` in the `MAST::zlm` call):
  ```
  ~ g3 * sex * anno3.scvi + GEM + percent.mt + nFeature_RNA
  ```
  Interaction terms are expanded before fitting; include nuisance covariates directly.

- **runMUSCAT_v5 / runMUSCAT2_v1** (design matrix built internally):
  ```
  contrast = "groupIS - groupSAH + batch2 - batch1"
  ```
  Nested or crossed factors should be encoded via explicit contrasts. `runMUSCAT2_v1` keeps the same syntax but adds NA cleaning before aggregation.

- **runNEBULA2_v1_with_formula** (full support for multi-way interactions and nested random effects):
  ```
  ~ g3 * sex * anno3.scvi + GEM + percent.mt +
    g3:anno3.scvi + sex:anno3.scvi +
    (1 | GEM/hos_no) + (1 | patient_batch)
  ```
  - Nested term `(1 | GEM/hos_no)` is converted to fixed `GEM` plus random `hos_no`.
  - Multiple random-effect terms are parsed; the first one (or explicit `patient_col`) becomes the NEBULA `id`.
  - Interactions are kept verbatim in the design matrix, and complete-separation warnings are issued via contingency checks.

## Migration Notes
1. **Moved to `analysis.R`:** `runMAST_v1`, `runMUSCAT_v5`, `runMUSCAT2_v1`, `runNEBULA_v1`, `runNEBULA2_v1`, `runNEBULA2_v1_with_pseudobulk`, `runNEBULA2_v1_with_formula`, plus compatibility wrappers (`runMAST`, `runMUSCAT`, `runNEBULA`).
2. **Archived:** Older wrappers now live in `myR/R/test_to_delete.R` as `_legacy` variants.
3. **Deprecated entry points:** `myR/R/test_analysis.R`, `test.R`, `test_cursor.R`, and `test_working.R` no longer define DE helpers. They simply remain as scaffolding or comments.

## Testing & Debugging
- Quick regression: `Rscript myR/tests/test_formula1_light.R` (see `_wt/analysis/test_formula1_light.R` for reference) ensures `runNEBULA2_v1_with_formula` still handles nested random effects.
- Full analysis: `Rscript run_formula1_analysis.R` (paths in `_wt/analysis/`) exercises the entire NEBULA2 pipeline against `IS6_sex_added_251110.qs`.
- MUSCAT smoke test: `Rscript test_run_muscat2.R --config <...>` runs through pseudobulk + edgeR.

For new datasets, prefer the formula interface so you can include interactions and multi-level random intercepts explicitly.
