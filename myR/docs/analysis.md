# Differential Expression Function Guide

## Overview
- All DE helpers that mention **MAST**, **MUSCAT**, or **NEBULA** now live in `myR/R/analysis.R`.
- `runNEBULA` (formula/pseudobulk switches) and `runMUSCAT` (improved NA handling) provide the unified APIs.
- Older snippets from `test_analysis.R`, `test.R`, `test_cursor.R`, and `test_working.R` were suffixed with `_legacy` and archived in `myR/R/test_to_delete.R`.
- This consolidation ensures a single source of truth for NEBULA, MUSCAT, limma/edgeR/DESeq2 pipelines while keeping the older experiments available for reference.

## Function Inventory
| Function | Purpose | Key Inputs | Previous Location |
| --- | --- | --- | --- |
| `runMAST_v1` | Single-cell hurdle model via MAST | `formula`, `lrt_variable`, `min_cells_expr` | `myR/R/test_analysis.R` (pre 7fe20d3) |
| `runMUSCAT` | Unified MUSCAT pseudobulk helper (edgeR/DESeq2/limma backends) | `cluster_id`, `sample_id`, `group_id`, `contrast`, `method`, `remove_na_groups` | `myR/R/analysis.R` |
| `runNEBULA` | Unified NEBULA helper (classic, formula, pseudobulk modes) | `fixed_effects`, `covar_effects`, `formula`, `pseudobulk_args`, `patient_col` | `myR/R/analysis.R` |
| `run_deg_consensus` | Multi-method orchestrator (muscat/NEBULA + limma/edgeR/DESeq2) | `methods`, `contrast`, `cluster_id`, `sample_id`, `group_id`, `batch_id` | `myR/R/deg_consensus/run_deg_consensus.R` |
| `run_deg_consensus_with_plots` | Full pipeline (consensus + plots + .qs snapshots) | `methods`, `significance_mode`, `output_dir`, `make_plots` | `myR/R/deg_consensus/deg_consensus_pipeline.R` |
| `_legacy` wrappers | Deprecated entry points (`runMUSCAT (legacy mode)`, `runMUSCAT`, `runNEBULA2_*`) kept for backwards compatibility | Same as before (emit `.Deprecated` warnings)` | `myR/R/test_to_delete.R` |

Legacy wrappers that previously lived in `test.R` or `test_cursor.R` are now preserved as `run<NAME>_legacy()` functions inside `myR/R/test_to_delete.R` for historical debugging.

## Complex Formula Examples
The engines accept full lme4-style formulas. Tested templates:

- **runMAST_v1** (hurdle model, no random effects inside the formula — random intercepts handled via `random = ~1|patient` in the `MAST::zlm` call):
  ```
  ~ g3 * sex * anno3.scvi + GEM + percent.mt + nFeature_RNA
  ```
  Interaction terms are expanded before fitting; include nuisance covariates directly.

- **runMUSCAT** (design matrix built internally):
  ```
  contrast = "groupIS - groupSAH + batch2 - batch1"
  ```
  Nested or crossed factors should be encoded via explicit contrasts. `runMUSCAT` keeps the same syntax while `remove_na_groups` toggles the legacy behaviour.

- **runNEBULA(formula = …)** (full support for multi-way interactions and nested random effects):
  ```
  ~ g3 * sex * anno3.scvi + GEM + percent.mt +
    g3:anno3.scvi + sex:anno3.scvi +
    (1 | GEM/hos_no) + (1 | patient_batch)
  ```
  - Nested term `(1 | GEM/hos_no)` is converted to fixed `GEM` plus random `hos_no`.
  - Multiple random-effect terms are parsed; the first one (or explicit `patient_col`) becomes the NEBULA `id`.
  - Interactions are kept verbatim in the design matrix, and complete-separation warnings are issued via contingency checks.

## Migration Notes
1. **Unified in `analysis.R`:** `runMAST_v1`, `runMUSCAT`, and `runNEBULA` (with formula/pseudobulk switches).
2. **New consensus module:** Everything under `myR/R/deg_consensus/` originates from the dedicated `_wt/deg-consensus` worktree and is now shipped inside the package. The helpers (`run_deg_consensus*`, method-specific wrappers, matrices/plots) are callable directly after `library(myR)`.
3. **Wrappers:** Deprecated names (`runMUSCAT (legacy mode)`, `runMUSCAT`, `runNEBULA2_*`) simply forward to the new helpers and emit warnings.
4. **Archived:** Older experiments remain in `myR/R/test_to_delete.R` as `_legacy` helpers.
5. **Deprecated entry points:** `myR/R/test_analysis.R`, `test.R`, `test_cursor.R`, and `test_working.R` no longer define DE helpers. They simply remain as scaffolding or comments.

## Testing & Debugging
- Quick regression: `Rscript myR/tests/test_formula1_light.R` (see `_wt/analysis/scripts/analysis/test_formula1_light.R`) ensures `runNEBULA(formula = ...)` still handles nested random effects.
- Full analysis: `Rscript scripts/analysis/run_formula1_analysis.R` exercises the unified `runNEBULA` pipeline against `IS6_sex_added_251110.qs`.
- MUSCAT smoke test: `Rscript test_run_muscat2.R --config <...>` now calls `runMUSCAT` for pseudobulk + edgeR.
- Consensus pipeline: `run_deg_consensus_with_plots()` saves raw/standardised results plus QA plots (PCA/UMAP/heatmap) alongside `.qs` copies in `/data/user3/sobj`.

For new datasets, prefer the formula interface so you can include interactions and multi-level random intercepts explicitly.

## Formula Guard & Debugging

`runNEBULA(formula = …)` now carries a built-in guard against complete separation and singular design matrices:

- `separation_min_cells` (default: 20) prunes factor levels that do not have enough cells before the design matrix is constructed. The pruning log is returned under `result$separation_guard$removal_log`.
- `separation_group_var` + `separation_min_cells_per_group` remove genes that fail to appear in every group (e.g., require ≥5 cells with counts > 0 for each `g3` level). Counts are inspected before `nebula::group_cell()`, and the number of dropped genes is stored in `result$separation_guard$gene_balance`.
- `simplify_interactions = TRUE` drops interaction terms that involve variables flagged by the separation detector. The final design formula is stored in `result$design_formula`, and removed interactions are listed in `result$separation_guard$simplified_interactions`.
- Diagnostics (pairwise contingency tables, QR rank, retry attempts) are echoed to the console and saved in `result$separation_guard$diagnostics`.

Script defaults:

- `run_formula1_analysis.R`: `separation_min_cells = 25`, `separation_group_var = "g3"`, and `separation_min_cells_per_group = 5`.
- `test_formula1_light.R`: `separation_min_cells = 15`, `separation_group_var = "g3"`, and `separation_min_cells_per_group = 3`.

When separation persists, inspect the printed contingency tables and either (1) drop the offending interaction or (2) remove/merge the sparse metadata levels upstream.

## DEG Consensus Pipeline

`myR/R/deg_consensus/deg_consensus_pipeline.R` exposes `run_deg_consensus_with_plots()` which orchestrates:

1. `run_deg_consensus()` across MUSCAT + limma/edgeR/DESeq2 helpers.
2. Standardisation (`standardize_deg_results()`), matrix construction (`build_deg_matrices()`), and consensus scoring.
3. Optional PCA/UMAP/volcano/heatmap outputs alongside `.qs` snapshots inside `/data/user3/sobj`.

The helper requires only a Seurat object plus a `contrast`. All intermediate files honour the `output_dir`/`prefix` arguments and are exported with conflict-safe filenames.

## Testing & Debugging

- Quick regression: `Rscript myR/tests/test_formula1_light.R` (see `_wt/analysis/scripts/analysis/test_formula1_light.R`) ensures `runNEBULA(formula = …)` still handles nested random effects **with** the separation guard enabled.
- Full analysis: `Rscript scripts/analysis/run_formula1_analysis.R` exercises the unified `runNEBULA` pipeline against `IS6_sex_added_251110.qs`. Use `scripts/analysis/check_formula1_status.sh` to tail `/tmp/formula1_full_analysis.log` and monitor running jobs + output `.qs` artefacts.
- Interactive triage: `Rscript scripts/analysis/run_formula1_interactive.R` is wired for stepwise debugging (pause between steps, inspect metadata tables, tweak formulas).
- Minimal example: `Rscript scripts/analysis/run_formula1_simple_example.R` illustrates how to call the formula interface on a tiny subset without touching `/data/user3/sobj`.
- MUSCAT smoke test: `Rscript test_run_muscat2.R --config <...>` now calls `runMUSCAT` for pseudobulk + edgeR.
- Consensus pipeline: `run_deg_consensus_with_plots()` saves raw/standardised results plus QA plots (PCA/UMAP/heatmap) alongside `.qs` copies in `/data/user3/sobj`.

All scripts must run from `/home/user3/GJC_KDW_250721` so that `start.R` wires up the environment before `devtools::load_all(myR)`.
