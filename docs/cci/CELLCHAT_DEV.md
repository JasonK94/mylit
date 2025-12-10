# CellChat Development & Workflow Documentation

## 1. Overview
This document details the development history, workflow architecture, and troubleshooting logs for the CellChat integration in this repository. The goal is to provide a robust pipeline for cell-cell interaction analysis, balancing statistical rigor (pseudobulk) with sensitivity (pooled analysis), and ensuring reproducibility.

## 2. Architecture & Modules

### 2.1 Core Wrapper (`myR/R/cci_cellchat_wrapper.R`)
- **`run_cellchat_analysis`**: The main driver function.
  - Handles Seurat object conversion to CellChat object.
  - Implements **two aggregation strategies**:
    1. **Pooled (Condition-based)**: Merges cells by condition (e.g., `g3`) before analysis. High power but potential false positives.
    2. **Sample-wise (Pseudobulk)**: Analyzes each sample (`hos_no`) independently, then averages. Rigorous but low sensitivity for rare interactions.
  - Supports **Checkpointing** (`step1_preprocess`, `step2_verify`, etc.) to resume long-running jobs.
  - Dynamic **Parallel Processing**: Automatically selects backend (`multicore` vs `multisession`) and core count based on cell numbers to strictly avoid overhead on small samples.

### 2.2 CLI Interface (`scripts/cellchat/run_cellchat_cli.R`)
- Provides a command-line interface for the wrapper.
- Key Options:
  - `-i`: Input Seurat object (.rds/.qs).
  - `-s`: Split column (e.g., `g3` for pooled, `hos_no` for sample-wise).
  - `--subset_split`: Filters specific groups from the split column (e.g., `"2,1"`).
  - `-m`: `min.cells` threshold.
  - `-d`: Database use (`Secreted Signaling`, `ECM-Receptor`, etc.).
  - `--subset_aggregate`: (Deprecated/Internal) Used for aggregating sample-wise results.

### 2.3 Plotting Modules
- **`scripts/cellchat/plot_cellchat_individual.R`**: Plots standard visualizations (Circle, Bubble) for a single (or merged) CellChat object.
- **`scripts/cellchat/plot_cellchat_comparison.R`**: Performs differential analysis (Circle Diff, Info Flow) between two CellChat objects.
- **`myR/R/cci_cellchat_plotting.R`**: Shared plotting functions.

## 3. Workflow Strategies

### Strategy A: Pooled Analysis (Recommended for Discovery)
- **Concept**: Treat all cells from patients in Group A as one "Super Sample".
- **Pros**: High sensitivity. Detects interactions even if rare in individual patients.
- **Cons**: Can yield **too many interactions** (noisy). P-values may be inflated due to pseudoreplication.
- **Mitigation**: 
  - Strict `min.cells` (e.g., 50-100).
  - Filter results by interaction probability quantile (e.g., Top 5%).
  - Focus on "Secreted Signaling" DB only.

### Strategy B: Sample-wise Analysis (Strict/Pseudobulk)
- **Concept**: Analyze Patient 1, Patient 2... independently. Then merge.
- **Pros**: Statistically sound. Removes patient-specific noise.
- **Cons**: If sample size is small (<1000 cells/patient) or cell types are rare, **ZERO interactions** may be detected due to failure to pass `min.cells` or statistical thresholds in individual runs.

## 4. Development History & Troubleshooting

### 2025-12-10: Pipeline Refinement & Optimization
- **Issue 1: "Too Many Interactions" (Initial Pooled Run)**
  - **Symptom**: Circle plots were hairballs (dense connections).
  - **Cause**: Default CellChat settings on 40k cells generate excessive significant pairs.
  - **Fix**: 
    - Implemented `db.use` filtering.
    - Added `min.cells` control via CLI.
  
- **Issue 2: "Zero Interactions" (Sample-wise Run)**
  - **Symptom**: `proper_full` run resulted in empty plots.
  - **Cause**: Splitting by `hos_no` left some cell types with <20 cells. `min.cells=20` filtered them out completely. `computeCommunProb` failed to find significant interactions in small subsets.
  - **Lesson**: Sample-wise analysis requires large datasets per patient. For current data, **Pooled Analysis** is necessary.

- **Issue 3: Parallel Processing Overhead**
  - **Symptom**: Step 1 (Preprocessing) took hours for small samples.
  - **Cause**: `future` package overhead when `n_cells < 1000`.
  - **Fix**: Implemented dynamic core allocation. If `cells < 1000`, force `sequential` (1 core). Speedup >10x.

- **Issue 4: CLI Flag Error (`subset_split`)**
  - **Symptom**: `Error: long flag "subset_split" is invalid`.
  - **Cause**: Added logic to CLI but forgot to add `make_option` definition.
  - **Fix**: Added proper `make_option` entry.

## 5. Future Improvements
1. **Top N Filtering**: Implement a post-processing filter to show only Top N strongest interactions in plots.
2. **Differential Plot Automation**: Auto-generate "Control vs Disease" plots after analysis.
3. **HTML Reporting**: Wrap plots into an HTML summary.

---
**Author**: Antigravity (AI Agent)
**Date**: 2025-12-10
