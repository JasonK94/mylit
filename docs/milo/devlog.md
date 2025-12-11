# Milo Pipeline Development Log

## 2025-12-11: Plot Caching Optimization and Density Visualization

### 1. Issue: Massive Plot File Size & Serialization Errors

#### Problem
- `run_milo_pipeline` execution resulted in extremely large `_05_plots.qs` files (30GB~50GB).
- Loading the package via `start.R` caused `renv` staging errors (`R6`, `scales` not found) and `View()` crashes in RStudio.
- **Root Cause**: The `ggplot` objects created within the pipeline captured the entire parent environment in their closure. This included the large `milo` (SingleCellExperiment) and `Seurat` objects, causing them to be serialized along with the plot.

#### Failed Attempts
- `qs::qsave` instead of `saveRDS`: Improved speed but did not solve the size issue.
- Simple environment cleaning: Removing objects from `plot_env` was insufficient because `patchwork` objects and nested layers maintained references.

#### Solution
- **Pipeline Refactoring**: Split plotting into `.milo_prepare_plot_data` (creates lightweight `_04_plot_inputs.qs`) and `.milo_render_plots`.
- **Environment Detachment**: Implemented a rigorous `.milo_clean_plot_environments` function that recursively replaces the environment of `ggplot` and `patchwork` objects with `emptyenv()`.
- **Output Format Change**: Switched `_05_plot_bundle` from a `.qs` object to a directory of direct PNG/PDF exports. This completely avoids serializing the R environment.

### 2. Feature: Weighted Kernel Density Plot

#### Initiative
- Users needed a way to visualize `logFC` or `SpatialFDR` distribution on the UMAP embedding, preserving the sign of change (up/down regulation).
- Existing point-based plots were hard to interpret in dense regions.

#### Implementation (`plots_density.R`)
- Created `plot_milo_density` function using `MASS::kde2d` (modified for weighting) logic.
- **Mode `density`**: Standard weighted density (absolute values).
- **Mode `average`**: Spatially smoothed weighted average (`sum(w*k) / sum(k)`), allowing visualization of positive/negative regions (e.g. logFC) using a diverging color scale.

#### Challenges & Fixes
1. **Coordinate Mapping Bug**:
   - Initial implementation of `geom_raster` dataframe had `x` and `y` repetition order swapped relative to `z` matrix (R uses column-major order).
   - *Fix*: Changed `expand.grid` logic to match `as.vector(matrix)`.

2. **Over-smoothing**:
   - `MASS::bandwidth.nrd` returned bandwidths too large for the UMAP structure (~20% of range), causing the density to look like a single blob in the center.
   - *Fix*: Applied a scaling factor (0.4) to the calculated bandwidth to preserve local structure.

3. **Masking Threshold**:
   - In `average` mode, low-density regions were masked using a 1% threshold relative to max density. Due to sharp peaks, this masked >99% of the plot.
   - *Fix*: Lowered threshold to 0.1% (`1e-3`) to reveal sparse regions.

4. **Stats Filter Error**:
   - `stats::filter` failed when the generated kernel size exceeded the grid size (due to large bandwidth or small grid).
   - *Fix*: Capped kernel radius at `floor((n-1)/2)`.

### 3. Summary of Cache Structure Changes
- **Old**: `_04_plots.rds` (Giant serialized object)
- **New**:
  - `_04_plot_inputs.qs`: Data only (Milo, DA results, metrics) - ~1GB
  - `_05_plot_bundle_*.png/pdf`: Final visualization images

