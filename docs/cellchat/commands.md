# CellChat Analysis Guide & Commands

## ✅ 1. Proper Aggregation Method (Recommended)
> **Like Pseudobulk DEG Analysis**: Rigorous and statistically sound.

### Understanding the Proper Method
**Before (WRONG - Pooled Naively):**
```
All cells from 15 patients → One big CellChat analysis
Problem: Pseudo-replication, inflated interactions, p-value inflation.
```

**After (CORRECT - Sample-wise & Aggregate):**
```
Step 1: Analyze each patient independently (split by 'hos_no')
  Patient1 (g3=2) → CellChat1
  ...
  Patient15 (g3=2) → CellChat15

Step 2: Merge by condition (aggregate by 'g3')
  g3="2": mergeCellChat([CellChat1...15]) → "Stroke" merged object
  g3="1": mergeCellChat([CellChat16...23]) → "Control" merged object

Step 3: Compare
  Stroke vs Control
```

### Usage Examples

#### A. Full Proper Analysis (All DBs, All Patients)
```bash
# Analyze all 3 DB types, with proper patient-level analysis
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling,ECM-Receptor,Cell-Cell Contact" \
  -m 20 \
  -p 0.05 \
  -o /data/user3/sobj/cci/cellchat/proper_full \
  -c 8
```
*   `-p 0.05`: Probability threshold for identifying expressed genes.
*   `-c 8`: Number of cores for parallel processing.

#### B. Secreted Signaling Only (Most Common)
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling" \
  -m 20 \
  -o /data/user3/sobj/cci/cellchat/proper_secreted \
  -c 8
```

---

## 2. Alternative: Pooled Analysis (High Sensitivity)
> **Use Case**: When individual samples have too few cells (<500) to detect interactions, or for initial discovery.
> **Note**: This treats all cells in a condition as one "Super Sample".

```bash
# Pooled Analysis (Split by condition 'g3', Filter specific groups "2,1")
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s g3 \
  --subset_split "2,1" \
  -d "Secreted Signaling" \
  -m 50 \
  -o /data/user3/sobj/cci/cellchat/pooled_secreted_strict \
  -c 8
```
*Tip: Use higher `-m` (min.cells, e.g., 50) and strict probability thresholds to control false positives.*

---
## 3. Comparison
```bash
Rscript scripts/cellchat/run_cellchat_comparison_analysis.R \
  --file1 /data/user3/sobj/cci/cellchat/permissive_run2/merged/X1/cellchat_merged.qs \
  --file2 /data/user3/sobj/cci/cellchat/permissive_run2/merged/X2/cellchat_merged.qs \
  --name1 "Control" \
  --name2 "Stroke" \
  -o /data/user3/sobj/cci/cellchat/permissive_run2/comparison \
  --renv /home/user3/GJC_KDW_250721/renv/
```

## 3. Plotting & Visualization

### Individual Plots with Filtering (New!)
Use `-t` (Top Links) and `-k` (Top Interactions per Link) to clean up hairy plots.

```bash
# -t 50: Keep only Top 50 strongest Cell-Cell Links
# -k 20: Keep only Top 20 interactions per Link (cleaning Bubble plots)
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/plot_cellchat_individual.R \
  -i /data/user3/sobj/cci/cellchat/pooled_secreted_strict/samples/X2/cellchat.qs \
  -o /data/user3/sobj/cci/cellchat/pooled_secreted_strict/plots_refined \
  -t 50 \
  -k 20
```

### Comparison Plots
Use `plot_cellchat_comparison.R` to compare Control vs Stroke.
```bash
Rscript scripts/cellchat/plot_cellchat_comparison.R \
  --file1 /data/user3/sobj/cci/cellchat/permissive_run2/merged/X1/cellchat_merged.qs \
  --file2 /data/user3/sobj/cci/cellchat/permissive_run2/merged/X2/cellchat_merged.qs \
  --name1 "Control" \
  --name2 "Stroke" \
  -o /data/user3/sobj/cci/cellchat/permissive_run2/merge_comparison \
  --renv /home/user3/GJC_KDW_250721/renv/
```

---

## 4. Inspection Tool
Inspect the internal structure (`net`, `netP`, `LRsig`) of a CellChat object.

```bash
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/inspect_cellchat.R \
  -i /data/user3/sobj/cci/cellchat/pooled_secreted_strict/samples/X2/cellchat.qs \
  -n 20
```

---

## Appendix: Why Proper Method is Better

| Aspect | Old (Pooled) Method | Proper (Sample-wise) Method |
|--------|---------------------|-----------------------------|
| **Statistical rigor** | Low (Pseudo-replication) | High (Respects bio-replicates) |
| **Comparability** | Hard to compare P-values | Can use diffInteractions properly |
| **Sensitivity** | Very High (can be noisy) | Moderate (depends on sample size) |
| **Analogy** | Bulk RNA Normalization | Pseudobulk DEG Analysis |

## Appendix: NicheNet File Organization
Files related to NicheNet analysis have been renamed for clarity:
- `cci_nichenet_run.R` (formerly `cci_run_analysis.R`)
- `cci_nichenet_utils.R` (formerly `utils_cci.R`)
- `cci_nichenet_wrapper.R` (formerly `CCI.R`)