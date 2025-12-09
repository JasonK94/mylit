# CellChat Analysis - Proper Aggregation Method

## ✅ NEW: Correct Implementation (Like Pseudobulk DEG)

### Understanding the Proper Method

**Before (WRONG - what we were doing):**
```
All cells from 15 patients → One big CellChat analysis
Problem: Pseudo-replication, inflated interactions
```

**After (CORRECT - now implemented):**
```
Step 1: Analyze each patient independently
  Patient1 (g3=2) → CellChat1
  Patient2 (g3=2) → CellChat2
  ...
  Patient15 (g3=2) → CellChat15

Step 2: Merge by condition
  g3="2": mergeCellChat([CellChat1...15]) → "Stroke" merged object
  g3="1": mergeCellChat([CellChat16...23]) → "Control" merged object

Step 3: Compare
  compareInteractions(Stroke, Control)
```

This is analogous to **pseudobulk DEG analysis** - more rigorous and statistically sound!

## Usage Examples

### Example 1: Full Proper Analysis (15 stroke + 8 control patients)
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
  -c 16
```

**What this does:**
- `-s hos_no`: Each of 23 patients analyzed independently
- `-a g3`: Merge 15 stroke patients, 8 control patients separately  
- `-d ...`: Run 3 separate analyses (Secreted, ECM, Cell-Cell)
- Result: 3 directories, each with:
  - `samples/`: Individual patient CellChat objects (23 total)
  - `merged/X2/`: Merged stroke CellChat (15 patients)
  - `merged/X1/`: Merged control CellChat (8 patients)

### Example 2: Stroke Patients Only, All 3 DB Types
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2" \
  -d "Secreted Signaling,ECM-Receptor,Cell-Cell Contact" \
  -m 20 \
  -o /data/user3/sobj/cci/cellchat/proper_stroke
```

### Example 3: Just Secreted Signaling (Most Common in Papers)
```bash
# This is what most papers do
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling" \
  -m 20 \
  -o /data/user3/sobj/cci/cellchat/proper_secreted
```

## Parameters Explained

- **`-s/--split_by`** (REQUIRED for proper method): Sample ID (e.g., `hos_no`)
  - Each sample gets independent CellChat analysis
  - Like pseudobulk: respect sample structure
  
- **`-a/--aggregate_by`** (Recommended): Condition (e.g., `g3`)
  - Merges split samples by condition using `mergeCellChat()`
  - Enables condition comparison
  
- **`-d/--db_use`**: Comma-separated DB types
  - If multiple specified, runs each as separate analysis
  - Most papers focus on "Secreted Signaling" only
  
- **`-p/--prob_threshold`**: Probability cutoff (default: 0.05)
  - Lower = stricter filtering
  - Usually don't need to change
  
- **`-m/--min_cells`**: Min cells per cell type (default: 10)
  - Cell type inclusion threshold
  - Not the main interaction filter

## Output Structure

```
output_dir/
├── Secreted.Signaling/
│   ├── samples/
│   │   ├── Patient1/cellchat.qs
│   │   ├── Patient2/cellchat.qs
│   │   └── ...
│   └── merged/
│       ├── X2/cellchat_merged.qs  (stroke, 15 patients)
│       └── X1/cellchat_merged.qs  (control, 8 patients)
├── ECM.Receptor/
│   └── (same structure)
└── Cell.Cell.Contact/
    └── (same structure)
```

## Why This is Better

| Aspect | Old Method | New Method |
|--------|------------|------------|
| **Statistical rigor** | Pools all cells | Respects sample structure |
| **Replication** | Pseudo-replication | True biological replicates |
| **Comparability** | Can't compare conditions | Can compare via mergeCellChat |
| **Publications** | Rarely accepted | Standard method |
| **Analogy** | Bulk RNA normalization | Pseudobulk DEG |

## Expected Results

With proper method + Secreted Signaling only:
- **Per patient**: ~50-100 interactions
- **Merged (stroke)**: ~100-200 interactions  
- **Much cleaner** bubble plots
- **Statistically sound** for comparison