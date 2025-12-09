# CellChat Analysis - Improved Commands

## Current Status
- X2 (g3="2") completed: 624 interactions still seems high
- X1 (g3="1") incomplete: stopped during computeCommunProb

## Recommended: Use stricter filtering

### Option 1: Higher min.cells threshold
```bash
# Analyze g3="2" with stricter filtering
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -b g3 \
  --subset_split "2" \
  -m 50 \
  -o /data/user3/sobj/cci/cellchat/run4_g3split_anno3_strict

# Analyze both with strict filtering
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -b g3 \
  --subset_split "2,1" \
  -m 50 \
  -o /data/user3/sobj/cci/cellchat/run4_g3split_anno3_strict
```

### Option 2: Per-patient analysis (most conservative)
```bash
# This will create one analysis per patient - much cleaner signals
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -b hos_no \
  -m 30 \
  -o /data/user3/sobj/cci/cellchat/run4_per_patient
```

### Option 3: Subset to specific DB
Add this to the CLI later - use only "Secreted Signaling" interactions

## Analysis
Current X2 results:
- 24,866 cells (g3="2" only - correctly split!)
- 25 cell types
- 624 interactions (still quite high)
- 45 pathways

The split worked correctly, but you're right that the signal is still strong. This is because:
1. Still analyzing 24,866 cells together (multiple patients combined within g3="2")
2. All interaction types included (Secreted, ECM, Cell-Cell Contact)

Better approach: Split by individual patient (hos_no) for cleanest signal.