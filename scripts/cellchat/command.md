# CellChat Analysis Commands

## Issue Analysis
The original analysis combined all patient samples together, resulting in:
- Over-aggregated signals (too strong)
- Bubble plots that are too dense to interpret

## Solution: Split by Sample/Condition

### Option 1: Split by condition (g3 status)
```bash
# Analyze only g3="2" patients
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -b g3 \
  --subset_split "2" \
  -m 20 \
  -o /data/user3/sobj/cci/cellchat/run_g3split_anno3

# Analyze both g3="2" and g3="1"
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -b g3 \
  --subset_split "2,1" \
  -m 20 \
  -o /data/user3/sobj/cci/cellchat/run_g3split_anno3
```

### Option 2: Split by patient (hos_no)
```bash
# Will create separate analysis for each patient
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -b hos_no \
  -m 20 \
  -o /data/user3/sobj/cci/cellchat/run_patient_split_anno3
```

### Option 3: Higher filtering threshold (no splitting)
```bash
# Use higher min.cells to reduce noise
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -m 50 \
  -o /data/user3/sobj/cci/cellchat/run_highfilt_anno3
```

## Parameters
- `-i`: Input Seurat object file
- `-g`: Cell type grouping column
- `-b`: Split by this metadata column (e.g., g3, hos_no)
- `--subset_split`: Comma-separated list of specific split groups to analyze
- `-m`: Minimum cells per cell type group (higher = stricter filtering, less noise)
- `-o`: Output directory
- `-c`: Number of cores (default: 16)