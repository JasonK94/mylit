# CellChat Analysis Scenarios & Commands

This document outlines different analysis strategies for CellChat, from strict per-sample analysis to pooled approaches, and when to use each.

---

## 1. Pooled Analysis (Recommended for Statistical Power)

**Use case**: When individual samples have too few cells to detect robust interactions. Pools cells by condition (e.g., all stroke patients together).
**Pros**: High statistical power, detects subtle signals.
**Cons**: Can mask patient heterogeneity (pseudo-replication issue).
**Control**: Use higher `min.cells` (e.g., 50-100) and restrict DB to "Secreted Signaling" to reduce noise.

```bash
# Run Pooled Analysis (Split by condition 'g3')
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s g3 \
  --subset_split "2,1" \
  -d "Secreted Signaling" \
  -m 50 \
  -o /data/user3/sobj/cci/cellchat/pooled_secreted_strict \
  -c 1
```

---

## 2. Proper Sample-wise Analysis (Strict & Rigorous)

**Use case**: Ideally the standard method. Analyzes each patient independently, then merges.
**Pros**: Statistically rigorous, respects biological replicates.
**Cons**: Requires sufficient cells per patient. If samples are small (<500 cells), interactions may be lost (Zero interactions).

```bash
# Run Sample-wise Analysis (Split by patient 'hos_no', Aggregate by 'g3')
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling" \
  -m 10 \
  -o /data/user3/sobj/cci/cellchat/proper_sample_wise \
  -c 1
```

---

## 3. High-Sensitivity Analysis (For Rare Populations)

**Use case**: When looking for rare cell types or subtle interactions.
**Settings**: Low `min.cells`, potentially all DB types.

```bash
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s g3 \
  -d "Secreted Signaling,ECM-Receptor" \
  -m 5 \
  -o /data/user3/sobj/cci/cellchat/sensitive_run \
  -c 1
```

---

## 4. Plotting & Comparison

After running analysis, generate comparison plots.

```bash
# Plot individual results
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/plot_cellchat_individual.R \
  -i /data/user3/sobj/cci/cellchat/pooled_secreted_strict/Secreted.Signaling \
  -o /data/user3/sobj/cci/cellchat/pooled_secreted_strict/plots

# (Future) Differential plotting script
```

# Plotting with Filtering
# -t 50: Keep Top 50 strongest cell-cell links
# -k 20: Keep Top 20 interactions per link
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/plot_cellchat_individual.R \
  -i /data/user3/sobj/cci/cellchat/pooled_secreted_strict_c1/samples/X2/cellchat.qs \
  -o /data/user3/sobj/cci/cellchat/pooled_secreted_strict_c1/plots_refined \
  -t 50 \
  -k 20

# Inspect Internal Data
Rscript /home/user3/data_user3/git_repo/_wt/cellchat/scripts/cellchat/inspect_cellchat.R \
  -i /data/user3/sobj/cci/cellchat/pooled_secreted_strict_c1/samples/X2/cellchat.qs \
  -n 20