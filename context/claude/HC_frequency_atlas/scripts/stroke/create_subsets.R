#!/usr/bin/env Rscript
# ==============================================================================
# Create analysis subsets after annotation QC
# Depends on: annotation_qc.R results for cluster decisions
# ==============================================================================
# Input:  4_strokev8_annotated.qs
# Output: 5_strokev8_clean.qs      (LQ/RBC removed, PLT kept)
#         5_1_hc_is.qs             (HC + IS only, for Layer 1)
#         5_2_is_g3.qs             (IS only with valid g3, for Layer 2 & FGS)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
})

cat("=== Create Subsets Script ===\n")
cat("Started:", format(Sys.time()), "\n\n")

base_dir <- "/data/user3/sobj/stroke_hc_v8_2"
qc_dir <- file.path(base_dir, "results/annotation_qc")

# --- 1. Load annotated object ---
cat("[1/5] Loading annotated object...\n")
sobj <- qread(file.path(base_dir, "4_strokev8_annotated.qs"))
cat("  Original:", ncol(sobj), "cells\n")

# --- 2. Review QC results and decide on removals ---
cat("\n[2/5] Reviewing QC results...\n")

# Load suspicious co-expression data
suspicious <- read.csv(file.path(qc_dir, "suspicious_coexpression.csv"))
cat("  Suspicious co-expression entries:\n")
print(suspicious %>% select(cluster, pair, pct_both, n_cells) %>% head(20))

# Load centroid distances for ambiguous clusters
nn_info <- read.csv(file.path(qc_dir, "cluster_nearest_neighbors.csv"))
cat("\n  Ambiguous clusters (19, 25) nearest neighbors:\n")
print(nn_info %>% filter(cluster %in% c(19, 25)))

# --- 3. Define removal criteria ---
cat("\n[3/5] Applying removal and re-annotation criteria...\n")

# =====================================================================
# QC DECISION LOG (based on annotation_qc.R results, 2026-02-12):
#
# REMOVE (doublets / low quality / contaminants):
#   - Cluster 14 (T/NK_Low_quality_14, 6498 cells): Low quality
#   - Cluster 24 (Mono_Low_quality_24, 2810 cells): Low quality
#   - Cluster 31 (B_Low_quality_31, 519 cells): Low quality
#   - Cluster 27 (Erythrocyte, 1050 cells): RBC contamination
#   - Cluster 34 (Erythrocyte, 246 cells): RBC contamination
#   - Cluster 13 (Macrophage-like, 7435 cells): T/Mono DOUBLET
#       Evidence: CD3D+CD14=38%, CD3D+LYZ=48%, NKG7+CD14=46%
#   - Cluster 22 (Plasma_cell, 3020 cells): Multi-lineage DOUBLET
#       Evidence: MS4A1+CD14=29%, CD3D+MS4A1=24%, NKG7+CD14=23%
#
# KEEP & RE-ANNOTATE:
#   - Cluster 19 (CD4_S100A8_CD14 → CD4+T_S100A8_stress, 4735 cells):
#       Evidence: CD3D+CD14=0.02%, CD3D=69%, CD14=0.04%. NOT doublets.
#       S100A8 expression likely stress/activation response in T cells.
#   - Cluster 25 (CD4_S100A8_CSF3R → T/NK_activated, 2747 cells):
#       Evidence: CD3D+CD14=0%, NKG7=34%, CD3D=54%. NOT doublets.
#       Mixed T/NK phenotype with CSF3R expression.
#
# KEEP AS-IS:
#   - Cluster 28 (Platelet/PLA, 824 cells): Keep for PMA analysis
#   - Cluster 35 (Mast_cell, 22 cells): Keep (rare but real)
# =====================================================================

# Clusters to remove
remove_clusters <- c(14, 24, 31, 27, 34, 13, 22)
remove_anno1 <- c(
  "T/NK_Low_quality_14",
  "Mono_Low_quality_24",
  "B_Low_quality_31",
  "Erythrocyte",
  "Macrophage-like"   # cluster 13 - doublet
)

# Apply removals
cells_to_remove <- colnames(sobj)[
  sobj$anno1 %in% remove_anno1 |
  (sobj$seurat_clusters == 22)   # Plasma doublet cluster (anno1="Plasma_cell" but cluster 22 specifically)
]
cat("  Removing", length(cells_to_remove), "cells:\n")
for (a in unique(sobj$anno1[colnames(sobj) %in% cells_to_remove])) {
  cat("    ", a, ":", sum(sobj$anno1[colnames(sobj) %in% cells_to_remove] == a), "\n")
}
cat("    Cluster 22 (Plasma doublet):", sum(sobj$seurat_clusters[colnames(sobj) %in% cells_to_remove] == 22), "\n")

sobj_clean <- subset(sobj, cells = setdiff(colnames(sobj), cells_to_remove))
cat("  After removal:", ncol(sobj_clean), "cells\n")

# Re-annotate ambiguous clusters
sobj_clean$anno1_v2 <- sobj_clean$anno1
sobj_clean$anno1_v2[sobj_clean$seurat_clusters == 19] <- "CD4+T_S100A8_stress"
sobj_clean$anno1_v2[sobj_clean$seurat_clusters == 25] <- "T/NK_activated"
cat("  Re-annotated cluster 19 → CD4+T_S100A8_stress:", sum(sobj_clean$seurat_clusters == 19), "cells\n")
cat("  Re-annotated cluster 25 → T/NK_activated:", sum(sobj_clean$seurat_clusters == 25), "cells\n")

# Update anno2 to reflect changes
sobj_clean$anno2_v2 <- sobj_clean$anno2
# Cluster 19 stays Tc, cluster 25 stays Tc → no change needed for anno2

# --- 4. Create subsets ---
cat("\n[4/5] Creating subsets...\n")

# 5_strokev8_clean.qs: All clean cells (HC + all Stroke)
qsave(sobj_clean, file.path(base_dir, "5_strokev8_clean.qs"))
cat("  Saved 5_strokev8_clean.qs:", ncol(sobj_clean), "cells\n")

# 5_1_hc_is.qs: HC + IS only
sobj_hc_is <- subset(sobj_clean, index_injury_hand %in% c("HC", "IS"))
qsave(sobj_hc_is, file.path(base_dir, "5_1_hc_is.qs"))
cat("  Saved 5_1_hc_is.qs:", ncol(sobj_hc_is), "cells\n")
cat("    HC:", sum(sobj_hc_is$cohort == "HC"), "cells\n")
cat("    IS:", sum(sobj_hc_is$index_injury_hand == "IS"), "cells\n")

# 5_2_is_g3.qs: IS only with valid g3 (for FGS/TML)
sobj_is <- subset(sobj_clean, index_injury_hand == "IS" & !is.na(g3) & g3 %in% c(1, 2))
# Ensure g3 is a clean factor
sobj_is$g3 <- factor(sobj_is$g3, levels = c(1, 2))
qsave(sobj_is, file.path(base_dir, "5_2_is_g3.qs"))
cat("  Saved 5_2_is_g3.qs:", ncol(sobj_is), "cells\n")
cat("    g3==1 (good):", sum(sobj_is$g3 == 1), "cells,",
    length(unique(sobj_is$patient_name[sobj_is$g3 == 1])), "patients\n")
cat("    g3==2 (bad):", sum(sobj_is$g3 == 2), "cells,",
    length(unique(sobj_is$patient_name[sobj_is$g3 == 2])), "patients\n")

# --- 5. Summary ---
cat("\n[5/5] Summary...\n")
cat("  anno1 distribution in IS g3 subset:\n")
print(sort(table(sobj_is$anno1), decreasing = TRUE))
cat("\n  anno2 distribution in IS g3 subset:\n")
print(sort(table(sobj_is$anno2), decreasing = TRUE))

cat("\nCompleted:", format(Sys.time()), "\n")
