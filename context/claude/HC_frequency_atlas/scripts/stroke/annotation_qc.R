#!/usr/bin/env Rscript
# ==============================================================================
# Annotation QC: Centroid distances, canonical markers, co-expression analysis
# ==============================================================================
# Output: /data/user3/sobj/stroke_hc_v8_2/results/annotation_qc/
#
# Strategy:
# 1. Cluster centroid distances in scVI latent space
# 2. Canonical marker absolute expression (pct.exp + avg_exp)
# 3. Mutually exclusive marker co-expression for doublet detection
# 4. g3 / cohort group distribution

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
  library(Matrix)
})

cat("=== Annotation QC Script ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# --- 1. Load data ---
cat("[1/6] Loading annotated object...\n")
sobj <- qread("/data/user3/sobj/stroke_hc_v8_2/4_strokev8_annotated.qs")
cat("  Cells:", ncol(sobj), "Genes:", nrow(sobj), "\n")

outdir <- "/data/user3/sobj/stroke_hc_v8_2/results/annotation_qc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- 2. Cluster centroid distances in scVI latent space ---
cat("\n[2/6] Computing cluster centroid distances (scVI latent space)...\n")
scvi_emb <- Embeddings(sobj, "integrated.scvi")  # cells x 30 latent dims
cluster_ids <- sobj$seurat_clusters

# Compute centroids
centroids <- tapply(seq_len(nrow(scvi_emb)), cluster_ids, function(idx) {
  colMeans(scvi_emb[idx, , drop = FALSE])
})
centroid_mat <- do.call(rbind, centroids)
rownames(centroid_mat) <- names(centroids)

# Pairwise euclidean distances
centroid_dist <- as.matrix(dist(centroid_mat))

# Save
write.csv(centroid_dist, file.path(outdir, "centroid_distances_scvi.csv"))
cat("  Saved centroid distance matrix:", nrow(centroid_dist), "x", ncol(centroid_dist), "\n")

# For each cluster, find nearest neighbor
nn_info <- data.frame(
  cluster = rownames(centroid_dist),
  anno1 = sapply(rownames(centroid_dist), function(cl) {
    unique(sobj$anno1[sobj$seurat_clusters == cl])[1]
  }),
  anno2 = sapply(rownames(centroid_dist), function(cl) {
    unique(sobj$anno2[sobj$seurat_clusters == cl])[1]
  }),
  n_cells = as.numeric(table(cluster_ids)[rownames(centroid_dist)]),
  nearest_cluster = sapply(rownames(centroid_dist), function(cl) {
    dists <- centroid_dist[cl, ]
    dists[cl] <- Inf
    names(which.min(dists))
  }),
  nearest_distance = sapply(rownames(centroid_dist), function(cl) {
    dists <- centroid_dist[cl, ]
    dists[cl] <- Inf
    min(dists)
  }),
  stringsAsFactors = FALSE
)
nn_info <- nn_info[order(as.numeric(nn_info$cluster)), ]
write.csv(nn_info, file.path(outdir, "cluster_nearest_neighbors.csv"), row.names = FALSE)
cat("  Nearest neighbor info saved\n")

# --- 3. Canonical marker expression (absolute, not relative) ---
cat("\n[3/6] Evaluating canonical marker expression...\n")

canonical_markers <- list(
  T_cell = c("CD3D", "CD3E", "CD3G"),
  CD4 = c("CD4", "IL7R"),
  CD8 = c("CD8A", "CD8B"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  NK = c("NKG7", "GNLY", "KLRF1", "NCAM1"),
  B_cell = c("MS4A1", "CD79A"),
  Plasma = c("JCHAIN", "MZB1"),
  Mono_classical = c("CD14", "LYZ", "S100A8", "S100A9"),
  Mono_nonclassical = c("FCGR3A", "MS4A7"),
  DC = c("FCER1A", "CD1C", "CLEC9A"),
  pDC = c("LILRA4", "CLEC4C"),
  Platelet = c("PPBP", "PF4"),
  RBC = c("HBB", "HBA1"),
  Mast = c("TPSAB1", "CPA3"),
  Proliferating = c("MKI67", "TOP2A"),
  ISG = c("ISG15", "MX1", "IFI44L"),
  Neutrophil_like = c("CSF3R", "CXCR2"),
  Macrophage = c("C1QA", "APOE")
)

# Filter to genes present in the object
all_markers <- unique(unlist(canonical_markers))
present_markers <- all_markers[all_markers %in% rownames(sobj)]
missing_markers <- setdiff(all_markers, present_markers)
if (length(missing_markers) > 0) {
  cat("  Missing markers:", paste(missing_markers, collapse = ", "), "\n")
}

# DotPlot data extraction (pct.exp + avg.exp per cluster)
cat("  Computing per-cluster expression stats for", length(present_markers), "markers...\n")

# Use Seurat's internal method for efficiency
dp <- DotPlot(sobj, features = present_markers, group.by = "seurat_clusters")
dp_data <- dp$data
write.csv(dp_data, file.path(outdir, "canonical_marker_dotplot_data.csv"), row.names = FALSE)
cat("  DotPlot data saved\n")

# Also compute per anno1
dp_anno1 <- DotPlot(sobj, features = present_markers, group.by = "anno1")
write.csv(dp_anno1$data, file.path(outdir, "canonical_marker_dotplot_anno1.csv"), row.names = FALSE)

# --- 4. Doublet co-expression analysis ---
cat("\n[4/6] Evaluating mutually exclusive marker co-expression (doublet detection)...\n")

doublet_pairs <- list(
  T_Mono = c("CD3D", "CD14"),
  T_B = c("CD3D", "MS4A1"),
  T_Mono2 = c("CD3D", "LYZ"),
  B_Mono = c("MS4A1", "CD14"),
  T_PLT = c("CD3D", "PPBP"),
  Mono_PLT = c("CD14", "PPBP"),   # PMA - keep these!
  NK_Mono = c("NKG7", "CD14"),
  T_pDC = c("CD3D", "LILRA4"),
  T_RBC = c("CD3D", "HBB"),
  Leukocyte_RBC = c("PTPRC", "HBB")
)

# For each cluster, compute fraction of cells co-expressing both markers
coexp_results <- list()
for (pair_name in names(doublet_pairs)) {
  genes <- doublet_pairs[[pair_name]]
  if (!all(genes %in% rownames(sobj))) next
  
  expr <- GetAssayData(sobj, layer = "data")[genes, , drop = FALSE]
  # Co-expression: both > 0
  both_pos <- (expr[1, ] > 0) & (expr[2, ] > 0)
  gene1_pos <- expr[1, ] > 0
  gene2_pos <- expr[2, ] > 0
  
  per_cluster <- tapply(seq_along(both_pos), sobj$seurat_clusters, function(idx) {
    n <- length(idx)
    data.frame(
      n_cells = n,
      pct_both = round(100 * sum(both_pos[idx]) / n, 2),
      pct_gene1 = round(100 * sum(gene1_pos[idx]) / n, 2),
      pct_gene2 = round(100 * sum(gene2_pos[idx]) / n, 2)
    )
  })
  
  df <- do.call(rbind, per_cluster)
  df$cluster <- rownames(df)
  df$pair <- pair_name
  df$gene1 <- genes[1]
  df$gene2 <- genes[2]
  coexp_results[[pair_name]] <- df
}

coexp_df <- do.call(rbind, coexp_results)
rownames(coexp_df) <- NULL
write.csv(coexp_df, file.path(outdir, "doublet_coexpression.csv"), row.names = FALSE)
cat("  Co-expression analysis saved for", length(doublet_pairs), "marker pairs\n")

# Flag suspicious clusters (high co-expression of mutually exclusive markers)
# Exclude Mono_PLT as PMA is expected
suspicious <- coexp_df %>%
  filter(pair != "Mono_PLT") %>%
  filter(pct_both > 15) %>%
  arrange(desc(pct_both))
write.csv(suspicious, file.path(outdir, "suspicious_coexpression.csv"), row.names = FALSE)
cat("  Suspicious clusters (>15% co-expression):", nrow(suspicious), "entries\n")

# --- 5. Group distribution ---
cat("\n[5/6] Computing group distributions...\n")

# Per-cluster cohort distribution
cohort_dist <- sobj@meta.data %>%
  group_by(seurat_clusters, anno1, anno2, cohort) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  mutate(total = HC + Stroke, pct_HC = round(100 * HC / total, 1), pct_Stroke = round(100 * Stroke / total, 1))
write.csv(cohort_dist, file.path(outdir, "cluster_cohort_distribution.csv"), row.names = FALSE)

# Per-cluster g3 distribution (IS only)
g3_dist <- sobj@meta.data %>%
  filter(index_injury_hand == "IS", !is.na(g3)) %>%
  group_by(seurat_clusters, anno1, anno2, g3) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = g3, values_from = n, values_fill = 0)
colnames(g3_dist)[colnames(g3_dist) %in% c("1", "2")] <- c("g3_good", "g3_bad")
g3_dist <- g3_dist %>%
  mutate(total_IS = g3_good + g3_bad,
         pct_good = round(100 * g3_good / total_IS, 1),
         pct_bad = round(100 * g3_bad / total_IS, 1))
write.csv(g3_dist, file.path(outdir, "cluster_g3_distribution.csv"), row.names = FALSE)

# Patient-level summary
patient_summary <- sobj@meta.data %>%
  group_by(patient_name, cohort, index_injury_hand, g3, sex, age, GEM) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(cohort, g3)
write.csv(patient_summary, file.path(outdir, "patient_summary.csv"), row.names = FALSE)

cat("  IS patients with g3:", nrow(patient_summary %>% filter(!is.na(g3), g3 %in% c(1,2))), "\n")
cat("  g3==1 (good):", sum(patient_summary$g3 == 1, na.rm = TRUE), "patients\n")
cat("  g3==2 (bad):", sum(patient_summary$g3 == 2, na.rm = TRUE), "patients\n")
cat("  HC:", sum(patient_summary$cohort == "HC"), "subjects\n")

# --- 6. Summary report ---
cat("\n[6/6] Generating summary...\n")

# Overall cell counts
cat("\n=== OVERALL STATS ===\n")
cat("Total cells:", ncol(sobj), "\n")
cat("HC cells:", sum(sobj$cohort == "HC"), "\n")
cat("Stroke cells:", sum(sobj$cohort == "Stroke"), "\n")
cat("IS cells with g3:", sum(sobj$index_injury_hand == "IS" & !is.na(sobj$g3) & sobj$g3 %in% c(1,2)), "\n")

cat("\n=== ANNOTATION DISTRIBUTION ===\n")
anno_table <- table(sobj$anno1)
print(sort(anno_table, decreasing = TRUE))

cat("\nCompleted:", format(Sys.time()), "\n")
cat("All outputs saved to:", outdir, "\n")
