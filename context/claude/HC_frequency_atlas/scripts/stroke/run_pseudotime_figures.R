#!/usr/bin/env Rscript
# ==============================================================================
# Pseudotime FeaturePlot overlays, density distributions, and signature-sorted stacked bars
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(monocle3)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

cat("=== Pseudotime & Signature Figures ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ <- file.path(BASE, "trajectory_v2")
OUT <- file.path(BASE, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PART 1: Pseudotime FeaturePlot overlays on scVI UMAP
# ==============================================================================
cat("[1/4] Pseudotime FeaturePlot overlays...\n")

compartments <- c("mono", "cd4", "cd8")
compartment_labels <- c("Myeloid", "CD4+ T", "CD8+ T/NK")

for (i in seq_along(compartments)) {
  comp <- compartments[i]
  comp_label <- compartment_labels[i]
  cat(sprintf("  %s...\n", comp))

  sobj_path <- file.path(TRAJ, comp, "subset.qs")
  cds_path <- file.path(TRAJ, comp, "monocle3_cds.rds")

  if (!file.exists(sobj_path) || !file.exists(cds_path)) {
    cat(sprintf("  SKIP: Missing files for %s\n", comp))
    next
  }

  sobj <- qread(sobj_path)
  cds <- readRDS(cds_path)

  # Extract pseudotime from monocle3 CDS
  pt <- pseudotime(cds)
  pt[is.infinite(pt)] <- NA

  # Add pseudotime to Seurat object
  cells_with_pt <- intersect(names(pt), colnames(sobj))
  sobj$pseudotime <- NA_real_
  sobj$pseudotime[match(cells_with_pt, colnames(sobj))] <- pt[cells_with_pt]

  # Get scVI UMAP coordinates
  umap_coords <- Embeddings(sobj, "umap.scvi")
  plot_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    pseudotime = sobj$pseudotime,
    cohort = sobj$cohort,
    anno1 = sobj$anno1
  )

  # 1a. Pseudotime FeaturePlot
  p1 <- ggplot(plot_df[!is.na(plot_df$pseudotime), ], aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.1, alpha = 0.5) +
    scale_color_viridis_c(option = "inferno", na.value = "grey80") +
    labs(title = sprintf("Pseudotime: %s", comp_label), color = "PT") +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))

  # 1b. Split by cohort
  p2 <- ggplot(plot_df[!is.na(plot_df$pseudotime), ], aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.1, alpha = 0.5) +
    scale_color_viridis_c(option = "inferno", na.value = "grey80") +
    facet_wrap(~cohort) +
    labs(title = sprintf("Pseudotime: %s (by Cohort)", comp_label), color = "PT") +
    theme_minimal() +
    theme(legend.position = "right", strip.text = element_text(face = "bold"))

  ggsave(file.path(OUT, sprintf("22_pseudotime_umap_%s.png", comp)), p1,
         width = 8, height = 6, dpi = 200)
  ggsave(file.path(OUT, sprintf("23_pseudotime_umap_%s_cohort.png", comp)), p2,
         width = 14, height = 6, dpi = 200)
  cat(sprintf("  Saved: 22/23_pseudotime_umap_%s\n", comp))

  # 1c. Pseudotime density by cohort
  p3 <- ggplot(plot_df[!is.na(plot_df$pseudotime), ], aes(pseudotime, fill = cohort)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("HC" = "#2166AC", "Stroke" = "#B2182B")) +
    labs(title = sprintf("Pseudotime Distribution: %s", comp_label),
         x = "Pseudotime", y = "Density") +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, sprintf("24_pseudotime_density_%s_cohort.png", comp)), p3,
         width = 8, height = 5, dpi = 200)

  # 1d. Pseudotime density by g3 (IS only)
  if ("g3" %in% names(sobj@meta.data)) {
    plot_df_g3 <- plot_df
    plot_df_g3$g3 <- sobj$g3
    plot_df_g3 <- plot_df_g3[!is.na(plot_df_g3$pseudotime) & !is.na(plot_df_g3$g3), ]
    if (nrow(plot_df_g3) > 100) {
      plot_df_g3$g3 <- factor(plot_df_g3$g3, labels = c("Good (g3=1)", "Bad (g3=2)"))
      p4 <- ggplot(plot_df_g3, aes(pseudotime, fill = g3)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = c("Good (g3=1)" = "#4DAF4A", "Bad (g3=2)" = "#E41A1C")) +
        labs(title = sprintf("Pseudotime Distribution: %s (IS only)", comp_label),
             x = "Pseudotime", y = "Density") +
        theme_minimal() +
        theme(legend.position = "right")
      ggsave(file.path(OUT, sprintf("25_pseudotime_density_%s_g3.png", comp)), p4,
             width = 8, height = 5, dpi = 200)
    }
  }

  cat(sprintf("  Saved: 24/25_pseudotime_density_%s\n", comp))

  # 1e. Pseudotime boxplot by celltype
  plot_df_ct <- plot_df[!is.na(plot_df$pseudotime), ]
  if (nrow(plot_df_ct) > 100) {
    ct_order <- plot_df_ct %>%
      group_by(anno1) %>%
      summarise(med = median(pseudotime, na.rm = TRUE)) %>%
      arrange(med) %>%
      pull(anno1)
    plot_df_ct$anno1 <- factor(plot_df_ct$anno1, levels = ct_order)

    p5 <- ggplot(plot_df_ct, aes(anno1, pseudotime, fill = cohort)) +
      geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
      scale_fill_manual(values = c("HC" = "#2166AC", "Stroke" = "#B2182B")) +
      coord_flip() +
      labs(title = sprintf("Pseudotime by Cell Type: %s", comp_label),
           x = "", y = "Pseudotime") +
      theme_minimal()
    ggsave(file.path(OUT, sprintf("26_pseudotime_boxplot_%s.png", comp)), p5,
           width = 10, height = max(4, length(ct_order) * 0.4), dpi = 200)
    cat(sprintf("  Saved: 26_pseudotime_boxplot_%s\n", comp))
  }

  rm(sobj, cds); gc(verbose = FALSE)
}

# ==============================================================================
# PART 2: Stacked bar sorted by FGS signature score
# ==============================================================================
cat("\n[2/4] Stacked bars sorted by signature score...\n")

# Load the clean object
sobj <- qread(file.path(BASE, "5_strokev8_clean.qs"))
Idents(sobj) <- "anno1"

# Load FGS genes for module scoring
fgs_genes <- read.csv(file.path(BASE, "fgs/whole_IS_50_cmgi_genes.csv"))
top50 <- head(fgs_genes$gene, 50)
top_in_data <- intersect(top50, rownames(sobj))

# Add module score
sobj <- AddModuleScore(sobj, features = list(FGS_top50 = top_in_data), name = "FGS_score")

# Patient-level aggregation
patient_data <- sobj@meta.data %>%
  group_by(patient_name, cohort) %>%
  summarise(
    fgs_score = mean(FGS_score1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(fgs_score) %>%
  mutate(rank = row_number())

# Cell type proportions per patient
freq_data <- sobj@meta.data %>%
  group_by(patient_name, anno1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_name) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# Merge with patient ranking
freq_ranked <- freq_data %>%
  left_join(patient_data %>% select(patient_name, rank, cohort, fgs_score), by = "patient_name")

# Order patients by FGS score
freq_ranked$patient_name <- factor(freq_ranked$patient_name,
                                    levels = patient_data$patient_name)

# Define anno1 order: monocytes at bottom, T cells at top
# (User requested: Mono at bottom, T cell at top)
anno1_order <- c(
  "CD14+ Monocyte", "CD16+ Monocyte", "Inflammatory Monocyte",
  "CD4_S100A8_CSF3R", "CD4_S100A8_CD14",
  "cDC1", "cDC2", "pDC",
  "ISG+ Myeloid", "Proliferating",
  "Plasma_cell", "B_cell",
  "NK_cell",
  "MAIT",
  "ISG+ T_cell",
  "CD8+ Trm", "CD8+ T_Cytotoxic",
  "Treg",
  "CD4+ T_Naive/Memory",
  "Platelet/PLA"
)
# Use only those present in data
available_anno1 <- intersect(anno1_order, unique(freq_ranked$anno1))
missing_anno1 <- setdiff(unique(freq_ranked$anno1), anno1_order)
final_order <- c(available_anno1, missing_anno1)
freq_ranked$anno1 <- factor(freq_ranked$anno1, levels = rev(final_order))

# Color palette
n_types <- length(final_order)
colors <- scales::hue_pal()(n_types)
names(colors) <- final_order

p_stacked <- ggplot(freq_ranked, aes(x = patient_name, y = freq, fill = anno1)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors) +
  labs(title = "Cell Type Composition by Patient (sorted by FGS Signature Score)",
       subtitle = "Low FGS score (left) → High FGS score (right)",
       x = "Patients (sorted by FGS score)", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 7),
    plot.title = element_text(face = "bold")
  )

# Add cohort annotation bar at bottom
cohort_bar <- ggplot(patient_data, aes(x = factor(patient_name, levels = patient_name),
                                         y = 1, fill = cohort)) +
  geom_tile() +
  scale_fill_manual(values = c("HC" = "#2166AC", "Stroke" = "#B2182B")) +
  theme_void() +
  theme(legend.position = "none")

# Add FGS score line at top
score_line <- ggplot(patient_data, aes(x = factor(patient_name, levels = patient_name),
                                         y = fgs_score, group = 1)) +
  geom_line(color = "darkred", linewidth = 0.5) +
  geom_point(aes(color = cohort), size = 0.8) +
  scale_color_manual(values = c("HC" = "#2166AC", "Stroke" = "#B2182B")) +
  labs(y = "FGS Score") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), legend.position = "none")

combined <- score_line / p_stacked / cohort_bar +
  plot_layout(heights = c(1, 6, 0.3))

ggsave(file.path(OUT, "27_stacked_bar_fgs_score_sorted.png"), combined,
       width = 16, height = 10, dpi = 200)
ggsave(file.path(OUT, "27_stacked_bar_fgs_score_sorted.pdf"), combined,
       width = 16, height = 10)
cat("  Saved: 27_stacked_bar_fgs_score_sorted\n")

# Same but for anno2
freq_data2 <- sobj@meta.data %>%
  group_by(patient_name, anno2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_name) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  left_join(patient_data %>% select(patient_name, rank, cohort, fgs_score), by = "patient_name")

freq_data2$patient_name <- factor(freq_data2$patient_name,
                                    levels = patient_data$patient_name)
anno2_order <- c("CD14_Mono", "CD16_Mono", "DC", "B_cell", "Plasma",
                 "NK", "CD4_T", "CD8_T")
available_anno2 <- intersect(anno2_order, unique(freq_data2$anno2))
missing_anno2 <- setdiff(unique(freq_data2$anno2), anno2_order)
freq_data2$anno2 <- factor(freq_data2$anno2, levels = rev(c(available_anno2, missing_anno2)))

p_stacked2 <- ggplot(freq_data2, aes(x = patient_name, y = freq, fill = anno2)) +
  geom_bar(stat = "identity", width = 1) +
  labs(title = "Cell Type Composition (anno2, sorted by FGS Score)",
       x = "Patients (sorted by FGS score)", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right", plot.title = element_text(face = "bold"))

combined2 <- score_line / p_stacked2 / cohort_bar +
  plot_layout(heights = c(1, 6, 0.3))

ggsave(file.path(OUT, "28_stacked_bar_anno2_fgs_sorted.png"), combined2,
       width = 14, height = 8, dpi = 200)
cat("  Saved: 28_stacked_bar_anno2_fgs_sorted\n")

# ==============================================================================
# PART 3: Key gene FeaturePlots on scVI UMAP
# ==============================================================================
cat("\n[3/4] Key gene FeaturePlots...\n")

key_genes <- c("S100A8", "S100A9", "VCAN", "FCN1", "GZMB", "NKG7",
               "IL1B", "CXCL8", "ISG15", "FOXP3", "HIF1A", "CCL2")
key_in_data <- intersect(key_genes, rownames(sobj))

if (length(key_in_data) >= 4) {
  # Use FeaturePlot with scVI UMAP
  DefaultAssay(sobj) <- "RNA"
  p_feat <- FeaturePlot(sobj, features = key_in_data, reduction = "umap.scvi",
                         ncol = 4, pt.size = 0.05, order = TRUE) &
    scale_color_viridis_c(option = "magma") &
    theme(plot.title = element_text(size = 10))

  ggsave(file.path(OUT, "29_featureplot_key_genes.png"), p_feat,
         width = 20, height = ceiling(length(key_in_data) / 4) * 5, dpi = 200)
  cat(sprintf("  Saved: 29_featureplot_key_genes (%d genes)\n", length(key_in_data)))
}

# ==============================================================================
# PART 4: FGS score FeaturePlot on scVI UMAP
# ==============================================================================
cat("\n[4/4] FGS score FeaturePlot...\n")

p_fgs_umap <- FeaturePlot(sobj, features = "FGS_score1", reduction = "umap.scvi",
                            pt.size = 0.05, order = TRUE) +
  scale_color_viridis_c(option = "magma") +
  labs(title = "FGS Signature Score (Top 50 genes)", color = "Score")

p_fgs_split <- FeaturePlot(sobj, features = "FGS_score1", reduction = "umap.scvi",
                             pt.size = 0.05, order = TRUE, split.by = "cohort") &
  scale_color_viridis_c(option = "magma")

ggsave(file.path(OUT, "30_fgs_score_umap.png"), p_fgs_umap,
       width = 8, height = 6, dpi = 200)
ggsave(file.path(OUT, "31_fgs_score_umap_cohort.png"), p_fgs_split,
       width = 14, height = 6, dpi = 200)
cat("  Saved: 30/31_fgs_score_umap\n")

rm(sobj); gc(verbose = FALSE)
cat("\n=== PSEUDOTIME & SIGNATURE FIGURES COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
