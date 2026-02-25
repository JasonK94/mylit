#!/usr/bin/env Rscript
# Figure 3: DEG & FGS Prognostic Signature
# Panels A, B for the stroke scRNAseq paper
# Output: /data/user3/sobj/stroke_hc_v8_2/figures/v3_paper/fig3_*.png

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
})

outdir <- "/data/user3/sobj/stroke_hc_v8_2/figures/v3_paper"

# Publication theme
theme_pub <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

cohort_cols <- c("HC" = "#2166AC", "Stroke" = "#B2182B")
g3_cols <- c("1" = "#4DAF4A", "2" = "#E41A1C")

# ============================================================
# Panel A: DEG count per cell type (consensus)
# ============================================================
cat("Generating Panel A: DEG count bar chart...\n")

deg <- qread("/data/user3/sobj/stroke_hc_v8_2/deg/consensus/L1_cohort_anno1_1/results_consensus.qs")

# Extract consensus scores
cs <- deg$consensus_scores

# Count significant consensus DEGs per cluster
# Using meta_p_adj < 0.05 and agreement > 0 (i.e., at least some method agreement)
deg_counts <- cs %>%
  filter(meta_p_adj < 0.05 & n_significant >= 2) %>%
  mutate(
    direction = ifelse(weighted_beta > 0, "Up in Stroke", "Down in Stroke")
  ) %>%
  group_by(cluster_id, direction) %>%
  summarise(n_genes = n(), .groups = "drop")

# Also get total per cluster
deg_total <- deg_counts %>%
  group_by(cluster_id) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  arrange(desc(total))

# Clean cluster names
name_clean <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\+", "+", x)
  x <- gsub("/", "/", x)
  x
}

deg_counts$cell_type <- name_clean(deg_counts$cluster_id)
deg_total$cell_type <- name_clean(deg_total$cluster_id)

# Order by total count
ct_order <- deg_total %>% arrange(total) %>% pull(cell_type)
deg_counts$cell_type <- factor(deg_counts$cell_type, levels = ct_order)

p_a <- ggplot(deg_counts, aes(x = n_genes, y = cell_type, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c("Up in Stroke" = "#B2182B", "Down in Stroke" = "#2166AC"),
    name = "Direction"
  ) +
  # Add total count labels
  geom_text(
    data = deg_total,
    aes(x = total + 5, y = cell_type, label = total),
    inherit.aes = FALSE, size = 3.2, hjust = 0, color = "grey30"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  xlab("Number of consensus DEGs") +
  ylab("") +
  ggtitle("Consensus DEGs per cell type (L1: HC vs Stroke)") +
  theme_pub +
  theme(
    legend.position = c(0.75, 0.2),
    legend.background = element_rect(fill = "white", color = "grey80"),
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(outdir, "fig3_A_deg_count.png"),
  p_a, width = 10, height = 8, dpi = 300, bg = "white"
)
cat("  Saved fig3_A_deg_count.png\n")

# ============================================================
# Panel B: FGS score UMAP
# ============================================================
cat("Loading Seurat object for FGS UMAP...\n")

# Load FGS sobj which has the score
fgs_sobj <- qread("/data/user3/sobj/stroke_hc_v8_2/unbiased/fgs/whole_IS_50_sobj.qs")
cat("  FGS sobj loaded:", ncol(fgs_sobj), "cells (IS only)\n")

# Check available score columns
score_cols <- grep("fgs|score", colnames(fgs_sobj@meta.data), value = TRUE, ignore.case = TRUE)
cat("  Score columns:", paste(score_cols, collapse = ", "), "\n")

# Use fgs_meta_score
if ("fgs_meta_score" %in% colnames(fgs_sobj@meta.data)) {
  cat("  Using fgs_meta_score\n")

  p_b_fgs <- FeaturePlot(
    fgs_sobj,
    features = "fgs_meta_score",
    reduction = "umap.scvi",
    raster = FALSE,
    pt.size = 0.1,
    order = TRUE
  ) +
    scale_color_gradientn(
      colors = c("#2166AC", "#F7F7F7", "#B2182B"),
      name = "FGS score",
      limits = c(-3, 3),
      oob = squish
    ) +
    ggtitle("FGS prognostic signature score (IS patients)") +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_pub +
    theme(
      legend.position = "right",
      axis.line = element_line(linewidth = 0.5)
    )

  ggsave(
    file.path(outdir, "fig3_B_fgs_umap.png"),
    p_b_fgs, width = 10, height = 8, dpi = 300, bg = "white"
  )
  cat("  Saved fig3_B_fgs_umap.png\n")
} else {
  cat("  WARNING: fgs_meta_score not found. Skipping Panel B.\n")
}

rm(fgs_sobj); gc()

# ============================================================
# Panel C: FGS score by g3 outcome (violin)
# ============================================================
cat("Generating Panel C: FGS score by g3 outcome...\n")

# Reload FGS sobj for g3 analysis
fgs_sobj <- qread("/data/user3/sobj/stroke_hc_v8_2/unbiased/fgs/whole_IS_50_sobj.qs")

# Filter to g3 patients only
g3_data <- fgs_sobj@meta.data %>%
  filter(!is.na(g3) & g3 %in% c(1, 2)) %>%
  mutate(
    g3_label = ifelse(g3 == 1, "Good (g3=1)", "Bad (g3=2)"),
    g3_factor = factor(g3)
  )

if (nrow(g3_data) > 0 && "fgs_meta_score" %in% colnames(g3_data)) {
  # Wilcoxon test
  wt <- wilcox.test(fgs_meta_score ~ g3_factor, data = g3_data)
  p_label <- ifelse(wt$p.value < 2.2e-16, "p < 2.2e-16",
                     sprintf("p = %.2e", wt$p.value))

  p_c <- ggplot(g3_data, aes(x = g3_label, y = fgs_meta_score, fill = g3_factor)) +
    geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(
      values = c("1" = "#4DAF4A", "2" = "#E41A1C"),
      guide = "none"
    ) +
    annotate("text", x = 1.5, y = max(g3_data$fgs_meta_score, na.rm = TRUE) + 0.3,
             label = p_label, size = 4, fontface = "italic") +
    annotate("segment", x = 1, xend = 2,
             y = max(g3_data$fgs_meta_score, na.rm = TRUE) + 0.1,
             yend = max(g3_data$fgs_meta_score, na.rm = TRUE) + 0.1,
             linewidth = 0.5) +
    xlab("") +
    ylab("FGS prognostic score") +
    ggtitle("FGS score by 3-month outcome") +
    theme_pub

  ggsave(
    file.path(outdir, "fig3_C_fgs_violin_g3.png"),
    p_c, width = 6, height = 6, dpi = 300, bg = "white"
  )
  cat("  Saved fig3_C_fgs_violin_g3.png\n")
}

rm(fgs_sobj); gc()

# ============================================================
# Panel D: FGS score by g3 outcome - patient-level
# ============================================================
cat("Generating Panel D: Patient-level FGS score...\n")

fgs_sobj <- qread("/data/user3/sobj/stroke_hc_v8_2/unbiased/fgs/whole_IS_50_sobj.qs")

patient_scores <- fgs_sobj@meta.data %>%
  filter(!is.na(g3) & g3 %in% c(1, 2)) %>%
  group_by(patient_name, g3) %>%
  summarise(
    mean_fgs = mean(fgs_meta_score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    g3_label = ifelse(g3 == 1, "Good (g3=1)", "Bad (g3=2)"),
    g3_factor = factor(g3)
  )

if (nrow(patient_scores) > 0) {
  wt2 <- wilcox.test(mean_fgs ~ g3_factor, data = patient_scores)
  p_label2 <- sprintf("p = %.3f", wt2$p.value)

  p_d <- ggplot(patient_scores, aes(x = g3_label, y = mean_fgs, fill = g3_factor)) +
    geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
    geom_jitter(aes(size = n_cells), width = 0.15, alpha = 0.7, shape = 21, stroke = 0.5) +
    scale_fill_manual(
      values = c("1" = "#4DAF4A", "2" = "#E41A1C"),
      guide = "none"
    ) +
    scale_size_continuous(name = "# cells", range = c(2, 6)) +
    annotate("text", x = 1.5, y = max(patient_scores$mean_fgs) + 0.15,
             label = p_label2, size = 4, fontface = "italic") +
    annotate("segment", x = 1, xend = 2,
             y = max(patient_scores$mean_fgs) + 0.05,
             yend = max(patient_scores$mean_fgs) + 0.05,
             linewidth = 0.5) +
    xlab("") +
    ylab("Mean FGS score per patient") +
    ggtitle("Patient-level FGS prognostic score") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85))

  ggsave(
    file.path(outdir, "fig3_D_fgs_patient_g3.png"),
    p_d, width = 6, height = 6, dpi = 300, bg = "white"
  )
  cat("  Saved fig3_D_fgs_patient_g3.png\n")
}

rm(fgs_sobj); gc()

# ============================================================
# References to external validation panels (E-F)
# ============================================================
cat("\n--- External validation panels (reference only) ---\n")
ext_dir <- "/data/user3/sobj/stroke_hc_v8_2/deg/external_validation"
roc_file <- file.path(ext_dir, "03_validation_roc_curves.png")
auc_file <- file.path(ext_dir, "01_validation_auc_comparison.png")

if (file.exists(roc_file)) {
  file.copy(roc_file, file.path(outdir, "fig3_EF_validation_roc_curves.png"), overwrite = TRUE)
  cat("  Copied ROC curves -> fig3_EF_validation_roc_curves.png\n")
}
if (file.exists(auc_file)) {
  file.copy(auc_file, file.path(outdir, "fig3_EF_validation_auc_comparison.png"), overwrite = TRUE)
  cat("  Copied AUC comparison -> fig3_EF_validation_auc_comparison.png\n")
}

cat("Figure 3 complete!\n")
