#!/usr/bin/env Rscript
# Figure 4: Monocyte Trajectory & Gene Dynamics
# Publication-quality panels for stroke scRNAseq paper
# ============================================================================

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(qs)

# -- Paths -------------------------------------------------------------------
base_dir    <- "/data/user3/sobj/stroke_hc_v8_2"
out_dir     <- file.path(base_dir, "figures/v3_paper")
traj_dir    <- file.path(base_dir, "trajectory/v3")
analysis_dir <- file.path(traj_dir, "analysis")
v1_fig_dir  <- file.path(base_dir, "figures/v1_initial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Colors ------------------------------------------------------------------
col_cohort <- c("HC" = "#2166AC", "Stroke" = "#B2182B")
col_g3     <- c("1" = "#4DAF4A", "2" = "#E41A1C")

# -- Theme -------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 10),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 10),
    strip.text    = element_text(size = 11, face = "bold"),
    plot.margin   = margin(10, 10, 10, 10)
  )

cat("=== Figure 4: Monocyte Trajectory & Gene Dynamics ===\n\n")

# ============================================================================
# PANEL A: Monocyte UMAP with pseudotime
# ============================================================================
cat("--- Panel A: Monocyte pseudotime UMAP ---\n")

# Load mono subset from v2 (has umap.scvi coordinates + metadata)
mono_subset <- qread(file.path(base_dir, "trajectory/v2/mono/subset.qs"))
cat("  Loaded mono subset:", ncol(mono_subset), "cells\n")

# Check for pseudotime metadata in Seurat object
pt_cols <- grep("pseudotime|slingshot|sling", colnames(mono_subset@meta.data),
                value = TRUE, ignore.case = TRUE)
cat("  Pseudotime columns in Seurat:", paste(pt_cols, collapse = ", "), "\n")

# Get UMAP coordinates from Seurat
umap_data <- as.data.frame(Seurat::Embeddings(mono_subset, reduction = "umap.scvi"))
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
umap_data$cohort <- mono_subset@meta.data$cohort
umap_data$cell_id <- rownames(umap_data)

# Try to get pseudotime from Seurat metadata first
pt_col <- NULL
if ("slingPseudotime_1" %in% colnames(mono_subset@meta.data)) {
  pt_col <- "slingPseudotime_1"
  umap_data$pseudotime <- mono_subset@meta.data[[pt_col]]
} else if ("pseudotime" %in% colnames(mono_subset@meta.data)) {
  pt_col <- "pseudotime"
  umap_data$pseudotime <- mono_subset@meta.data[[pt_col]]
} else if (length(pt_cols) > 0) {
  pt_col <- pt_cols[1]
  umap_data$pseudotime <- mono_subset@meta.data[[pt_col]]
}

# If not in Seurat, load from monocle3 CDS
if (is.null(pt_col)) {
  cat("  Pseudotime not in Seurat metadata; loading monocle3 CDS...\n")
  library(monocle3)
  cds_file <- file.path(base_dir, "trajectory/v2/mono/monocle3_cds.rds")
  if (file.exists(cds_file)) {
    cds <- readRDS(cds_file)
    pt_vals <- monocle3::pseudotime(cds)
    cat("  CDS cells:", length(pt_vals), " | Seurat cells:", nrow(umap_data), "\n")
    # Match by cell names
    common <- intersect(names(pt_vals), umap_data$cell_id)
    cat("  Common cells:", length(common), "\n")
    umap_data$pseudotime <- pt_vals[umap_data$cell_id]
    rm(cds); gc()
  } else {
    cat("  WARNING: monocle3 CDS not found at:", cds_file, "\n")
  }
}

# Also try slingshot SCE for pseudotime
if (all(is.na(umap_data$pseudotime))) {
  cat("  Trying slingshot SCE...\n")
  sce_file <- file.path(base_dir, "trajectory/v2/mono/slingshot_sce.rds")
  if (file.exists(sce_file)) {
    sce <- readRDS(sce_file)
    if ("slingPseudotime_1" %in% colnames(SingleCellExperiment::colData(sce))) {
      pt_vals <- SingleCellExperiment::colData(sce)$slingPseudotime_1
      names(pt_vals) <- colnames(sce)
      common <- intersect(names(pt_vals), umap_data$cell_id)
      cat("  Slingshot common cells:", length(common), "\n")
      umap_data$pseudotime <- pt_vals[umap_data$cell_id]
    }
    rm(sce); gc()
  }
}

rm(mono_subset); gc()

# Plot if pseudotime available
if (!all(is.na(umap_data$pseudotime))) {
  umap_plot <- umap_data[!is.na(umap_data$pseudotime), ]
  cat("  Cells with pseudotime:", nrow(umap_plot), "\n")
  umap_plot <- umap_plot[order(umap_plot$pseudotime), ]

  p_a <- ggplot(umap_plot, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_viridis_c(option = "inferno", name = "Pseudotime") +
    labs(title = "Monocyte trajectory (pseudotime)",
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub +
    theme(legend.position = "right") +
    guides(color = guide_colorbar(barwidth = 1, barheight = 8))

  ggsave(file.path(out_dir, "fig4_A_pseudotime_umap.png"), p_a,
         width = 7, height = 5.5, dpi = 300, bg = "white")
  cat("  Saved: fig4_A_pseudotime_umap.png\n")

  # Also create split by cohort version
  p_a2 <- ggplot(umap_plot, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_viridis_c(option = "inferno", name = "Pseudotime") +
    facet_wrap(~cohort) +
    labs(title = "Monocyte pseudotime by cohort",
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub +
    theme(legend.position = "right") +
    guides(color = guide_colorbar(barwidth = 1, barheight = 8))

  ggsave(file.path(out_dir, "fig4_A2_pseudotime_umap_split.png"), p_a2,
         width = 12, height = 5.5, dpi = 300, bg = "white")
  cat("  Saved: fig4_A2_pseudotime_umap_split.png\n")
} else {
  # Fall back to existing plot
  src_a <- file.path(v1_fig_dir, "22_pseudotime_umap_mono.png")
  if (file.exists(src_a)) {
    file.copy(src_a, file.path(out_dir, "fig4_A_pseudotime_umap.png"), overwrite = TRUE)
    cat("  Copied existing pseudotime UMAP from v1\n")
  } else {
    cat("  WARNING: Could not generate Panel A (no pseudotime data)\n")
  }
}

rm(umap_data); gc()

# ============================================================================
# PANELS B-C: Pseudotime density (reference existing high-quality plots)
# ============================================================================
cat("\n--- Panels B-C: Pseudotime density ---\n")
# These already exist in v1_initial and are of sufficient quality.
# We copy them to v3_paper with standardized naming.
existing_density <- list(
  "fig4_B_density_cohort.png" = file.path(v1_fig_dir, "24_pseudotime_density_mono_cohort.png"),
  "fig4_C_density_g3.png"     = file.path(v1_fig_dir, "25_pseudotime_density_mono_g3.png")
)

for (dest_name in names(existing_density)) {
  src <- existing_density[[dest_name]]
  if (file.exists(src)) {
    file.copy(src, file.path(out_dir, dest_name), overwrite = TRUE)
    cat("  Copied:", dest_name, "(from existing v1 figure)\n")
  } else {
    cat("  WARNING: Missing source:", src, "\n")
  }
}

# ============================================================================
# PANELS D-E: Gene dynamics curves (reference existing high-quality plots)
# ============================================================================
cat("\n--- Panels D-E: Gene dynamics curves ---\n")
existing_curves <- list(
  "fig4_D_gene_dynamics_cohort.png" = file.path(analysis_dir, "diff_curves_mono_cohort.png"),
  "fig4_E_gene_dynamics_g3.png"     = file.path(analysis_dir, "diff_curves_mono_g3.png")
)

for (dest_name in names(existing_curves)) {
  src <- existing_curves[[dest_name]]
  if (file.exists(src)) {
    file.copy(src, file.path(out_dir, dest_name), overwrite = TRUE)
    cat("  Copied:", dest_name, "(from existing analysis)\n")
  } else {
    cat("  WARNING: Missing source:", src, "\n")
  }
}

# Also copy top genes plots for supplementary reference
for (comp in c("mono_cohort", "mono_g3")) {
  src <- file.path(analysis_dir, paste0("top_genes_", comp, ".png"))
  if (file.exists(src)) {
    file.copy(src, file.path(out_dir, paste0("fig4_supp_top_genes_", comp, ".png")),
              overwrite = TRUE)
    cat("  Copied supplementary: top_genes_", comp, "\n")
  }
}

# ============================================================================
# PANEL F: Effect size heatmap (rABC values)
# ============================================================================
cat("\n--- Panel F: Effect size heatmap ---\n")

# Read effect sizes
es <- read.csv(file.path(analysis_dir, "effect_sizes_all.csv"), stringsAsFactors = FALSE)
cat("  Total entries:", nrow(es), "\n")

# Filter to successful fits only
es_ok <- es %>% filter(status == "success")
cat("  Successful fits:", nrow(es_ok), "\n")

# Create compartment_analysis label
es_ok <- es_ok %>%
  mutate(
    condition = paste0(compartment, "/", analysis),
    # Determine significance: interaction_padj < 0.05
    is_sig = !is.na(interaction_padj) & interaction_padj < 0.05
  )

# Get unique genes and conditions
conditions <- c("mono/cohort", "mono/g3", "cd4/cohort", "cd4/g3", "cd8/cohort", "cd8/g3")

# For each gene, compute max rABC across conditions
gene_max <- es_ok %>%
  group_by(gene) %>%
  summarise(max_rABC = max(rABC, na.rm = TRUE),
            n_sig = sum(is_sig, na.rm = TRUE),
            mean_rABC = mean(rABC, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(max_rABC))

# Get top 30 genes by max rABC
top_genes <- head(gene_max$gene, 30)
cat("  Top 30 genes selected\n")

# Build matrix: genes x conditions
mat <- matrix(NA, nrow = length(top_genes), ncol = length(conditions),
              dimnames = list(top_genes, conditions))
sig_mat <- matrix(FALSE, nrow = length(top_genes), ncol = length(conditions),
                  dimnames = list(top_genes, conditions))

for (i in seq_len(nrow(es_ok))) {
  g <- es_ok$gene[i]
  c <- es_ok$condition[i]
  if (g %in% top_genes && c %in% conditions) {
    mat[g, c] <- es_ok$rABC[i]
    sig_mat[g, c] <- es_ok$is_sig[i]
  }
}

# Read concordance data to identify 3-method-confirmed genes
conc_mono_cohort <- read.csv(file.path(base_dir, "trajectory/v3/lamian/mono/cohort/concordance_merged.csv"),
                             stringsAsFactors = FALSE)
# Also read pseudobulk results for the 3-method overlap
pb_mono_cohort <- read.csv(file.path(base_dir, "trajectory/v3/pseudobulk/mono/cohort/pb_gamm_results.csv"),
                           stringsAsFactors = FALSE)

# Identify genes confirmed by GAMM + Lamian (both sig in mono/cohort)
# GAMM sig: interaction_padj < 0.05
gamm_sig_genes <- es_ok %>%
  filter(compartment == "mono", analysis == "cohort", is_sig) %>%
  pull(gene)

# Lamian sig: fdr.chisq.overall < 0.05
lamian_sig_genes <- conc_mono_cohort %>%
  filter(!is.na(fdr.chisq.overall) & fdr.chisq.overall < 0.05) %>%
  pull(gene)

# PB GAMM sig
pb_sig_genes <- if (nrow(pb_mono_cohort) > 0) {
  pb_col <- intersect(c("fdr", "p.adjust", "padj", "interaction_padj"),
                      colnames(pb_mono_cohort))
  if (length(pb_col) > 0) {
    pb_mono_cohort %>%
      filter(.data[[pb_col[1]]] < 0.05) %>%
      pull(gene)
  } else {
    character(0)
  }
} else {
  character(0)
}

# 3-method confirmed: GAMM + Lamian + PB
three_method_genes <- intersect(intersect(gamm_sig_genes, lamian_sig_genes), pb_sig_genes)
# If PB data sparse, fall back to 2-method: GAMM + Lamian
two_method_genes <- intersect(gamm_sig_genes, lamian_sig_genes)

cat("  GAMM sig genes (mono/cohort):", length(gamm_sig_genes), "\n")
cat("  Lamian sig genes (mono/cohort):", length(lamian_sig_genes), "\n")
cat("  PB GAMM sig genes:", length(pb_sig_genes), "\n")
cat("  3-method confirmed:", length(three_method_genes), "-", paste(three_method_genes, collapse = ", "), "\n")
cat("  2-method confirmed:", length(two_method_genes), "-", paste(two_method_genes, collapse = ", "), "\n")

# Use the larger set for annotation
confirmed_genes <- if (length(three_method_genes) >= 5) three_method_genes else two_method_genes
n_method_label <- if (length(three_method_genes) >= 5) "3-method" else "2-method (GAMM+Lamian)"

# Replace NAs with 0 for clustering (genes not tested in a condition get 0)
mat_for_plot <- mat
mat_for_plot[is.na(mat_for_plot)] <- 0

# Create color function - cap at 95th percentile for better contrast
rABC_cap <- quantile(mat, 0.95, na.rm = TRUE)
col_fun <- colorRamp2(
  c(0, rABC_cap / 2, rABC_cap),
  c("#2166AC", "white", "#B2182B")
)

# Column labels
col_labels <- c(
  "mono/cohort" = "Mono\nHC vs Stroke",
  "mono/g3"     = "Mono\ng3=1 vs g3=2",
  "cd4/cohort"  = "CD4 T\nHC vs Stroke",
  "cd4/g3"      = "CD4 T\ng3=1 vs g3=2",
  "cd8/cohort"  = "CD8 T\nHC vs Stroke",
  "cd8/g3"      = "CD8 T\ng3=1 vs g3=2"
)

# Row annotation: confirmed by multiple methods
is_confirmed <- top_genes %in% confirmed_genes
row_ha <- rowAnnotation(
  "Confirmed" = anno_simple(
    ifelse(is_confirmed, "Yes", "No"),
    col = c("Yes" = "#FFD700", "No" = "grey90"),
    width = unit(4, "mm")
  ),
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 9)
)

# Column annotation: compartment
comp_colors <- c("mono" = "#E41A1C", "cd4" = "#377EB8", "cd8" = "#4DAF4A")
comp_labels <- gsub("/.*", "", conditions)
col_ha <- HeatmapAnnotation(
  "Compartment" = comp_labels,
  "Layer" = gsub(".*/", "", conditions),
  col = list(
    "Compartment" = comp_colors,
    "Layer" = c("cohort" = "#984EA3", "g3" = "#FF7F00")
  ),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(4, "mm"),
  gap = unit(1, "mm")
)

# Cell function: add * for significant, mark NAs
cell_fun <- function(j, i, x, y, width, height, fill) {
  if (is.na(mat[i, j])) {
    # Mark originally NA cells (gene not tested in this condition)
    grid.rect(x, y, width, height,
              gp = gpar(fill = "grey90", col = NA))
    grid.text("-", x, y, gp = gpar(fontsize = 7, col = "grey60"))
  } else if (!is.na(sig_mat[i, j]) && sig_mat[i, j]) {
    grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
  }
}

# Draw heatmap
png(file.path(out_dir, "fig4_F_effect_size_heatmap.png"),
    width = 8, height = 10, units = "in", res = 300, bg = "white")

ht <- Heatmap(
  mat_for_plot,
  name = "rABC",
  col = col_fun,
  na_col = "grey95",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  row_dend_width = unit(15, "mm"),
  column_labels = col_labels[conditions],
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10,
                       fontface = ifelse(top_genes %in% confirmed_genes, "bold", "plain")),
  cell_fun = cell_fun,
  top_annotation = col_ha,
  right_annotation = row_ha,
  row_title = "Gene",
  column_title = "Relative ABC (rABC) — Gene Dynamics Effect Size",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    title = "rABC",
    at = seq(0, ceiling(rABC_cap), by = 0.5),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  border = TRUE
)

draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))

# Add footer annotation
decorate_heatmap_body("rABC", {
  grid.text(
    paste0("* p_adj < 0.05; bold gene names = ", n_method_label, " confirmed"),
    x = 0.5, y = unit(-3, "mm"), just = "top",
    gp = gpar(fontsize = 8, fontface = "italic")
  )
})

dev.off()
cat("  Saved: fig4_F_effect_size_heatmap.png\n")

# ============================================================================
# PANEL G: Compartment summary & concordance
# ============================================================================
cat("\n--- Panel G: Compartment summary & concordance ---\n")

# Read compartment summary
comp_summ <- read.csv(file.path(analysis_dir, "compartment_summary.csv"), stringsAsFactors = FALSE)

# Read concordance summary
conc_summ <- read.csv(file.path(base_dir, "trajectory/v3/lamian/concordance_summary.csv"),
                      stringsAsFactors = FALSE)

# Create a clean bar plot of n_sig genes and Jaccard index per compartment/analysis
comp_summ$label <- paste0(comp_summ$compartment, "\n", comp_summ$analysis)
comp_summ$label <- factor(comp_summ$label,
                          levels = paste0(rep(c("mono", "cd4", "cd8"), each = 2), "\n",
                                         rep(c("cohort", "g3"), 3)))

# Merge concordance data
conc_summ$label <- paste0(conc_summ$compartment, "\n", conc_summ$analysis)
comp_merged <- merge(comp_summ, conc_summ[, c("label", "both_sig", "jaccard")],
                     by = "label", all.x = TRUE)

# Panel G1: Number of significant genes
p_g1 <- ggplot(comp_merged, aes(x = label, y = n_sig, fill = label)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_sig), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = rep(c("#984EA3", "#FF7F00"), 3),
                    guide = "none") +
  labs(title = "Significant gene dynamics (interaction p_adj < 0.05)",
       x = "", y = "Number of significant genes") +
  theme_pub +
  theme(axis.text.x = element_text(size = 9))

# Panel G2: Jaccard concordance (GAMM vs Lamian)
conc_data <- conc_summ %>%
  filter(!is.na(jaccard)) %>%
  mutate(label_clean = paste0(compartment, "/", analysis))

p_g2 <- ggplot(conc_data, aes(x = reorder(label_clean, -jaccard), y = jaccard,
                               fill = jaccard)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), vjust = -0.3, size = 3.5) +
  scale_fill_gradient(low = "#FEE0D2", high = "#CB181D", guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "GAMM-Lamian concordance (Jaccard index)",
       x = "", y = "Jaccard index") +
  theme_pub +
  coord_cartesian(ylim = c(0, 0.85))

# Combine G1 and G2
library(patchwork)
p_g <- p_g1 + p_g2 + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = list(c("G1", "G2")),
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(out_dir, "fig4_G_compartment_concordance.png"), p_g,
       width = 12, height = 5, dpi = 300, bg = "white")
cat("  Saved: fig4_G_compartment_concordance.png\n")

# Also reference existing concordance plot
conc_src <- file.path(base_dir, "trajectory/v3/lamian/concordance_jaccard.png")
if (file.exists(conc_src)) {
  file.copy(conc_src, file.path(out_dir, "fig4_G_ref_concordance_jaccard.png"), overwrite = TRUE)
  cat("  Copied reference: concordance_jaccard.png\n")
}

# Also reference compartment summary
comp_src <- file.path(analysis_dir, "compartment_summary.png")
if (file.exists(comp_src)) {
  file.copy(comp_src, file.path(out_dir, "fig4_G_ref_compartment_summary.png"), overwrite = TRUE)
  cat("  Copied reference: compartment_summary.png\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Figure 4 Complete ===\n")
cat("Generated panels:\n")
cat("  A:  fig4_A_pseudotime_umap.png (new)\n")
cat("  A2: fig4_A2_pseudotime_umap_split.png (new)\n")
cat("  B:  fig4_B_density_cohort.png (referenced from v1)\n")
cat("  C:  fig4_C_density_g3.png (referenced from v1)\n")
cat("  D:  fig4_D_gene_dynamics_cohort.png (referenced from v3 analysis)\n")
cat("  E:  fig4_E_gene_dynamics_g3.png (referenced from v3 analysis)\n")
cat("  F:  fig4_F_effect_size_heatmap.png (NEW - publication quality)\n")
cat("  G:  fig4_G_compartment_concordance.png (NEW - publication quality)\n")
cat("\nOutput directory:", out_dir, "\n")
