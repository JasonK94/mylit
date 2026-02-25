#!/usr/bin/env Rscript
# Figure 5: Cell-Cell Communication
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
library(patchwork)
library(grid)

# -- Paths -------------------------------------------------------------------
base_dir    <- "/data/user3/sobj/stroke_hc_v8_2"
out_dir     <- file.path(base_dir, "figures/v3_paper")
cci_dir     <- file.path(base_dir, "cci/plots")
cc_dir      <- file.path(cci_dir, "cellchat_L1_cohort_v2")
integ_dir   <- file.path(cci_dir, "integrated_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Colors ------------------------------------------------------------------
col_cohort <- c("HC" = "#2166AC", "Stroke" = "#B2182B")
method_colors <- c(
  "CellChat" = "#E41A1C",
  "MNN"      = "#377EB8",
  "LIANA"    = "#4DAF4A",
  "CC"       = "#E41A1C",
  "CC+MNN"   = "#984EA3",
  "CC+LIANA" = "#FF7F00",
  "MNN+LIANA"= "#A65628",
  "CC+MNN+LIANA" = "#333333"
)
n_methods_colors <- c("1" = "#FEE0D2", "2" = "#FC9272", "3" = "#DE2D26")

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

cat("=== Figure 5: Cell-Cell Communication ===\n\n")

# ============================================================================
# PANEL A: CellChat differential interaction heatmap (reference existing)
# ============================================================================
cat("--- Panel A: CellChat diff interaction heatmap ---\n")
# Use diff_heatmap_count (count-based, more interpretable)
src_a <- file.path(cc_dir, "05_diff_heatmap_count.png")
if (file.exists(src_a)) {
  file.copy(src_a, file.path(out_dir, "fig5_A_cellchat_diff_heatmap.png"), overwrite = TRUE)
  cat("  Copied: fig5_A_cellchat_diff_heatmap.png (from CellChat v2)\n")
} else {
  cat("  WARNING: Missing source:", src_a, "\n")
}

# ============================================================================
# PANEL B: CellChat rankNet (reference existing)
# ============================================================================
cat("\n--- Panel B: CellChat rankNet ---\n")
src_b <- file.path(cc_dir, "07_rankNet_stacked.png")
if (file.exists(src_b)) {
  file.copy(src_b, file.path(out_dir, "fig5_B_cellchat_rankNet.png"), overwrite = TRUE)
  cat("  Copied: fig5_B_cellchat_rankNet.png (from CellChat v2)\n")
} else {
  cat("  WARNING: Missing source:", src_b, "\n")
}

# ============================================================================
# PANEL C: Top consensus pathways (NEW - publication quality)
# ============================================================================
cat("\n--- Panel C: Consensus pathway bar plot ---\n")

pw <- read.csv(file.path(integ_dir, "pathway_evidence_L1.csv"), stringsAsFactors = FALSE)
cat("  Loaded", nrow(pw), "pathways\n")

# Select top 20 pathways by combined_score
pw_top <- pw %>%
  arrange(desc(combined_score)) %>%
  head(20) %>%
  mutate(
    pathway = factor(pathway, levels = rev(pathway)),
    n_methods_f = factor(n_methods_support, levels = c("1", "2", "3"))
  )

# Determine which methods support each pathway
pw_top <- pw_top %>%
  mutate(
    has_cc   = !is.na(cc_score) & cc_score > 0,
    has_mnn  = !is.na(mnn_mean_score) & mnn_mean_score > 0,
    has_liana = !is.na(liana_mean_neg_log_rank) & liana_mean_neg_log_rank > 0,
    method_label = paste0(
      ifelse(has_cc, "CC", ""),
      ifelse(has_cc & (has_mnn | has_liana), "+", ""),
      ifelse(has_mnn, "MNN", ""),
      ifelse(has_mnn & has_liana, "+", ""),
      ifelse(has_liana, "LIANA", "")
    )
  )

# Main bar plot
p_c1 <- ggplot(pw_top, aes(x = combined_score, y = pathway)) +
  geom_col(aes(fill = n_methods_f), width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = n_methods_colors,
    name = "Methods\nsupporting",
    drop = FALSE
  ) +
  labs(
    title = "Top 20 CCI pathways by combined evidence score",
    subtitle = "Integration of CellChat, MultiNicheNet, and LIANA (L1: HC vs Stroke)",
    x = "Combined evidence score",
    y = ""
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position.inside = c(0.85, 0.25),
    legend.position = "inside",
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

# Add method indicator strip on the right
# Create a small table of method support for the dot panel
dot_data <- pw_top %>%
  select(pathway, has_cc, has_mnn, has_liana) %>%
  pivot_longer(cols = c(has_cc, has_mnn, has_liana),
               names_to = "method", values_to = "present") %>%
  mutate(
    method = recode(method, "has_cc" = "CellChat", "has_mnn" = "MNN", "has_liana" = "LIANA"),
    method = factor(method, levels = c("CellChat", "MNN", "LIANA"))
  )

p_c2 <- ggplot(dot_data, aes(x = method, y = pathway)) +
  geom_point(aes(color = present), size = 3) +
  scale_color_manual(values = c("TRUE" = "#333333", "FALSE" = "grey85"),
                     guide = "none") +
  labs(x = "", y = "", title = "Method\nsupport") +
  theme_pub +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(10, 5, 10, 0)
  )

# Combine
p_c <- p_c1 + p_c2 + plot_layout(widths = c(4, 1))

ggsave(file.path(out_dir, "fig5_C_pathway_barplot.png"), p_c,
       width = 11, height = 8, dpi = 300, bg = "white")
cat("  Saved: fig5_C_pathway_barplot.png\n")

# ============================================================================
# PANEL D: Consensus interactions dot plot (NEW - publication quality)
# ============================================================================
cat("\n--- Panel D: Consensus interactions dot plot ---\n")

cons <- read.csv(file.path(integ_dir, "consensus_interactions_L1_ge2methods.csv"),
                 stringsAsFactors = FALSE)
cat("  Loaded", nrow(cons), "consensus interactions (>=2 methods)\n")

# Create descriptive label
cons <- cons %>%
  mutate(
    interaction_label = paste0(ligand, " -> ", receptor),
    pair_label = paste0(sender, " -> ", receiver),
    n_methods_f = factor(n_methods, levels = c("2", "3"))
  )

# Sort by n_methods then by interaction_key
cons <- cons %>%
  arrange(desc(n_methods), interaction_label)

# Dot plot: sender-receiver pairs, size by n_methods
p_d <- ggplot(cons, aes(x = interaction_label, y = pair_label)) +
  geom_point(aes(size = n_methods, color = methods), alpha = 0.8) +
  scale_size_continuous(range = c(4, 8), breaks = c(2, 3), name = "N methods") +
  scale_color_manual(
    values = c("CC+LIANA" = "#FF7F00", "CC+MNN" = "#984EA3", "MNN+LIANA" = "#A65628"),
    name = "Method\ncombination"
  ) +
  labs(
    title = "Consensus CCI interactions (>=2 methods concordant, L1)",
    subtitle = paste0(nrow(cons), " interactions between cell type pairs"),
    x = "Ligand -> Receptor",
    y = "Sender -> Receiver"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

ggsave(file.path(out_dir, "fig5_D_consensus_interactions.png"), p_d,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("  Saved: fig5_D_consensus_interactions.png\n")

# Also create a variant focusing on the sender-receiver heatmap view
# Count interactions per sender-receiver pair
sr_counts <- cons %>%
  group_by(sender, receiver) %>%
  summarise(
    n_interactions = n(),
    ligands = paste(unique(ligand), collapse = ", "),
    receptors = paste(unique(receptor), collapse = ", "),
    .groups = "drop"
  )

p_d2 <- ggplot(sr_counts, aes(x = receiver, y = sender, fill = n_interactions)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_interactions), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#FEE0D2", high = "#CB181D", name = "N consensus\ninteractions") +
  labs(
    title = "Consensus interactions: sender-receiver pairs (L1, >=2 methods)",
    x = "Receiver cell type",
    y = "Sender cell type"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_equal()

ggsave(file.path(out_dir, "fig5_D2_sender_receiver_heatmap.png"), p_d2,
       width = 6, height = 5, dpi = 300, bg = "white")
cat("  Saved: fig5_D2_sender_receiver_heatmap.png\n")

# ============================================================================
# PANEL E: Pathway direction (reference existing)
# ============================================================================
cat("\n--- Panel E: Pathway direction ---\n")
src_e <- file.path(integ_dir, "09_pathway_direction_L1.png")
if (file.exists(src_e)) {
  file.copy(src_e, file.path(out_dir, "fig5_E_pathway_direction.png"), overwrite = TRUE)
  cat("  Copied: fig5_E_pathway_direction.png\n")
} else {
  cat("  WARNING: Missing source:", src_e, "\n")
}

# ============================================================================
# PANEL F: Method correlation (reference existing)
# ============================================================================
cat("\n--- Panel F: Method correlation ---\n")
src_f <- file.path(integ_dir, "04_pairwise_method_scatter_L1.png")
if (file.exists(src_f)) {
  file.copy(src_f, file.path(out_dir, "fig5_F_method_correlation.png"), overwrite = TRUE)
  cat("  Copied: fig5_F_method_correlation.png\n")
} else {
  cat("  WARNING: Missing source:", src_f, "\n")
}

# ============================================================================
# SUPPLEMENTARY: Additional referenced panels
# ============================================================================
cat("\n--- Supplementary: Additional CCI panels ---\n")

supp_refs <- list(
  "fig5_supp_upset_L1.png"     = file.path(integ_dir, "01_upset_interaction_concordance_L1.png"),
  "fig5_supp_upset_gene_L1.png" = file.path(integ_dir, "01b_upset_gene_concordance_L1.png"),
  "fig5_supp_consensus_top30.png" = file.path(integ_dir, "02_consensus_interactions_L1_top30.png"),
  "fig5_supp_method_summary.png"  = file.path(integ_dir, "08_method_summary_panel.png"),
  "fig5_supp_cellchat_bubble.png" = file.path(cc_dir, "09_bubble_comparison.png"),
  "fig5_supp_cellchat_scatter.png"= file.path(cc_dir, "12_signaling_scatter.png"),
  "fig5_supp_diff_circle_count.png" = file.path(cc_dir, "03_diff_circle_count.png")
)

for (dest_name in names(supp_refs)) {
  src <- supp_refs[[dest_name]]
  if (file.exists(src)) {
    file.copy(src, file.path(out_dir, dest_name), overwrite = TRUE)
    cat("  Copied:", dest_name, "\n")
  } else {
    cat("  WARNING: Missing:", src, "\n")
  }
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Figure 5 Complete ===\n")
cat("Generated panels:\n")
cat("  A:  fig5_A_cellchat_diff_heatmap.png (referenced from CellChat v2)\n")
cat("  B:  fig5_B_cellchat_rankNet.png (referenced from CellChat v2)\n")
cat("  C:  fig5_C_pathway_barplot.png (NEW - publication quality)\n")
cat("  D:  fig5_D_consensus_interactions.png (NEW - publication quality)\n")
cat("  D2: fig5_D2_sender_receiver_heatmap.png (NEW - publication quality)\n")
cat("  E:  fig5_E_pathway_direction.png (referenced from integrated summary)\n")
cat("  F:  fig5_F_method_correlation.png (referenced from integrated summary)\n")
cat("\nSupplementary panels also copied.\n")
cat("Output directory:", out_dir, "\n")
