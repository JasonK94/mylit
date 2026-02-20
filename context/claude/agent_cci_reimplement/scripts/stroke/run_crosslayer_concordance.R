#!/usr/bin/env Rscript
# ==============================================================================
# Cross-Layer Concordance Analysis
# Compares DEG results across Layer 1 (HC vs Stroke) and Layer 2 (g3 within IS)
# Also integrates FGS gene signatures
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(qs)
})

cat("=== Cross-Layer Concordance Analysis ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "figures")
DEG_BASE <- file.path(BASE, "deg_consensus")

# ==============================================================================
# PART 1: Load DEG L1 results
# ==============================================================================
cat("[1/5] Loading DEG L1 results...\n")

deg_l1_dir <- file.path(DEG_BASE, "L1_cohort_anno1_1")

# Load consensus results
l1_consensus <- qread(file.path(deg_l1_dir, "results_consensus.qs"))

# Load per-method results
l1_edger <- read.csv(file.path(deg_l1_dir, "method_results", "muscat-edgeR_results.csv"))
l1_deseq2 <- read.csv(file.path(deg_l1_dir, "method_results", "muscat-DESeq2_results.csv"))
l1_nebula <- read.csv(file.path(deg_l1_dir, "method_results", "nebula_results.csv"))

# Standardize column names: muscat uses p_adj.loc, NEBULA uses FDR
standardize_deg <- function(df) {
  if ("p_adj.loc" %in% names(df) && !"FDR" %in% names(df)) {
    df$FDR <- df$p_adj.loc
  }
  if ("p_val" %in% names(df) && !"pvalue" %in% names(df)) {
    df$pvalue <- df$p_val
  }
  df
}

l1_edger <- standardize_deg(l1_edger)
l1_deseq2 <- standardize_deg(l1_deseq2)
l1_nebula <- standardize_deg(l1_nebula)

cat(sprintf("  edgeR: %d genes, DESeq2: %d genes, NEBULA: %d genes\n",
            nrow(l1_edger), nrow(l1_deseq2), nrow(l1_nebula)))

# ==============================================================================
# PART 2: DEG summary statistics per cluster
# ==============================================================================
cat("[2/5] Summarizing DEG results per cluster...\n")

summarize_degs <- function(df, method_name, fdr_cut = 0.05, lfc_cut = 0.25) {
  df %>%
    filter(!is.na(FDR) & !is.na(logFC)) %>%
    group_by(cluster_id) %>%
    summarise(
      total_tested = n(),
      sig_up = sum(FDR < fdr_cut & logFC > lfc_cut),
      sig_down = sum(FDR < fdr_cut & logFC < -lfc_cut),
      sig_total = sig_up + sig_down,
      method = method_name,
      .groups = "drop"
    )
}

deg_summary <- bind_rows(
  summarize_degs(l1_edger, "muscat-edgeR"),
  summarize_degs(l1_deseq2, "muscat-DESeq2"),
  summarize_degs(l1_nebula, "NEBULA")
)

# Bar plot of DEG counts per cluster per method
p_deg_count <- ggplot(deg_summary, aes(x = reorder(cluster_id, sig_total), y = sig_total, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "DEG Count per Cluster (L1: Stroke vs HC)",
       subtitle = "FDR < 0.05, |logFC| > 0.25",
       x = "", y = "Number of DEGs", fill = "Method") +
  theme_minimal()

ggsave(file.path(OUT, "32_deg_count_per_cluster_L1.png"), p_deg_count,
       width = 12, height = 8, dpi = 200)
cat("  Saved: 32_deg_count_per_cluster_L1\n")

# ==============================================================================
# PART 3: Method concordance (logFC scatter between methods)
# ==============================================================================
cat("[3/5] Method concordance...\n")

# Merge edgeR and NEBULA by gene+cluster
merged_en <- inner_join(
  l1_edger %>% select(gene, cluster_id, logFC_edgeR = logFC, FDR_edgeR = FDR),
  l1_nebula %>% select(gene, cluster_id, logFC_nebula = logFC, FDR_nebula = FDR),
  by = c("gene", "cluster_id")
) %>% filter(!is.na(logFC_edgeR) & !is.na(logFC_nebula))

# Top clusters by cell count
top_clusters <- merged_en %>%
  count(cluster_id) %>%
  arrange(desc(n)) %>%
  head(6) %>%
  pull(cluster_id)

p_scatter_methods <- ggplot(merged_en %>% filter(cluster_id %in% top_clusters),
                             aes(logFC_edgeR, logFC_nebula)) +
  geom_point(aes(color = ifelse(FDR_edgeR < 0.05 & FDR_nebula < 0.05, "Both sig",
                                 ifelse(FDR_edgeR < 0.05 | FDR_nebula < 0.05, "One sig", "NS"))),
             size = 0.3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Both sig" = "red", "One sig" = "orange", "NS" = "grey70")) +
  facet_wrap(~cluster_id, scales = "free") +
  labs(title = "Method Concordance: muscat-edgeR vs NEBULA",
       x = "logFC (muscat-edgeR)", y = "logFC (NEBULA)", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT, "33_method_concordance_edgeR_nebula.png"), p_scatter_methods,
       width = 14, height = 10, dpi = 200)
cat("  Saved: 33_method_concordance_edgeR_nebula\n")

# ==============================================================================
# PART 4: Cross-layer overlap (L1 vs L2 DEGs)
# ==============================================================================
cat("[4/5] Cross-layer concordance...\n")

# Check if L2 results exist
l2_dir_candidates <- list.files(DEG_BASE, pattern = "L2_g3_anno1", full.names = TRUE)
l2_dir <- NULL
for (d in l2_dir_candidates) {
  if (file.exists(file.path(d, "method_results", "nebula_results.csv"))) {
    l2_dir <- d
    break
  }
}

if (!is.null(l2_dir)) {
  cat("  L2 DEG results found, computing cross-layer concordance...\n")

  l2_nebula <- read.csv(file.path(l2_dir, "method_results", "nebula_results.csv"))
  l2_nebula <- standardize_deg(l2_nebula)
  cat(sprintf("  L2 NEBULA: %d genes\n", nrow(l2_nebula)))

  # Find shared clusters
  shared_clusters <- intersect(unique(l1_nebula$cluster_id), unique(l2_nebula$cluster_id))
  cat(sprintf("  Shared clusters: %d\n", length(shared_clusters)))

  # Merge L1 and L2 NEBULA results
  cross <- inner_join(
    l1_nebula %>% filter(cluster_id %in% shared_clusters) %>%
      select(gene, cluster_id, logFC_L1 = logFC, FDR_L1 = FDR),
    l2_nebula %>% filter(cluster_id %in% shared_clusters) %>%
      select(gene, cluster_id, logFC_L2 = logFC, FDR_L2 = FDR),
    by = c("gene", "cluster_id")
  ) %>% filter(!is.na(logFC_L1) & !is.na(logFC_L2))

  cat(sprintf("  Cross-layer gene pairs: %d\n", nrow(cross)))

  # Cross-layer logFC scatter
  cross$category <- case_when(
    cross$FDR_L1 < 0.05 & cross$FDR_L2 < 0.05 ~ "Both layers",
    cross$FDR_L1 < 0.05 ~ "L1 only (HC vs IS)",
    cross$FDR_L2 < 0.05 ~ "L2 only (g3)",
    TRUE ~ "NS"
  )

  for (cl in head(shared_clusters, 6)) {
    cl_data <- cross %>% filter(cluster_id == cl)
    if (nrow(cl_data) < 10) next

    p_cross <- ggplot(cl_data, aes(logFC_L1, logFC_L2, color = category)) +
      geom_point(size = 0.5, alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("Both layers" = "red", "L1 only (HC vs IS)" = "blue",
                                     "L2 only (g3)" = "green", "NS" = "grey80")) +
      labs(title = sprintf("Cross-Layer: %s", cl),
           x = "logFC (L1: Stroke vs HC)", y = "logFC (L2: g3=2 vs g3=1)") +
      theme_minimal()
    # Don't save individual, combine below
  }

  # Global cross-layer scatter (all clusters)
  p_cross_global <- ggplot(cross, aes(logFC_L1, logFC_L2, color = category)) +
    geom_point(size = 0.2, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Both layers" = "red", "L1 only (HC vs IS)" = "blue",
                                   "L2 only (g3)" = "green", "NS" = "grey80")) +
    facet_wrap(~cluster_id, scales = "free") +
    labs(title = "Cross-Layer DEG Concordance: L1 (Stroke vs HC) × L2 (g3=2 vs g3=1)",
         x = "logFC (L1)", y = "logFC (L2)") +
    theme_minimal() +
    theme(legend.position = "bottom", strip.text = element_text(size = 7))

  ggsave(file.path(OUT, "34_cross_layer_deg_scatter.png"), p_cross_global,
         width = 18, height = 14, dpi = 200)
  cat("  Saved: 34_cross_layer_deg_scatter\n")

  # Cross-layer concordance statistics
  concordance_stats <- cross %>%
    filter(FDR_L1 < 0.05 | FDR_L2 < 0.05) %>%
    group_by(cluster_id) %>%
    summarise(
      both_sig = sum(FDR_L1 < 0.05 & FDR_L2 < 0.05),
      l1_only = sum(FDR_L1 < 0.05 & FDR_L2 >= 0.05),
      l2_only = sum(FDR_L1 >= 0.05 & FDR_L2 < 0.05),
      concordant = sum(FDR_L1 < 0.05 & FDR_L2 < 0.05 & sign(logFC_L1) == sign(logFC_L2)),
      discordant = sum(FDR_L1 < 0.05 & FDR_L2 < 0.05 & sign(logFC_L1) != sign(logFC_L2)),
      .groups = "drop"
    )

  write.csv(concordance_stats, file.path(OUT, "35_cross_layer_concordance_stats.csv"), row.names = FALSE)
  cat("  Saved: 35_cross_layer_concordance_stats.csv\n")

} else {
  cat("  L2 DEG results not yet available. Will skip cross-layer scatter.\n")
  cat("  Run this script again after L2 DEG completes.\n")
}

# ==============================================================================
# PART 5: FGS gene overlap with DEGs
# ==============================================================================
cat("[5/5] FGS × DEG overlap...\n")

fgs_genes <- read.csv(file.path(BASE, "fgs/whole_IS_50_cmgi_genes.csv"))
fgs_top50 <- head(fgs_genes$gene, 50)
fgs_all <- fgs_genes$gene

# Check FGS genes in L1 DEG results
fgs_in_deg <- l1_nebula %>%
  filter(gene %in% fgs_all) %>%
  filter(!is.na(FDR) & !is.na(logFC))

fgs_sig <- fgs_in_deg %>%
  filter(FDR < 0.05, abs(logFC) > 0.25)

cat(sprintf("  FGS genes found in DEG: %d / %d\n", length(unique(fgs_in_deg$gene)), length(fgs_all)))
cat(sprintf("  FGS genes significant (FDR<0.05, |lfc|>0.25): %d\n", length(unique(fgs_sig$gene))))

# Heatmap: FGS gene logFC across clusters
fgs_lfc <- fgs_in_deg %>%
  filter(gene %in% fgs_top50) %>%
  select(gene, cluster_id, logFC) %>%
  pivot_wider(names_from = cluster_id, values_from = logFC) %>%
  as.data.frame()

if (nrow(fgs_lfc) > 3) {
  rownames(fgs_lfc) <- fgs_lfc$gene
  fgs_lfc$gene <- NULL
  mat <- as.matrix(fgs_lfc)
  mat[is.na(mat)] <- 0

  # Plot as ggplot heatmap
  fgs_long <- fgs_in_deg %>%
    filter(gene %in% fgs_top50) %>%
    mutate(sig = ifelse(FDR < 0.05, "*", ""))

  p_fgs_deg <- ggplot(fgs_long, aes(x = cluster_id, y = gene, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = sig), size = 3, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-2, 2),
                          oob = scales::squish) +
    labs(title = "FGS Top 50 Genes: logFC in L1 DEG (Stroke vs HC)",
         subtitle = "* = FDR < 0.05",
         x = "", y = "", fill = "logFC") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))

  ggsave(file.path(OUT, "36_fgs_genes_in_deg_heatmap.png"), p_fgs_deg,
         width = 16, height = 12, dpi = 200)
  cat("  Saved: 36_fgs_genes_in_deg_heatmap\n")
}

# Save summary table
fgs_deg_summary <- fgs_in_deg %>%
  group_by(gene) %>%
  summarise(
    n_clusters_tested = n(),
    n_sig_clusters = sum(FDR < 0.05 & abs(logFC) > 0.25),
    mean_logFC = mean(logFC, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(fgs_genes %>% select(gene, importance = abs_contribution), by = "gene") %>%
  arrange(desc(n_sig_clusters), min_FDR)

write.csv(fgs_deg_summary, file.path(OUT, "37_fgs_deg_overlap_summary.csv"), row.names = FALSE)
cat("  Saved: 37_fgs_deg_overlap_summary.csv\n")

cat("\n=== CROSS-LAYER CONCORDANCE COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
