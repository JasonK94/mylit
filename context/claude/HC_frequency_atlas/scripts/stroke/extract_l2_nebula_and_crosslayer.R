#!/usr/bin/env Rscript
# ==============================================================================
# Extract L2 NEBULA results from consensus.qs and run cross-layer concordance
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

cat("=== Extract L2 NEBULA & Cross-Layer Concordance ===\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "figures")

# Load L2 consensus results (contains NEBULA results even though consensus failed)
l2_res <- qread(file.path(BASE, "deg_consensus/L2_g3_anno1_1/results_consensus.qs"))

# Extract NEBULA results
cat("L2 result structure:\n")
cat("  Names:", paste(names(l2_res), collapse = ", "), "\n")

# Extract NEBULA results from the nested structure
if ("results" %in% names(l2_res) && "nebula" %in% names(l2_res$results)) {
  l2_nebula <- l2_res$results$nebula
  cat("  Found in l2_res$results$nebula\n")
} else if ("standardized_results" %in% names(l2_res) && "nebula" %in% names(l2_res$standardized_results)) {
  l2_nebula <- l2_res$standardized_results$nebula
  cat("  Found in l2_res$standardized_results$nebula\n")
} else {
  l2_nebula <- NULL
  cat("  Could not find NEBULA results\n")
}

if (is.null(l2_nebula) || (is.data.frame(l2_nebula) && nrow(l2_nebula) == 0)) {
  cat("ERROR: Could not extract L2 NEBULA results. Checking structure...\n")
  str(l2_res, max.level = 2)
  stop("Cannot proceed without L2 NEBULA results")
}

cat(sprintf("  L2 NEBULA extracted: %d rows\n", nrow(l2_nebula)))

# Standardize column names
if ("p_adj.loc" %in% names(l2_nebula) && !"FDR" %in% names(l2_nebula)) {
  l2_nebula$FDR <- l2_nebula$p_adj.loc
}
if ("pvalue" %in% names(l2_nebula) && !"FDR" %in% names(l2_nebula)) {
  l2_nebula$FDR <- p.adjust(l2_nebula$pvalue, method = "BH")
}

# Save as CSV for future use
dir.create(file.path(BASE, "deg_consensus/L2_g3_anno1_1/method_results"), showWarnings = FALSE)
write.csv(l2_nebula, file.path(BASE, "deg_consensus/L2_g3_anno1_1/method_results/nebula_results.csv"),
          row.names = FALSE)
cat("  Saved: L2 NEBULA results CSV\n")

# ---- Load L1 NEBULA results ----
l1_nebula <- read.csv(file.path(BASE, "deg_consensus/L1_cohort_anno1_1/method_results/nebula_results.csv"))
cat(sprintf("  L1 NEBULA: %d rows\n", nrow(l1_nebula)))

# ---- Cross-layer concordance ----
shared_clusters <- intersect(unique(l1_nebula$cluster_id), unique(l2_nebula$cluster_id))
cat(sprintf("  Shared clusters: %d\n", length(shared_clusters)))

cross <- inner_join(
  l1_nebula %>% filter(cluster_id %in% shared_clusters) %>%
    select(gene, cluster_id, logFC_L1 = logFC, FDR_L1 = FDR),
  l2_nebula %>% filter(cluster_id %in% shared_clusters) %>%
    select(gene, cluster_id, logFC_L2 = logFC, FDR_L2 = FDR),
  by = c("gene", "cluster_id")
) %>% filter(!is.na(logFC_L1) & !is.na(logFC_L2))

cat(sprintf("  Cross-layer gene pairs: %d\n", nrow(cross)))

cross$category <- case_when(
  cross$FDR_L1 < 0.05 & cross$FDR_L2 < 0.05 ~ "Both layers",
  cross$FDR_L1 < 0.05 ~ "L1 only (Stroke vs HC)",
  cross$FDR_L2 < 0.05 ~ "L2 only (g3=Bad vs Good)",
  TRUE ~ "NS"
)

cat("  Category counts:\n")
print(table(cross$category))

# Global cross-layer scatter
p_cross <- ggplot(cross, aes(logFC_L1, logFC_L2, color = category)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Both layers" = "red", "L1 only (Stroke vs HC)" = "#2166AC",
                                 "L2 only (g3=Bad vs Good)" = "#4DAF4A", "NS" = "grey80")) +
  facet_wrap(~cluster_id, scales = "free") +
  labs(title = "Cross-Layer DEG Concordance (NEBULA)",
       subtitle = "L1: Stroke vs HC | L2: Bad (g3=2) vs Good (g3=1) outcome",
       x = "logFC (L1: Stroke vs HC)", y = "logFC (L2: g3=2 vs g3=1)") +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))

ggsave(file.path(OUT, "34_cross_layer_deg_scatter.png"), p_cross,
       width = 18, height = 14, dpi = 200)
cat("  Saved: 34_cross_layer_deg_scatter\n")

# Concordance statistics
concordance_stats <- cross %>%
  filter(FDR_L1 < 0.05 | FDR_L2 < 0.05) %>%
  group_by(cluster_id) %>%
  summarise(
    both_sig = sum(FDR_L1 < 0.05 & FDR_L2 < 0.05),
    l1_only = sum(FDR_L1 < 0.05 & FDR_L2 >= 0.05),
    l2_only = sum(FDR_L1 >= 0.05 & FDR_L2 < 0.05),
    concordant = sum(FDR_L1 < 0.05 & FDR_L2 < 0.05 & sign(logFC_L1) == sign(logFC_L2)),
    discordant = sum(FDR_L1 < 0.05 & FDR_L2 < 0.05 & sign(logFC_L1) != sign(logFC_L2)),
    concordance_rate = round(concordant / max(1, concordant + discordant), 3),
    .groups = "drop"
  )

write.csv(concordance_stats, file.path(OUT, "35_cross_layer_concordance_stats.csv"), row.names = FALSE)
cat("  Saved: 35_cross_layer_concordance_stats\n")

# Top cross-layer concordant genes
top_concordant <- cross %>%
  filter(FDR_L1 < 0.05 & FDR_L2 < 0.05 & sign(logFC_L1) == sign(logFC_L2)) %>%
  mutate(combined_effect = abs(logFC_L1) + abs(logFC_L2)) %>%
  arrange(desc(combined_effect)) %>%
  head(50)

write.csv(top_concordant, file.path(OUT, "38_top_concordant_genes_L1_L2.csv"), row.names = FALSE)
cat(sprintf("  Saved: 38_top_concordant_genes (top %d)\n", nrow(top_concordant)))

cat("\n=== CROSS-LAYER COMPLETE ===\n")
