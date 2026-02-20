#!/usr/bin/env Rscript
# ==============================================================================
# Descriptive Figures Part 2: Frequency boxplots + FGS signature + key genes
# Continuation from Part 1 (UMAP, DotPlot, stacked bars already done)
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(ggpubr)
})

# Source myR plot functions
source("/home/user3/data_user3/git_repo/mylit/myR/R/utils_aggregation.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots_box.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots_heatmap.R")

cat("=== Descriptive Figures Part 2 ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "figures")

# ---- Load data ----
cat("[1/6] Loading data...\n")
sobj <- qread(file.path(BASE, "5_strokev8_clean.qs"))
cat("  Cells:", ncol(sobj), "\n")
Idents(sobj) <- "anno1"

# ---- Figure 5: Frequency boxplots ----
cat("\n[2/6] Frequency boxplots...\n")

tryCatch({
  p_freq_anno2 <- pb_freq(sobj, cluster = "anno2", group = "cohort",
                           sample = "patient_name")
  ggsave(file.path(OUT, "09_freq_boxplot_anno2_cohort.pdf"), p_freq_anno2,
         width = 12, height = 6)
  ggsave(file.path(OUT, "09_freq_boxplot_anno2_cohort.png"), p_freq_anno2,
         width = 12, height = 6, dpi = 200)
  cat("  Saved: 09_freq_boxplot_anno2_cohort\n")
}, error = function(e) {
  cat("  pb_freq failed:", e$message, "\n")
  cat("  Falling back to manual frequency calculation...\n")

  # Manual frequency boxplot
  meta <- sobj@meta.data
  freq_df <- meta %>%
    group_by(patient_name, cohort, anno2) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient_name) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()

  p <- ggplot(freq_df, aes(x = anno2, y = freq, fill = cohort)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
               size = 0.8, alpha = 0.5) +
    stat_compare_means(aes(group = cohort), method = "wilcox.test",
                       label = "p.signif", size = 3) +
    theme_bw() +
    labs(x = "Cell Type (anno2)", y = "Proportion", fill = "Cohort",
         title = "Cell Type Frequency: HC vs Stroke") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT, "09_freq_boxplot_anno2_cohort.pdf"), p, width = 12, height = 6)
  ggsave(file.path(OUT, "09_freq_boxplot_anno2_cohort.png"), p, width = 12, height = 6, dpi = 200)
  cat("  Saved: 09_freq_boxplot_anno2_cohort (fallback)\n")
})

# Layer 2 (g3)
tryCatch({
  sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))

  meta_is <- sobj_is@meta.data
  freq_g3 <- meta_is %>%
    mutate(g3 = factor(g3, levels = c("1", "2"), labels = c("Good", "Bad"))) %>%
    group_by(patient_name, g3, anno2) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient_name) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()

  p_g3 <- ggplot(freq_g3, aes(x = anno2, y = freq, fill = g3)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
               size = 0.8, alpha = 0.5) +
    stat_compare_means(aes(group = g3), method = "wilcox.test",
                       label = "p.signif", size = 3) +
    theme_bw() +
    labs(x = "Cell Type (anno2)", y = "Proportion", fill = "Outcome (g3)",
         title = "Cell Type Frequency: Good vs Bad (IS only)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT, "10_freq_boxplot_anno2_g3.pdf"), p_g3, width = 12, height = 6)
  ggsave(file.path(OUT, "10_freq_boxplot_anno2_g3.png"), p_g3, width = 12, height = 6, dpi = 200)
  cat("  Saved: 10_freq_boxplot_anno2_g3\n")

  # anno1 level
  freq_anno1_g3 <- meta_is %>%
    mutate(g3 = factor(g3, levels = c("1", "2"), labels = c("Good", "Bad"))) %>%
    group_by(patient_name, g3, anno1) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient_name) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()

  p_anno1_g3 <- ggplot(freq_anno1_g3, aes(x = anno1, y = freq, fill = g3)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    stat_compare_means(aes(group = g3), method = "wilcox.test",
                       label = "p.signif", size = 2.5) +
    theme_bw() +
    labs(x = "Cell Type (anno1)", y = "Proportion", fill = "Outcome (g3)",
         title = "Cell Type Frequency: Good vs Bad (IS, anno1)") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))

  ggsave(file.path(OUT, "12_freq_boxplot_anno1_g3.pdf"), p_anno1_g3, width = 20, height = 7)
  ggsave(file.path(OUT, "12_freq_boxplot_anno1_g3.png"), p_anno1_g3, width = 20, height = 7, dpi = 200)
  cat("  Saved: 12_freq_boxplot_anno1_g3\n")

  rm(sobj_is); gc(verbose = FALSE)
}, error = function(e) cat("  g3 freq plots error:", e$message, "\n"))

# anno1 × cohort
meta <- sobj@meta.data
freq_anno1 <- meta %>%
  group_by(patient_name, cohort, anno1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_name) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p_anno1_cohort <- ggplot(freq_anno1, aes(x = anno1, y = freq, fill = cohort)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
  stat_compare_means(aes(group = cohort), method = "wilcox.test",
                     label = "p.signif", size = 2.5) +
  theme_bw() +
  labs(x = "Cell Type (anno1)", y = "Proportion", fill = "Cohort",
       title = "Cell Type Frequency: HC vs Stroke (anno1)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))

ggsave(file.path(OUT, "11_freq_boxplot_anno1_cohort.pdf"), p_anno1_cohort, width = 20, height = 7)
ggsave(file.path(OUT, "11_freq_boxplot_anno1_cohort.png"), p_anno1_cohort, width = 20, height = 7, dpi = 200)
cat("  Saved: 11_freq_boxplot_anno1_cohort\n")

# ---- Figure 6: FGS signature plots ----
cat("\n[3/6] FGS signature analysis...\n")

fgs_genes <- read.csv(file.path(BASE, "fgs/whole_IS_50_cmgi_genes.csv"))
top30 <- head(fgs_genes$gene, 30)
top_in_data <- intersect(top30, rownames(sobj))
cat("  Top 30 FGS genes (in data):", length(top_in_data), "\n")

if (length(top_in_data) > 5) {
  # Add module score
  sobj <- AddModuleScore(sobj, features = list(top_in_data), name = "FGS_score", seed = 42)

  # Heatmap: FGS genes by anno1 (aggregated)
  tryCatch({
    p_hm <- plot_heatmap_genes(sobj, features = top_in_data,
                                group.by = "anno1", aggregate = TRUE,
                                normalize = "row")
    ggsave(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.pdf"), p_hm, width = 12, height = 10)
    ggsave(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.png"), p_hm, width = 12, height = 10, dpi = 200)
    cat("  Saved: 13_heatmap_fgs_top30_by_anno1\n")
  }, error = function(e) cat("  Heatmap error:", e$message, "\n"))

  # Heatmap split by cohort
  tryCatch({
    p_hm2 <- plot_heatmap_genes(sobj, features = top_in_data,
                                 group.by = "anno1", split.by = "cohort",
                                 aggregate = TRUE, normalize = "row")
    ggsave(file.path(OUT, "14_heatmap_fgs_top30_by_anno1_cohort.pdf"), p_hm2, width = 18, height = 10)
    ggsave(file.path(OUT, "14_heatmap_fgs_top30_by_anno1_cohort.png"), p_hm2, width = 18, height = 10, dpi = 200)
    cat("  Saved: 14_heatmap_fgs_top30_by_anno1_cohort\n")
  }, error = function(e) cat("  Heatmap split error:", e$message, "\n"))

  # FGS module score per compartment × cohort (pseudobulk boxplot)
  score_meta <- sobj@meta.data %>%
    group_by(patient_name, cohort, anno2) %>%
    summarise(FGS_score = mean(FGS_score1, na.rm = TRUE), .groups = "drop")

  p_score <- ggplot(score_meta, aes(x = anno2, y = FGS_score, fill = cohort)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
               size = 1, alpha = 0.5) +
    stat_compare_means(aes(group = cohort), method = "wilcox.test",
                       label = "p.signif", size = 3) +
    theme_bw() +
    labs(x = "Compartment", y = "FGS Module Score (pseudobulk)",
         title = "FGS Signature: HC vs Stroke") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT, "15_fgs_score_boxplot_anno2_cohort.pdf"), p_score, width = 12, height = 6)
  ggsave(file.path(OUT, "15_fgs_score_boxplot_anno2_cohort.png"), p_score, width = 12, height = 6, dpi = 200)
  cat("  Saved: 15_fgs_score_boxplot_anno2_cohort\n")

  # g3 (IS only)
  sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
  score_g3 <- sobj_is@meta.data %>%
    mutate(g3 = factor(g3, levels = c("1", "2"), labels = c("Good", "Bad"))) %>%
    group_by(patient_name, g3, anno2) %>%
    summarise(FGS_score = mean(FGS_score1, na.rm = TRUE), .groups = "drop")

  p_score_g3 <- ggplot(score_g3, aes(x = anno2, y = FGS_score, fill = g3)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
               size = 1, alpha = 0.5) +
    stat_compare_means(aes(group = g3), method = "wilcox.test",
                       label = "p.signif", size = 3) +
    theme_bw() +
    labs(x = "Compartment", y = "FGS Module Score (pseudobulk)",
         title = "FGS Signature: Good vs Bad (IS only)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT, "16_fgs_score_boxplot_anno2_g3.pdf"), p_score_g3, width = 12, height = 6)
  ggsave(file.path(OUT, "16_fgs_score_boxplot_anno2_g3.png"), p_score_g3, width = 12, height = 6, dpi = 200)
  cat("  Saved: 16_fgs_score_boxplot_anno2_g3\n")
  rm(sobj_is); gc(verbose = FALSE)
}

# ---- Figure 7: Key gene pseudobulk boxplots ----
cat("\n[4/6] Key gene pseudobulk boxplots...\n")

key_genes <- c("S100A8", "S100A9", "TXNIP", "TNF", "IL1B", "CXCL8",
               "ISG15", "IFI6", "GZMB", "FOXP3", "HLA-DRA", "CD74")
key_genes <- intersect(key_genes, rownames(sobj))

# Manual pseudobulk for key genes
pb_data <- AverageExpression(sobj, features = key_genes,
                              group.by = c("patient_name", "cohort", "anno2"),
                              return.seurat = FALSE)$RNA

# Reshape: rows = genes, cols = patient_cohort_anno2
if (!is.null(pb_data)) {
  pb_long <- as.data.frame(t(pb_data)) %>%
    tibble::rownames_to_column("group_id") %>%
    tidyr::separate(group_id, into = c("patient_name", "cohort", "anno2"),
                    sep = "_", extra = "merge") %>%
    tidyr::pivot_longer(cols = all_of(key_genes), names_to = "gene", values_to = "expr")

  p_key <- ggplot(pb_long, aes(x = anno2, y = expr, fill = cohort)) +
    geom_boxplot(outlier.size = 0.3, position = position_dodge(0.8)) +
    stat_compare_means(aes(group = cohort), method = "wilcox.test",
                       label = "p.signif", size = 2.5) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    theme_bw() +
    labs(x = "Compartment", y = "Expression (pseudobulk)",
         title = "Key Genes: HC vs Stroke") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 8))

  ggsave(file.path(OUT, "17_key_genes_boxplot_anno2_cohort.png"), p_key,
         width = 20, height = 15, dpi = 150)
  cat("  Saved: 17_key_genes_boxplot_anno2_cohort\n")
}

# g3 comparison
sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
pb_g3 <- AverageExpression(sobj_is, features = key_genes,
                            group.by = c("patient_name", "g3", "anno2"),
                            return.seurat = FALSE)$RNA

if (!is.null(pb_g3)) {
  pb_g3_long <- as.data.frame(t(pb_g3)) %>%
    tibble::rownames_to_column("group_id") %>%
    tidyr::separate(group_id, into = c("patient_name", "g3", "anno2"),
                    sep = "_", extra = "merge") %>%
    mutate(g3 = factor(g3, levels = c("1", "2"), labels = c("Good", "Bad"))) %>%
    tidyr::pivot_longer(cols = all_of(key_genes), names_to = "gene", values_to = "expr")

  p_key_g3 <- ggplot(pb_g3_long, aes(x = anno2, y = expr, fill = g3)) +
    geom_boxplot(outlier.size = 0.3, position = position_dodge(0.8)) +
    stat_compare_means(aes(group = g3), method = "wilcox.test",
                       label = "p.signif", size = 2.5) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    theme_bw() +
    labs(x = "Compartment", y = "Expression (pseudobulk)",
         title = "Key Genes: Good vs Bad (IS only)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 8))

  ggsave(file.path(OUT, "18_key_genes_boxplot_anno2_g3.png"), p_key_g3,
         width = 20, height = 15, dpi = 150)
  cat("  Saved: 18_key_genes_boxplot_anno2_g3\n")
}
rm(sobj_is); gc(verbose = FALSE)

# ---- Figure 8: MASC summary ----
cat("\n[5/6] MASC summary...\n")

masc_files <- list.files(file.path(BASE, "MASC"), pattern = "_results.qs$", full.names = TRUE)
masc_summary <- list()
for (f in masc_files) {
  name <- gsub("_results.qs$", "", basename(f))
  res <- qread(f)
  or_col <- grep("\\.OR$", colnames(res), value = TRUE)[1]
  pval_col <- "model.pvalue"
  fdr_col <- "model.pvalue.fdr"
  if (!is.null(or_col) && or_col %in% colnames(res)) {
    df <- data.frame(
      comparison = name,
      cluster = res$cluster,
      OR = res[[or_col]],
      pvalue = res[[pval_col]],
      fdr = res[[fdr_col]],
      stringsAsFactors = FALSE
    )
    masc_summary[[name]] <- df
  }
}
if (length(masc_summary) > 0) {
  masc_all <- do.call(rbind, masc_summary)
  write.csv(masc_all, file.path(OUT, "19_masc_summary_all.csv"), row.names = FALSE)
  cat("  Saved: 19_masc_summary_all.csv (", nrow(masc_all), "rows)\n")
}

# ---- Figure 9: FGS gene DotPlot by anno1 ----
cat("\n[6/6] FGS gene DotPlot...\n")

# Show where FGS top genes are expressed
fgs_top20 <- head(fgs_genes$gene, 20)
fgs_top20 <- intersect(fgs_top20, rownames(sobj))

p_fgs_dot <- DotPlot(sobj, features = fgs_top20, group.by = "anno1",
                      dot.scale = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8)) +
  coord_flip() +
  ggtitle("FGS Top 20 Genes Expression by Cell Type")

ggsave(file.path(OUT, "20_dotplot_fgs_top20.pdf"), p_fgs_dot, width = 14, height = 10)
ggsave(file.path(OUT, "20_dotplot_fgs_top20.png"), p_fgs_dot, width = 14, height = 10, dpi = 200)
cat("  Saved: 20_dotplot_fgs_top20\n")

rm(sobj); gc(verbose = FALSE)

cat("\n=== DESCRIPTIVE FIGURES PART 2 COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT, "\n")
cat("Files:\n")
for (f in sort(list.files(OUT))) cat("  ", f, "\n")
