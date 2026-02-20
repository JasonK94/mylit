#!/usr/bin/env Rscript
# Quick script to generate FGS heatmaps (fix for missing .prepare_plot_data)

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
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

# Source ALL myR utilities in correct order
source("/home/user3/data_user3/git_repo/mylit/myR/R/utils_data.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/utils_aggregation.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots_heatmap.R")

cat("=== FGS Heatmaps ===\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "figures")

# Load data
cat("Loading data...\n")
sobj <- qread(file.path(BASE, "5_strokev8_clean.qs"))
Idents(sobj) <- "anno1"

# FGS genes
fgs_genes <- read.csv(file.path(BASE, "fgs/whole_IS_50_cmgi_genes.csv"))
top30 <- head(fgs_genes$gene, 30)
top_in_data <- intersect(top30, rownames(sobj))

# Try plot_heatmap_genes
tryCatch({
  p <- plot_heatmap_genes(sobj, features = top_in_data,
                           group.by = "anno1", aggregate = TRUE,
                           normalize = "row")
  ggsave(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.pdf"), p, width = 12, height = 10)
  ggsave(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.png"), p, width = 12, height = 10, dpi = 200)
  cat("  Saved: 13_heatmap_fgs_top30_by_anno1\n")
}, error = function(e) {
  cat("  plot_heatmap_genes error:", e$message, "\n")
  cat("  Using ComplexHeatmap fallback...\n")

  # Manual heatmap with ComplexHeatmap
  avg <- AverageExpression(sobj, features = top_in_data, group.by = "anno1",
                            return.seurat = FALSE)$RNA
  mat <- as.matrix(avg)
  mat_z <- t(scale(t(mat)))  # row z-score

  pdf(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.pdf"), width = 12, height = 10)
  ht <- Heatmap(mat_z, name = "z-score",
                column_title = "FGS Top 30 Genes × anno1",
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = TRUE, show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
  draw(ht)
  dev.off()

  png(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.png"), width = 1200, height = 1000, res = 150)
  draw(ht)
  dev.off()
  cat("  Saved: 13_heatmap_fgs_top30_by_anno1 (ComplexHeatmap)\n")
})

# Heatmap split by cohort
tryCatch({
  avg_hc <- AverageExpression(subset(sobj, cohort == "HC"),
                               features = top_in_data, group.by = "anno1",
                               return.seurat = FALSE)$RNA
  avg_st <- AverageExpression(subset(sobj, cohort == "Stroke"),
                               features = top_in_data, group.by = "anno1",
                               return.seurat = FALSE)$RNA

  mat_hc <- t(scale(t(as.matrix(avg_hc))))
  mat_st <- t(scale(t(as.matrix(avg_st))))

  # Combine with annotation
  shared_cols <- intersect(colnames(mat_hc), colnames(mat_st))
  mat_combined <- cbind(mat_hc[, shared_cols], mat_st[, shared_cols])
  colnames(mat_combined) <- c(paste0("HC_", shared_cols), paste0("Stroke_", shared_cols))

  cohort_anno <- HeatmapAnnotation(
    Cohort = rep(c("HC", "Stroke"), each = length(shared_cols)),
    col = list(Cohort = c("HC" = "#2166AC", "Stroke" = "#B2182B"))
  )

  pdf(file.path(OUT, "14_heatmap_fgs_top30_by_anno1_cohort.pdf"), width = 18, height = 10)
  ht2 <- Heatmap(mat_combined, name = "z-score",
                  top_annotation = cohort_anno,
                  column_title = "FGS Top 30: HC vs Stroke by anno1",
                  column_split = rep(c("HC", "Stroke"), each = length(shared_cols)),
                  cluster_rows = TRUE, cluster_columns = FALSE,
                  show_row_names = TRUE, show_column_names = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 7),
                  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
  draw(ht2)
  dev.off()

  png(file.path(OUT, "14_heatmap_fgs_top30_by_anno1_cohort.png"), width = 1800, height = 1000, res = 150)
  draw(ht2)
  dev.off()
  cat("  Saved: 14_heatmap_fgs_top30_by_anno1_cohort\n")
}, error = function(e) cat("  Split heatmap error:", e$message, "\n"))

# Patient-level pseudobulk heatmap (top 20 FGS genes)
cat("\n  Patient-level heatmap...\n")
tryCatch({
  top20 <- head(top_in_data, 20)
  avg_patient <- AverageExpression(sobj, features = top20,
                                    group.by = c("patient_name"),
                                    return.seurat = FALSE)$RNA
  mat_patient <- t(scale(t(as.matrix(avg_patient))))

  # Get patient-level cohort info
  patient_meta <- sobj@meta.data %>%
    group_by(patient_name) %>%
    summarise(cohort = first(cohort), g3 = first(g3), .groups = "drop") %>%
    as.data.frame()
  rownames(patient_meta) <- patient_meta$patient_name

  col_order <- colnames(mat_patient)
  patient_cohort <- patient_meta[col_order, "cohort"]

  top_anno <- HeatmapAnnotation(
    Cohort = patient_cohort,
    col = list(Cohort = c("HC" = "#2166AC", "Stroke" = "#B2182B")),
    show_legend = TRUE
  )

  pdf(file.path(OUT, "21_heatmap_fgs_top20_patient_level.pdf"), width = 20, height = 8)
  ht3 <- Heatmap(mat_patient, name = "z-score",
                  top_annotation = top_anno,
                  column_title = "FGS Top 20: Patient-level Expression",
                  column_split = patient_cohort,
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  show_row_names = TRUE, show_column_names = FALSE,
                  row_names_gp = gpar(fontsize = 9),
                  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
  draw(ht3)
  dev.off()

  png(file.path(OUT, "21_heatmap_fgs_top20_patient_level.png"), width = 2000, height = 800, res = 150)
  draw(ht3)
  dev.off()
  cat("  Saved: 21_heatmap_fgs_top20_patient_level\n")
}, error = function(e) cat("  Patient heatmap error:", e$message, "\n"))

rm(sobj); gc(verbose = FALSE)
cat("\n=== HEATMAPS COMPLETE ===\n")
