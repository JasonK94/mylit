#!/usr/bin/env Rscript
# ==============================================================================
# Descriptive Figures for stroke_hc_v8_2
# Generates: UMAP, DotPlot, stacked bar, frequency boxplot, signature plots
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
})

# Source myR plot functions (order matters: utils first, then plots)
source("/home/user3/data_user3/git_repo/mylit/myR/R/utils_aggregation.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots_box.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/plots_heatmap.R")

cat("=== Descriptive Figures ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat("[1/8] Loading clean Seurat object...\n")
sobj <- qread(file.path(BASE, "5_strokev8_clean.qs"))
cat("  Cells:", ncol(sobj), " Genes:", nrow(sobj), "\n")
cat("  anno1 levels:", length(unique(sobj$anno1)), "\n")
cat("  anno2 levels:", length(unique(sobj$anno2)), "\n")

# Set default identity
Idents(sobj) <- "anno1"

# ---- Define palettes ----
anno2_colors <- c(
  "Tc" = "#E41A1C", "NKc" = "#FF7F00", "Mono" = "#377EB8",
  "DC" = "#4DAF4A", "Bc" = "#984EA3", "Platelet/PLA" = "#A65628",
  "Erythrocyte" = "#F781BF", "Mast_cell" = "#999999"
)

# ---- Figure 1A: UMAP by anno1, split by cohort ----
cat("\n[2/8] UMAP plots...\n")

p_umap_anno1 <- DimPlot(sobj, reduction = "umap.scvi", group.by = "anno1",
                         label = TRUE, label.size = 3, repel = TRUE,
                         pt.size = 0.1, raster = TRUE, raster.dpi = c(1024, 1024)) +
  theme(legend.position = "right", legend.text = element_text(size = 7)) +
  ggtitle("anno1 (21 cell types)")

ggsave(file.path(OUT, "01_umap_anno1.png"), p_umap_anno1,
       width = 14, height = 10, dpi = 150)
cat("  Saved: 01_umap_anno1.png\n")

p_umap_split <- DimPlot(sobj, reduction = "umap.scvi", group.by = "anno1",
                         split.by = "cohort", label = TRUE, label.size = 2.5,
                         repel = TRUE, pt.size = 0.05, raster = TRUE,
                         raster.dpi = c(1024, 1024)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6))

ggsave(file.path(OUT, "02_umap_anno1_split_cohort.png"), p_umap_split,
       width = 20, height = 10, dpi = 150)
cat("  Saved: 02_umap_anno1_split_cohort.png\n")

p_umap_anno2 <- DimPlot(sobj, reduction = "umap.scvi", group.by = "anno2",
                         label = TRUE, label.size = 4, repel = TRUE,
                         cols = anno2_colors, pt.size = 0.1,
                         raster = TRUE, raster.dpi = c(1024, 1024)) +
  ggtitle("anno2 (8 compartments)")

ggsave(file.path(OUT, "03_umap_anno2.png"), p_umap_anno2,
       width = 10, height = 8, dpi = 150)
cat("  Saved: 03_umap_anno2.png\n")

# Integration QC: color by GEM (batch)
p_umap_batch <- DimPlot(sobj, reduction = "umap.scvi", group.by = "GEM",
                         pt.size = 0.05, raster = TRUE,
                         raster.dpi = c(1024, 1024)) +
  ggtitle("Batch (GEM)")

ggsave(file.path(OUT, "04_umap_batch_GEM.png"), p_umap_batch,
       width = 10, height = 8, dpi = 150)
cat("  Saved: 04_umap_batch_GEM.png\n")

# ---- Figure 1B: DotPlot - canonical markers ----
cat("\n[3/8] DotPlot...\n")

canonical_markers <- c(
  # T cell
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B",
  # Naive/Memory
  "CCR7", "LEF1", "IL7R", "SELL",
  # Cytotoxic
  "NKG7", "GZMB", "GZMK", "PRF1", "GNLY",
  # Treg
  "FOXP3", "IL2RA", "CTLA4",
  # NK
  "NCAM1", "KLRD1", "KLRF1",
  # Monocyte
  "CD14", "LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7",
  # Inflammatory
  "IL1B", "CXCL8", "CCL2", "TNF",
  # DC
  "FCER1A", "CD1C", "CLEC9A", "IRF7", "IRF8", "LILRA4",
  # B cell
  "CD79A", "MS4A1", "CD19",
  # Plasma
  "MZB1", "JCHAIN", "IGHG1",
  # ISG
  "ISG15", "IFI6", "IFIT1", "MX1",
  # Platelet
  "PPBP", "PF4",
  # Proliferating
  "MKI67", "TOP2A"
)

# Filter to existing genes
canonical_markers <- intersect(canonical_markers, rownames(sobj))

p_dot <- DotPlot(sobj, features = canonical_markers, group.by = "anno1",
                 dot.scale = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 8)) +
  coord_flip() +
  ggtitle("Canonical Markers × anno1")

ggsave(file.path(OUT, "05_dotplot_canonical_markers.pdf"), p_dot,
       width = 14, height = 18)
ggsave(file.path(OUT, "05_dotplot_canonical_markers.png"), p_dot,
       width = 14, height = 18, dpi = 200)
cat("  Saved: 05_dotplot_canonical_markers\n")

# ---- Figure 2: Stacked bar - cell type composition ----
cat("\n[4/8] Stacked bar plots...\n")

# Proportional stacked bar by sample, sorted by Mono proportion
# cmb() from plots.R
p_cmb_sample <- cmb(sobj, identity = "anno2", group.by = "patient_name",
                     sort.by = "Mono")
ggsave(file.path(OUT, "06_stacked_bar_anno2_by_patient.pdf"), p_cmb_sample,
       width = 20, height = 6)
ggsave(file.path(OUT, "06_stacked_bar_anno2_by_patient.png"), p_cmb_sample,
       width = 20, height = 6, dpi = 200)
cat("  Saved: 06_stacked_bar_anno2_by_patient\n")

# Same but grouped by cohort with vertical line
p_cmb_cohort <- cmb(sobj, identity = "anno2", group.by = "patient_name",
                     split.by = "cohort", sort.by = "Mono")
ggsave(file.path(OUT, "07_stacked_bar_anno2_by_patient_cohort.pdf"), p_cmb_cohort,
       width = 20, height = 6)
ggsave(file.path(OUT, "07_stacked_bar_anno2_by_patient_cohort.png"), p_cmb_cohort,
       width = 20, height = 6, dpi = 200)
cat("  Saved: 07_stacked_bar_anno2_by_patient_cohort\n")

# anno1 level
p_cmb_anno1 <- cmb(sobj, identity = "anno1", group.by = "patient_name",
                    split.by = "cohort", sort.by = "CD14+ Monocyte")
ggsave(file.path(OUT, "08_stacked_bar_anno1_by_patient_cohort.pdf"), p_cmb_anno1,
       width = 20, height = 8)
ggsave(file.path(OUT, "08_stacked_bar_anno1_by_patient_cohort.png"), p_cmb_anno1,
       width = 20, height = 8, dpi = 200)
cat("  Saved: 08_stacked_bar_anno1_by_patient_cohort\n")

# ---- Figure 3: Frequency boxplots (cell type proportion per patient) ----
cat("\n[5/8] Frequency boxplots...\n")

# pb_freq from plots_box.R: cell type frequency per sample
p_freq_anno2 <- pb_freq(sobj, cluster = "anno2", group = "cohort",
                         sample = "patient_name")
ggsave(file.path(OUT, "09_freq_boxplot_anno2_cohort.pdf"), p_freq_anno2,
       width = 12, height = 6)
ggsave(file.path(OUT, "09_freq_boxplot_anno2_cohort.png"), p_freq_anno2,
       width = 12, height = 6, dpi = 200)
cat("  Saved: 09_freq_boxplot_anno2_cohort\n")

# For Layer 2 (g3), subset to IS patients only
sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
p_freq_anno2_g3 <- pb_freq(sobj_is, cluster = "anno2", group = "g3",
                            sample = "patient_name")
ggsave(file.path(OUT, "10_freq_boxplot_anno2_g3.pdf"), p_freq_anno2_g3,
       width = 12, height = 6)
ggsave(file.path(OUT, "10_freq_boxplot_anno2_g3.png"), p_freq_anno2_g3,
       width = 12, height = 6, dpi = 200)
cat("  Saved: 10_freq_boxplot_anno2_g3\n")

# anno1 level
p_freq_anno1 <- pb_freq(sobj, cluster = "anno1", group = "cohort",
                          sample = "patient_name")
ggsave(file.path(OUT, "11_freq_boxplot_anno1_cohort.pdf"), p_freq_anno1,
       width = 20, height = 8)
ggsave(file.path(OUT, "11_freq_boxplot_anno1_cohort.png"), p_freq_anno1,
       width = 20, height = 8, dpi = 200)
cat("  Saved: 11_freq_boxplot_anno1_cohort\n")

p_freq_anno1_g3 <- pb_freq(sobj_is, cluster = "anno1", group = "g3",
                             sample = "patient_name")
ggsave(file.path(OUT, "12_freq_boxplot_anno1_g3.pdf"), p_freq_anno1_g3,
       width = 20, height = 8)
ggsave(file.path(OUT, "12_freq_boxplot_anno1_g3.png"), p_freq_anno1_g3,
       width = 20, height = 8, dpi = 200)
cat("  Saved: 12_freq_boxplot_anno1_g3\n")
rm(sobj_is); gc(verbose = FALSE)

# ---- Figure 4: FGS signature gene expression ----
cat("\n[6/8] FGS signature analysis plots...\n")

# Load FGS genes
fgs_genes <- read.csv(file.path(BASE, "fgs/whole_IS_50_cmgi_genes.csv"))
top30_genes <- head(fgs_genes$gene, 30)
cat("  Top 30 FGS genes:", paste(head(top30_genes, 10), collapse = ", "), "...\n")

# Add module score for FGS top genes
top_genes_in_data <- intersect(top30_genes, rownames(sobj))
if (length(top_genes_in_data) > 5) {
  sobj <- AddModuleScore(sobj, features = list(top_genes_in_data),
                         name = "FGS_score", seed = 42)

  # Heatmap of FGS genes by anno1
  p_hm <- plot_heatmap_genes(sobj, features = top_genes_in_data,
                              group.by = "anno1", aggregate = TRUE,
                              normalize = "row")
  ggsave(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.pdf"), p_hm,
         width = 12, height = 10)
  ggsave(file.path(OUT, "13_heatmap_fgs_top30_by_anno1.png"), p_hm,
         width = 12, height = 10, dpi = 200)
  cat("  Saved: 13_heatmap_fgs_top30_by_anno1\n")

  # Heatmap split by cohort
  p_hm_split <- plot_heatmap_genes(sobj, features = top_genes_in_data,
                                    group.by = "anno1", split.by = "cohort",
                                    aggregate = TRUE, normalize = "row")
  ggsave(file.path(OUT, "14_heatmap_fgs_top30_by_anno1_cohort.pdf"), p_hm_split,
         width = 18, height = 10)
  ggsave(file.path(OUT, "14_heatmap_fgs_top30_by_anno1_cohort.png"), p_hm_split,
         width = 18, height = 10, dpi = 200)
  cat("  Saved: 14_heatmap_fgs_top30_by_anno1_cohort\n")

  # FGS module score boxplot by compartment × cohort
  p_score_cohort <- mybox_pb(sobj, features = "FGS_score1",
                              sample_col = "patient_name",
                              group.by = "anno2", split.by = "cohort",
                              violin = TRUE, add_stats = TRUE)
  ggsave(file.path(OUT, "15_fgs_score_boxplot_anno2_cohort.pdf"), p_score_cohort,
         width = 14, height = 6)
  ggsave(file.path(OUT, "15_fgs_score_boxplot_anno2_cohort.png"), p_score_cohort,
         width = 14, height = 6, dpi = 200)
  cat("  Saved: 15_fgs_score_boxplot_anno2_cohort\n")

  # FGS module score by g3 (IS only)
  sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
  p_score_g3 <- mybox_pb(sobj_is, features = "FGS_score1",
                          sample_col = "patient_name",
                          group.by = "anno2", split.by = "g3",
                          violin = TRUE, add_stats = TRUE)
  ggsave(file.path(OUT, "16_fgs_score_boxplot_anno2_g3.pdf"), p_score_g3,
         width = 14, height = 6)
  ggsave(file.path(OUT, "16_fgs_score_boxplot_anno2_g3.png"), p_score_g3,
         width = 14, height = 6, dpi = 200)
  cat("  Saved: 16_fgs_score_boxplot_anno2_g3\n")
  rm(sobj_is); gc(verbose = FALSE)
}

# ---- Figure 5: Key gene expression boxplots ----
cat("\n[7/8] Key gene boxplots (pseudobulk)...\n")

key_genes <- c("S100A8", "S100A9", "TXNIP", "TNF", "IL1B", "CXCL8",
               "ISG15", "IFI6", "GZMB", "FOXP3", "HLA-DRA", "CD74")
key_genes <- intersect(key_genes, rownames(sobj))

# Pseudobulk boxplot by anno2 × cohort
p_key_cohort <- mybox_pb(sobj, features = key_genes,
                          sample_col = "patient_name",
                          group.by = "anno2", split.by = "cohort",
                          violin = TRUE, add_stats = TRUE, ncol = 4)
ggsave(file.path(OUT, "17_key_genes_boxplot_anno2_cohort.png"), p_key_cohort,
       width = 24, height = 18, dpi = 150)
cat("  Saved: 17_key_genes_boxplot_anno2_cohort.png\n")

# Same for g3
sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
p_key_g3 <- mybox_pb(sobj_is, features = key_genes,
                      sample_col = "patient_name",
                      group.by = "anno2", split.by = "g3",
                      violin = TRUE, add_stats = TRUE, ncol = 4)
ggsave(file.path(OUT, "18_key_genes_boxplot_anno2_g3.png"), p_key_g3,
       width = 24, height = 18, dpi = 150)
cat("  Saved: 18_key_genes_boxplot_anno2_g3.png\n")
rm(sobj_is); gc(verbose = FALSE)

# ---- Figure 6: MASC summary ----
cat("\n[8/8] Consolidating MASC results...\n")

masc_files <- list.files(file.path(BASE, "MASC"), pattern = "_results.qs$", full.names = TRUE)
masc_summary <- list()
for (f in masc_files) {
  name <- gsub("_results.qs$", "", basename(f))
  res <- qread(f)
  # Extract key columns
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
  cat("  Saved: 19_masc_summary_all.csv (", nrow(masc_all), "rows )\n")
}

rm(sobj); gc(verbose = FALSE)

cat("\n=== DESCRIPTIVE FIGURES COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT, "\n")
cat("Files:\n")
cat(paste(" ", list.files(OUT), collapse = "\n"), "\n")
