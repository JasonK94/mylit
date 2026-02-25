#!/usr/bin/env Rscript
# Figure 1: Study Design & Immune Atlas Overview
# Panels B-E for the stroke scRNAseq paper
# Output: /data/user3/sobj/stroke_hc_v8_2/figures/v3_paper/fig1_*.png

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
  library(patchwork)
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
    plot.margin = margin(5, 5, 5, 5)
  )

# Colors
cohort_cols <- c("HC" = "#2166AC", "Stroke" = "#B2182B")

# 21 cell type colors (distinct palette)
anno1_cols <- c(
  "CD4+ T_Naive/Memory" = "#E41A1C",
  "CD8+ T_Cytotoxic"    = "#377EB8",
  "CD8+ Trm"            = "#4DAF4A",
  "MAIT"                = "#984EA3",
  "Treg"                = "#FF7F00",
  "ISG+ T_cell"         = "#FFFF33",
  "Proliferating"       = "#A65628",
  "NK_cell"             = "#F781BF",
  "CD14+ Monocyte"      = "#66C2A5",
  "CD16+ Monocyte"      = "#FC8D62",
  "Inflammatory Monocyte" = "#8DA0CB",
  "ISG+ Myeloid"        = "#E78AC3",
  "CD4_S100A8_CD14"     = "#A6D854",
  "CD4_S100A8_CSF3R"    = "#FFD92F",
  "B_cell"              = "#E5C494",
  "Plasma_cell"         = "#B3B3B3",
  "cDC1"                = "#1B9E77",
  "cDC2"                = "#D95F02",
  "pDC"                 = "#7570B3",
  "Mast_cell"           = "#E7298A",
  "Platelet/PLA"        = "#66A61E"
)

cat("Loading Seurat object...\n")
sobj <- qread("/data/user3/sobj/stroke_hc_v8_2/5_strokev8_clean.qs")
cat("Loaded:", ncol(sobj), "cells\n")

# ============================================================
# Panel B: UMAP by anno1 (21 cell types)
# ============================================================
cat("Generating Panel B: UMAP by anno1...\n")

p_b <- DimPlot(
  sobj,
  group.by = "anno1",
  reduction = "umap.scvi",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE,
  cols = anno1_cols,
  pt.size = 0.05
) +
  NoLegend() +
  ggtitle("Immune cell atlas (n = 205,277)") +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_pub +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(size = 0.5)
  )

ggsave(
  file.path(outdir, "fig1_B_umap_anno1.png"),
  p_b, width = 10, height = 8, dpi = 300, bg = "white"
)
cat("  Saved fig1_B_umap_anno1.png\n")

# ============================================================
# Panel C: UMAP split by cohort
# ============================================================
cat("Generating Panel C: UMAP split by cohort...\n")

p_c <- DimPlot(
  sobj,
  group.by = "anno1",
  split.by = "cohort",
  reduction = "umap.scvi",
  raster = FALSE,
  cols = anno1_cols,
  pt.size = 0.02,
  ncol = 2
) +
  NoLegend() +
  ggtitle("") +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_pub +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.line = element_line(size = 0.5)
  )

ggsave(
  file.path(outdir, "fig1_C_umap_split.png"),
  p_c, width = 16, height = 7, dpi = 300, bg = "white"
)
cat("  Saved fig1_C_umap_split.png\n")

# ============================================================
# Panel D: DotPlot of canonical markers x anno1
# ============================================================
cat("Generating Panel D: DotPlot canonical markers...\n")

markers <- c(
  # T cell
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B",
  # NK
  "NCAM1", "KLRF1", "GNLY", "NKG7",
  # B cell
  "CD19", "MS4A1", "CD79A",
  # Monocyte
  "CD14", "FCGR3A", "S100A8", "S100A9", "VCAN",
  # DC
  "CLEC9A", "CD1C", "LILRA4", "IRF7",
  # Platelet
  "PPBP", "PF4", "GP9",
  # Other
  "MKI67", "ISG15", "IFI44L"
)

# Order cell types logically
ct_order <- c(
  "CD4+ T_Naive/Memory", "Treg", "CD8+ T_Cytotoxic", "CD8+ Trm",
  "MAIT", "ISG+ T_cell", "Proliferating",
  "NK_cell",
  "CD14+ Monocyte", "CD16+ Monocyte", "Inflammatory Monocyte",
  "CD4_S100A8_CD14", "CD4_S100A8_CSF3R",
  "ISG+ Myeloid",
  "B_cell", "Plasma_cell",
  "cDC1", "cDC2", "pDC",
  "Mast_cell", "Platelet/PLA"
)

sobj$anno1 <- factor(sobj$anno1, levels = rev(ct_order))

p_d <- DotPlot(
  sobj,
  features = markers,
  group.by = "anno1",
  cols = c("lightgrey", "#B2182B"),
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle("Canonical marker expression") +
  xlab("") + ylab("") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(5, 10, 5, 5)
  )

ggsave(
  file.path(outdir, "fig1_D_dotplot.png"),
  p_d, width = 14, height = 8, dpi = 300, bg = "white"
)
cat("  Saved fig1_D_dotplot.png\n")

# ============================================================
# Panel E: Patient-level stacked bar (anno2, sorted by FGS score)
# ============================================================
cat("Generating Panel E: Patient-level stacked bar...\n")

# Reset anno1 factor for downstream
sobj$anno1 <- factor(sobj$anno1, levels = ct_order)

# Calculate per-patient anno2 proportions
md <- sobj@meta.data %>%
  select(patient_name, cohort, anno2) %>%
  filter(!is.na(anno2))

# Get patient-level proportions
patient_prop <- md %>%
  group_by(patient_name, cohort, anno2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Calculate per-patient CD14+ Mono fraction for sorting
mono_frac <- md %>%
  mutate(is_mono = anno2 == "Mono") %>%
  group_by(patient_name, cohort) %>%
  summarise(mono_frac = mean(is_mono), .groups = "drop")

# Try to get FGS meta score per patient if available
# FGS sobj only has IS patients, so use mono fraction as fallback for HC
fgs_sobj <- tryCatch(
  qread("/data/user3/sobj/stroke_hc_v8_2/unbiased/fgs/whole_IS_50_sobj.qs"),
  error = function(e) NULL
)

if (!is.null(fgs_sobj) && "fgs_meta_score" %in% colnames(fgs_sobj@meta.data)) {
  cat("  Using FGS meta score for IS patient sorting\n")
  fgs_patient <- fgs_sobj@meta.data %>%
    group_by(patient_name) %>%
    summarise(fgs_score = mean(fgs_meta_score, na.rm = TRUE), .groups = "drop")

  mono_frac <- mono_frac %>%
    left_join(fgs_patient, by = "patient_name") %>%
    mutate(sort_score = ifelse(!is.na(fgs_score), fgs_score, mono_frac))
  rm(fgs_sobj); gc()
} else {
  mono_frac$sort_score <- mono_frac$mono_frac
}

# Sort patients: HC by mono_frac, Stroke by FGS score
patient_order <- mono_frac %>%
  arrange(cohort, sort_score) %>%
  pull(patient_name)

patient_prop$patient_name <- factor(patient_prop$patient_name, levels = patient_order)

# anno2 order: Monocyte at bottom, T_cell at top
anno2_order <- c("Mono", "NKc", "Bc", "DC", "Platelet/PLA", "Mast_cell", "Tc")
# Keep only levels that exist
anno2_existing <- anno2_order[anno2_order %in% unique(patient_prop$anno2)]
# Add any missing
anno2_extra <- setdiff(unique(patient_prop$anno2), anno2_existing)
anno2_levels <- c(anno2_existing, anno2_extra)
patient_prop$anno2 <- factor(patient_prop$anno2, levels = anno2_levels)

# anno2 colors
anno2_cols <- c(
  "Tc"           = "#E41A1C",
  "NKc"          = "#F781BF",
  "Bc"           = "#E5C494",
  "Mono"         = "#66C2A5",
  "DC"           = "#7570B3",
  "Platelet/PLA" = "#66A61E",
  "Mast_cell"    = "#E7298A"
)

# Cohort color bar data
cohort_bar <- mono_frac %>%
  select(patient_name, cohort) %>%
  mutate(patient_name = factor(patient_name, levels = patient_order))

# Main stacked bar
p_bar <- ggplot(patient_prop, aes(x = patient_name, y = prop, fill = anno2)) +
  geom_bar(stat = "identity", width = 1, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = anno2_cols, name = "Compartment") +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  xlab("Patients (n = 100)") + ylab("Cell proportion") +
  theme_pub +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.spacing = unit(0, "lines")
  )

# Cohort color bar at top
p_cohort <- ggplot(cohort_bar, aes(x = patient_name, y = 1, fill = cohort)) +
  geom_tile() +
  scale_fill_manual(values = cohort_cols, name = "Cohort") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

# Combine with patchwork
p_e <- p_cohort / p_bar + plot_layout(heights = c(1, 20))

ggsave(
  file.path(outdir, "fig1_E_stacked_bar.png"),
  p_e, width = 14, height = 6, dpi = 300, bg = "white"
)
cat("  Saved fig1_E_stacked_bar.png\n")

cat("Figure 1 complete!\n")
