#!/usr/bin/env Rscript
# Quick fix: regenerate Fig1 Panel E with correct stack order
# Mono at bottom, Tc at top

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
  library(patchwork)
})

outdir <- "/data/user3/sobj/stroke_hc_v8_2/figures/v3_paper"

theme_pub <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

cohort_cols <- c("HC" = "#2166AC", "Stroke" = "#B2182B")

cat("Loading Seurat object...\n")
sobj <- qread("/data/user3/sobj/stroke_hc_v8_2/5_strokev8_clean.qs")
cat("Loaded:", ncol(sobj), "cells\n")

md <- sobj@meta.data %>%
  select(patient_name, cohort, anno2) %>%
  filter(!is.na(anno2))

patient_prop <- md %>%
  group_by(patient_name, cohort, anno2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

mono_frac <- md %>%
  mutate(is_mono = anno2 == "Mono") %>%
  group_by(patient_name, cohort) %>%
  summarise(mono_frac = mean(is_mono), .groups = "drop")

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

patient_order <- mono_frac %>%
  arrange(cohort, sort_score) %>%
  pull(patient_name)

patient_prop$patient_name <- factor(patient_prop$patient_name, levels = patient_order)

# anno2 order: Mono at bottom (first level), Tc at top (last level)
anno2_order <- c("Mono", "NKc", "Bc", "DC", "Platelet/PLA", "Mast_cell", "Tc")
anno2_existing <- anno2_order[anno2_order %in% unique(patient_prop$anno2)]
anno2_extra <- setdiff(unique(patient_prop$anno2), anno2_existing)
anno2_levels <- c(anno2_existing, anno2_extra)
patient_prop$anno2 <- factor(patient_prop$anno2, levels = anno2_levels)

anno2_cols <- c(
  "Tc"           = "#E41A1C",
  "NKc"          = "#F781BF",
  "Bc"           = "#E5C494",
  "Mono"         = "#66C2A5",
  "DC"           = "#7570B3",
  "Platelet/PLA" = "#66A61E",
  "Mast_cell"    = "#E7298A"
)

cohort_bar <- mono_frac %>%
  select(patient_name, cohort) %>%
  mutate(patient_name = factor(patient_name, levels = patient_order))

# Use position_stack(reverse = TRUE) so first factor level (Mono) is at bottom
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

p_cohort <- ggplot(cohort_bar, aes(x = patient_name, y = 1, fill = cohort)) +
  geom_tile() +
  scale_fill_manual(values = cohort_cols, name = "Cohort") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

p_e <- p_cohort / p_bar + plot_layout(heights = c(1, 20))

ggsave(
  file.path(outdir, "fig1_E_stacked_bar.png"),
  p_e, width = 14, height = 6, dpi = 300, bg = "white"
)
cat("  Saved fig1_E_stacked_bar.png (fixed stack order)\n")
