#!/usr/bin/env Rscript
# Figure 7: Integrative Summary
# Panels A-B for publication

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

outdir <- "/data/user3/sobj/stroke_hc_v8_2/figures/v3_paper"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
basedir <- "/data/user3/sobj/stroke_hc_v8_2"

# Color scheme
color_hc <- "#2166AC"
color_stroke <- "#B2182B"

# Publication theme
theme_pub <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================================
# Panel A: Multi-method concordance heatmap
# ============================================================================
cat("\n=== Panel A: Multi-method concordance heatmap ===\n")

# --- 1. MASC ---
cat("Loading MASC...\n")
masc_res <- qread(file.path(basedir, "frequency/MASC/MASC_anno1_cohort_v2_results.qs"))
cat("MASC class:", class(masc_res), "\n")
# Extract MASC results - it's a data.frame with cluster, model.pvalue, model.or
if(is.data.frame(masc_res)) {
  masc_df <- masc_res
} else if(is.list(masc_res)) {
  # Try to find the results data frame
  if("results" %in% names(masc_res)) masc_df <- masc_res$results
  else masc_df <- masc_res[[1]]
}
cat("MASC columns:", paste(colnames(masc_df), collapse=", "), "\n")
cat("MASC rows:", nrow(masc_df), "\n")

# Parse MASC - need cluster name, p-value, and direction (OR)
# The cluster column uses dots instead of special chars
# MASC uses dot-encoded names: CD14..Monocyte = CD14+ Monocyte, Platelet.PLA = Platelet/PLA
masc_name_map <- c(
  "B_cell" = "B_cell",
  "CD14..Monocyte" = "CD14+ Monocyte",
  "CD16..Monocyte" = "CD16+ Monocyte",
  "CD4_S100A8_CD14" = "CD4_S100A8_CD14",
  "CD4_S100A8_CSF3R" = "CD4_S100A8_CSF3R",
  "CD4..T_Naive.Memory" = "CD4+ T_Naive/Memory",
  "CD8..T_Cytotoxic" = "CD8+ T_Cytotoxic",
  "CD8..Trm" = "CD8+ Trm",
  "cDC1" = "cDC1",
  "cDC2" = "cDC2",
  "Inflammatory.Monocyte" = "Inflammatory Monocyte",
  "ISG..Myeloid" = "ISG+ Myeloid",
  "ISG..T_cell" = "ISG+ T_cell",
  "MAIT" = "MAIT",
  "Mast_cell" = "Mast_cell",
  "NK_cell" = "NK_cell",
  "pDC" = "pDC",
  "Plasma_cell" = "Plasma_cell",
  "Platelet.PLA" = "Platelet/PLA",
  "Proliferating" = "Proliferating",
  "Treg" = "Treg"
)

masc_parsed <- masc_df %>%
  mutate(
    cell_type = masc_name_map[cluster],
    masc_sig = model.pvalue < 0.05,
    masc_or = cohortStroke.OR,
    masc_dir = ifelse(masc_or > 1, "Stroke", "HC"),
    masc_effect = log2(masc_or)
  ) %>%
  select(cell_type, masc_sig, masc_dir, masc_effect, masc_pval = model.pvalue, masc_or)

cat("MASC parsed cell types:", paste(masc_parsed$cell_type, collapse=", "), "\n")

# --- 2. scCODA ---
cat("Loading scCODA...\n")
sccoda <- read.csv(file.path(basedir, "frequency/sccoda/sccoda_L1_cohort_credible.csv"),
                    stringsAsFactors = FALSE)
sccoda_df <- sccoda %>%
  rename(cell_type = Cell.Type) %>%
  mutate(sccoda_sig = Final.Parameter == "True" | Final.Parameter == TRUE)

# --- 3. MILO concordance ---
cat("Loading MILO...\n")
milo_conc <- read.csv(file.path(basedir, "frequency/milo/L1_concordance_3methods.csv"),
                       stringsAsFactors = FALSE)
milo_df <- milo_conc %>%
  select(cell_type, milo_sig = milo_da, milo_dir, milo_logFC = milo_logFC)

# --- 4. propeller ---
cat("Loading propeller...\n")
propeller <- read.csv(file.path(basedir, "frequency/propeller/propeller_L1_cohort.csv"),
                       stringsAsFactors = FALSE)
propeller_df <- propeller %>%
  rename(cell_type = BaselineProp.clusters) %>%
  mutate(
    propeller_sig = FDR < 0.05,
    propeller_dir = ifelse(PropRatio > 1, "HC", "Stroke"),
    propeller_effect = log2(PropRatio)
  ) %>%
  select(cell_type, propeller_sig, propeller_dir, propeller_effect, propeller_fdr = FDR)

# --- 5. MELD ---
cat("Loading MELD...\n")
meld <- read.csv(file.path(basedir, "frequency/meld/meld_L1_summary_by_celltype.csv"),
                  stringsAsFactors = FALSE)
meld_df <- meld %>%
  rename(cell_type = anno1) %>%
  mutate(
    meld_enriched = ifelse(mean_p_stroke > 0.6, "Stroke", ifelse(mean_p_stroke < 0.4, "HC", "NS")),
    meld_sig = meld_enriched != "NS",
    meld_effect = mean_p_stroke - 0.5  # deviation from 0.5
  ) %>%
  select(cell_type, meld_sig, meld_enriched, meld_effect, meld_p_stroke = mean_p_stroke)

# --- 6. Augur ---
cat("Loading Augur...\n")
augur_l1 <- read.csv(file.path(basedir, "frequency/augur/augur_L1_cohort_auc.csv"),
                      stringsAsFactors = FALSE)
augur_df <- augur_l1 %>%
  rename(cell_type = cell_type) %>%
  mutate(
    augur_sig = auc > 0.7,
    augur_effect = auc
  ) %>%
  select(cell_type, augur_sig, augur_effect)

# --- 7. DEG consensus ---
cat("Loading DEG consensus...\n")
deg_res <- qread(file.path(basedir, "deg/consensus/L1_cohort_anno1_1/results_consensus.qs"))
# DEG consensus uses cluster_id column with clean names (e.g., "CD14+ Monocyte")
# mean_beta > 0 = up in Stroke, mean_beta < 0 = down in Stroke
deg_counts <- deg_res$consensus_deg_list %>%
  group_by(cluster_id) %>%
  summarise(
    n_deg = n(),
    n_up = sum(mean_beta > 0, na.rm = TRUE),
    n_down = sum(mean_beta < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(cell_type = cluster_id) %>%
  mutate(
    deg_sig = n_deg > 0,
    deg_dir = ifelse(n_up > n_down, "Stroke", "HC"),
    deg_effect = log10(n_deg + 1)
  ) %>%
  select(cell_type, deg_sig, deg_dir, deg_effect, n_deg, n_up, n_down)

cat("DEG cell types:", paste(deg_counts$cell_type, collapse=", "), "\n")

# --- 8. FGS ---
cat("Loading FGS...\n")
fgs_summary <- read.csv(file.path(basedir, "unbiased/summary/fgs_per_celltype_summary.csv"),
                         stringsAsFactors = FALSE)
cat("FGS columns:", paste(colnames(fgs_summary), collapse=", "), "\n")
cat("FGS rows:", nrow(fgs_summary), "\n")

# --- Merge all methods ---
cat("\nMerging all methods...\n")

# Get all unique cell types across methods
all_cells <- sort(unique(c(
  masc_parsed$cell_type,
  sccoda_df$cell_type,
  milo_df$cell_type,
  propeller_df$cell_type,
  meld_df$cell_type,
  augur_df$cell_type,
  deg_counts$cell_type
)))

cat("All cell types (", length(all_cells), "):", paste(all_cells, collapse=", "), "\n")

# Create master data frame
master <- data.frame(cell_type = all_cells, stringsAsFactors = FALSE)

# Join each method
master <- master %>%
  left_join(masc_parsed %>% select(cell_type, masc_sig, masc_dir, masc_effect), by = "cell_type") %>%
  left_join(sccoda_df %>% select(cell_type, sccoda_sig), by = "cell_type") %>%
  left_join(milo_df %>% select(cell_type, milo_sig, milo_dir, milo_logFC), by = "cell_type") %>%
  left_join(propeller_df %>% select(cell_type, propeller_sig, propeller_dir, propeller_effect), by = "cell_type") %>%
  left_join(meld_df %>% select(cell_type, meld_sig, meld_enriched, meld_effect), by = "cell_type") %>%
  left_join(augur_df %>% select(cell_type, augur_sig, augur_effect), by = "cell_type") %>%
  left_join(deg_counts %>% select(cell_type, deg_sig, deg_dir, n_deg), by = "cell_type")

# Count number of methods significant per cell type
master <- master %>%
  rowwise() %>%
  mutate(
    n_sig = sum(c(
      isTRUE(masc_sig),
      isTRUE(sccoda_sig),
      isTRUE(milo_sig),
      isTRUE(propeller_sig),
      isTRUE(meld_sig),
      isTRUE(augur_sig),
      isTRUE(deg_sig)
    ), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(n_sig))

# Determine overall direction for annotation
master <- master %>%
  mutate(
    overall_dir = case_when(
      !is.na(milo_dir) ~ milo_dir,
      !is.na(masc_dir) ~ masc_dir,
      !is.na(propeller_dir) ~ propeller_dir,
      !is.na(meld_enriched) & meld_enriched != "NS" ~ meld_enriched,
      TRUE ~ "NS"
    )
  )

cat("\nMaster data (sorted by n_sig):\n")
print(master %>% select(cell_type, n_sig, overall_dir, masc_sig, sccoda_sig, milo_sig, propeller_sig, meld_sig, augur_sig, deg_sig))

# Build heatmap matrix: rows = cell types, columns = methods
methods <- c("MASC", "scCODA", "MILO", "propeller", "MELD", "Augur", "DEG")
hm_mat <- matrix(0, nrow = nrow(master), ncol = length(methods),
                  dimnames = list(master$cell_type, methods))

# Fill with signed significance: +1 = sig Stroke-enriched, -1 = sig HC-enriched, 0 = NS
for(i in seq_len(nrow(master))) {
  ct <- master$cell_type[i]
  # MASC
  if(isTRUE(master$masc_sig[i])) hm_mat[i, "MASC"] <- ifelse(master$masc_dir[i] == "Stroke", 1, -1)
  # scCODA
  if(isTRUE(master$sccoda_sig[i])) hm_mat[i, "scCODA"] <- ifelse(master$overall_dir[i] == "Stroke", 1, -1)
  # MILO
  if(isTRUE(master$milo_sig[i])) hm_mat[i, "MILO"] <- ifelse(master$milo_dir[i] == "IS" | master$milo_dir[i] == "Stroke", 1, -1)
  # propeller
  if(isTRUE(master$propeller_sig[i])) hm_mat[i, "propeller"] <- ifelse(master$propeller_dir[i] == "Stroke", 1, -1)
  # MELD
  if(isTRUE(master$meld_sig[i])) hm_mat[i, "MELD"] <- ifelse(master$meld_enriched[i] == "Stroke", 1, -1)
  # Augur (no direction, just separability)
  if(isTRUE(master$augur_sig[i])) hm_mat[i, "Augur"] <- 1
  # DEG (direction = majority of DEGs)
  if(isTRUE(master$deg_sig[i])) hm_mat[i, "DEG"] <- ifelse(master$deg_dir[i] == "Stroke", 1, -1)
}

# Effect size matrix for continuous color
eff_mat <- matrix(NA, nrow = nrow(master), ncol = length(methods),
                   dimnames = list(master$cell_type, methods))

for(i in seq_len(nrow(master))) {
  eff_mat[i, "MASC"] <- master$masc_effect[i]
  # scCODA is binary, use a fixed value if sig
  eff_mat[i, "scCODA"] <- ifelse(isTRUE(master$sccoda_sig[i]), hm_mat[i, "scCODA"] * 1, 0)
  eff_mat[i, "MILO"] <- ifelse(!is.na(master$milo_logFC[i]), master$milo_logFC[i], 0)
  eff_mat[i, "propeller"] <- ifelse(!is.na(master$propeller_effect[i]), master$propeller_effect[i], 0)
  eff_mat[i, "MELD"] <- ifelse(!is.na(master$meld_effect[i]), master$meld_effect[i], 0)
  eff_mat[i, "Augur"] <- ifelse(!is.na(master$augur_effect[i]), master$augur_effect[i] - 0.5, 0)
  eff_mat[i, "DEG"] <- ifelse(!is.na(master$n_deg[i]), log10(master$n_deg[i] + 1) * hm_mat[i, "DEG"], 0)
}

# Clip extreme values for visualization
eff_mat[is.na(eff_mat)] <- 0

# Create the heatmap
# Use a discrete color approach: blue (HC) - white (NS) - red (Stroke)
col_fun <- colorRamp2(c(-1, 0, 1), c(color_hc, "white", color_stroke))

# Annotation: n_sig and direction
ha_right <- rowAnnotation(
  "# Methods" = anno_barplot(
    master$n_sig,
    width = unit(2.5, "cm"),
    gp = gpar(fill = ifelse(master$overall_dir == "Stroke", color_stroke,
                             ifelse(master$overall_dir == "HC", color_hc, "grey70")),
              col = NA),
    axis_param = list(gp = gpar(fontsize = 8))
  ),
  "Direction" = master$overall_dir,
  col = list("Direction" = c("HC" = color_hc, "Stroke" = color_stroke, "IS" = color_stroke, "NS" = "grey80")),
  annotation_legend_param = list("Direction" = list(title = "Enriched in")),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10)
)

# Row order by n_sig
row_order <- order(master$n_sig, decreasing = TRUE)

ht <- Heatmap(
  hm_mat[row_order, ],
  name = "Significance",
  col = col_fun,
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 30,
  right_annotation = ha_right[row_order, ],
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- hm_mat[row_order, ][i, j]
    if(val != 0) {
      grid.circle(x, y, r = unit(3, "mm"),
                  gp = gpar(fill = ifelse(val > 0, color_stroke, color_hc),
                            col = "black", lwd = 0.5, alpha = 0.8))
    } else {
      grid.circle(x, y, r = unit(3, "mm"),
                  gp = gpar(fill = "grey95", col = "grey70", lwd = 0.3))
    }
  },
  heatmap_legend_param = list(
    title = "Significant",
    at = c(-1, 0, 1),
    labels = c("HC-enriched", "NS", "Stroke-enriched"),
    legend_height = unit(3, "cm")
  ),
  column_title = "Multi-Method Concordance: Cell Type Differential Abundance",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

png(file.path(outdir, "fig7_A_multimethod_heatmap.png"), width = 12, height = 10, units = "in", res = 300)
draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()
cat("Saved: fig7_A_multimethod_heatmap.png\n")


# ============================================================================
# Panel B: Augur L1 vs L2 scatter
# ============================================================================
cat("\n=== Panel B: Augur L1 vs L2 scatter ===\n")

augur_combined <- read.csv(file.path(basedir, "frequency/augur/augur_combined_auc.csv"),
                            stringsAsFactors = FALSE)

# Split into L1 and L2
augur_l1 <- augur_combined %>% filter(layer == "L1 (HC vs Stroke)") %>% select(cell_type, auc_l1 = auc)
augur_l2 <- augur_combined %>% filter(layer == "L2 (Good vs Bad)") %>% select(cell_type, auc_l2 = auc)

augur_scatter <- augur_l1 %>%
  inner_join(augur_l2, by = "cell_type")

cat("Augur scatter data (", nrow(augur_scatter), "cell types):\n")
print(augur_scatter %>% arrange(desc(auc_l1)))

# Highlight key cell types
highlight_types <- c("CD14+ Monocyte", "CD16+ Monocyte", "Inflammatory Monocyte",
                      "cDC1", "NK_cell", "CD4_S100A8_CSF3R")

augur_scatter <- augur_scatter %>%
  mutate(
    highlight = cell_type %in% highlight_types |
                auc_l1 > 0.8 | auc_l2 > 0.6,
    is_mono = grepl("Mono", cell_type),
    cell_label = gsub("_", " ", cell_type),
    cell_label = gsub("\\+", "+", cell_label)
  )

p_augur <- ggplot(augur_scatter, aes(x = auc_l1, y = auc_l2)) +
  # Reference lines
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey70") +
  # Random performance line
  annotate("text", x = 0.52, y = 0.495, label = "Random", size = 3, color = "grey50", fontface = "italic") +
  # Points
  geom_point(aes(fill = is_mono), size = 4, shape = 21, stroke = 0.5, color = "black") +
  # Labels
  geom_text_repel(
    aes(label = cell_label),
    size = 3,
    max.overlaps = 25,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    min.segment.length = 0.1,
    seed = 42
  ) +
  scale_fill_manual(
    values = c("TRUE" = color_stroke, "FALSE" = "grey60"),
    labels = c("TRUE" = "Monocyte", "FALSE" = "Other"),
    name = "Cell Type"
  ) +
  labs(
    title = "Cell Type Separability: Augur AUC",
    subtitle = "Layer 1 (HC vs Stroke) vs Layer 2 (Good vs Bad Outcome)",
    x = "L1 AUC (HC vs Stroke)",
    y = "L2 AUC (Good vs Bad)"
  ) +
  theme_pub +
  theme(legend.position = c(0.15, 0.85)) +
  coord_cartesian(xlim = c(0.55, 0.9), ylim = c(0.47, 0.65)) +
  scale_x_continuous(breaks = seq(0.5, 0.9, 0.1)) +
  scale_y_continuous(breaks = seq(0.45, 0.65, 0.05))

png(file.path(outdir, "fig7_B_augur_scatter.png"), width = 9, height = 7, units = "in", res = 300)
print(p_augur)
dev.off()
cat("Saved: fig7_B_augur_scatter.png\n")


cat("\n=== Figure 7 complete ===\n")
cat("Output files:\n")
cat("  ", file.path(outdir, "fig7_A_multimethod_heatmap.png"), "\n")
cat("  ", file.path(outdir, "fig7_B_augur_scatter.png"), "\n")
