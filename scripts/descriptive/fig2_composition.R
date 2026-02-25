#!/usr/bin/env Rscript
# Figure 2: Cell Composition Changes
# Panels A, D, E for the stroke scRNAseq paper
# Output: /data/user3/sobj/stroke_hc_v8_2/figures/v3_paper/fig2_*.png

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
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
    plot.margin = margin(10, 10, 10, 10)
  )

cohort_cols <- c("HC" = "#2166AC", "Stroke" = "#B2182B")

# ============================================================
# Panel A: MASC Forest Plot (L1, cohort)
# ============================================================
cat("Generating Panel A: MASC forest plot...\n")

masc <- qread("/data/user3/sobj/stroke_hc_v8_2/frequency/MASC/MASC_anno1_cohort_v2_results.qs")

# Fix cluster names: CD14..Monocyte -> CD14+ Monocyte, etc.
name_map <- c(
  "B_cell"                = "B cell",
  "CD14..Monocyte"        = "CD14+ Monocyte",
  "CD16..Monocyte"        = "CD16+ Monocyte",
  "CD4_S100A8_CD14"       = "CD4 S100A8/CD14",
  "CD4_S100A8_CSF3R"      = "CD4 S100A8/CSF3R",
  "CD4..T_Naive.Memory"   = "CD4+ T Naive/Memory",
  "CD8..T_Cytotoxic"      = "CD8+ T Cytotoxic",
  "CD8..Trm"              = "CD8+ Trm",
  "cDC1"                  = "cDC1",
  "cDC2"                  = "cDC2",
  "Inflammatory.Monocyte" = "Inflammatory Monocyte",
  "ISG..Myeloid"          = "ISG+ Myeloid",
  "ISG..T_cell"           = "ISG+ T cell",
  "MAIT"                  = "MAIT",
  "Mast_cell"             = "Mast cell",
  "NK_cell"               = "NK cell",
  "pDC"                   = "pDC",
  "Plasma_cell"           = "Plasma cell",
  "Platelet.PLA"          = "Platelet/PLA",
  "Proliferating"         = "Proliferating",
  "Treg"                  = "Treg"
)

masc$cell_type <- name_map[masc$cluster]
# If any unmapped, keep original
masc$cell_type[is.na(masc$cell_type)] <- masc$cluster[is.na(masc$cell_type)]

masc$sig <- masc$model.pvalue.fdr < 0.05
masc$direction <- ifelse(masc$cohortStroke.OR > 1, "Stroke-enriched", "HC-enriched")

# Log2 OR for display
masc$log2OR <- log2(masc$cohortStroke.OR)
masc$log2OR_lower <- log2(masc$cohortStroke.OR.95pct.ci.lower)
masc$log2OR_upper <- log2(masc$cohortStroke.OR.95pct.ci.upper)

# Cap extreme values for visualization
cap_val <- 10
masc$log2OR_plot <- pmax(pmin(masc$log2OR, cap_val), -cap_val)
masc$log2OR_lower_plot <- pmax(pmin(masc$log2OR_lower, cap_val), -cap_val)
masc$log2OR_upper_plot <- pmax(pmin(masc$log2OR_upper, cap_val), -cap_val)

# Sort by log2OR
masc <- masc %>% arrange(log2OR_plot)
masc$cell_type <- factor(masc$cell_type, levels = masc$cell_type)

# Color by significance and direction
masc$fill_col <- case_when(
  !masc$sig ~ "grey60",
  masc$direction == "Stroke-enriched" ~ "#B2182B",
  TRUE ~ "#2166AC"
)

p_a <- ggplot(masc, aes(x = log2OR_plot, y = cell_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = log2OR_lower_plot, xmax = log2OR_upper_plot),
    height = 0.3, linewidth = 0.5, color = "grey30"
  ) +
  geom_point(
    aes(fill = fill_col),
    shape = 21, size = 3, stroke = 0.5
  ) +
  scale_fill_identity() +
  # Add significance stars
  geom_text(
    data = masc %>% filter(sig),
    aes(x = log2OR_upper_plot + 0.3, label = ifelse(model.pvalue.fdr < 0.001, "***",
                                                      ifelse(model.pvalue.fdr < 0.01, "**", "*"))),
    hjust = 0, size = 4, color = "black"
  ) +
  scale_x_continuous(
    limits = c(-cap_val - 0.5, cap_val + 1.5),
    breaks = seq(-10, 10, 2)
  ) +
  xlab(expression(log[2]*"(Odds Ratio: Stroke / HC)")) +
  ylab("") +
  ggtitle("MASC: Differential abundance (L1: HC vs Stroke)") +
  annotate("text", x = -cap_val + 0.5, y = 0.5, label = "HC-enriched",
           color = "#2166AC", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = cap_val - 0.5, y = 0.5, label = "Stroke-enriched",
           color = "#B2182B", fontface = "bold", size = 3.5, hjust = 1) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 13)
  )

ggsave(
  file.path(outdir, "fig2_A_masc_forest.png"),
  p_a, width = 10, height = 8, dpi = 300, bg = "white"
)
cat("  Saved fig2_A_masc_forest.png\n")

# ============================================================
# Panel D: Multi-method concordance heatmap
# ============================================================
cat("Generating Panel D: Multi-method concordance heatmap...\n")

conc <- read.csv(
  "/data/user3/sobj/stroke_hc_v8_2/frequency/milo/L1_concordance_masc_sccoda_milo.csv",
  stringsAsFactors = FALSE
)

# Prepare matrix for heatmap
# Build a tidy format: cell_type, method, value (+1 enriched in Stroke, -1 HC, 0 NS)
heatmap_data <- conc %>%
  mutate(
    MILO = case_when(
      !milo_sig ~ 0,
      milo_dir == "IS_enriched" ~ 1,
      milo_dir == "HC_enriched" ~ -1,
      TRUE ~ 0
    ),
    MASC = case_when(
      is.na(masc_sig) | !masc_sig ~ 0,
      masc_dir == "IS_enriched" ~ 1,
      masc_dir == "HC_enriched" ~ -1,
      TRUE ~ 0
    ),
    scCODA = case_when(
      !sccoda_sig ~ 0,
      # scCODA doesn't have direction in this table; infer from MILO/MASC consensus
      milo_dir == "IS_enriched" ~ 1,
      milo_dir == "HC_enriched" ~ -1,
      TRUE ~ 0
    )
  ) %>%
  select(cell_type, MILO, MASC, scCODA, n_sig, consensus) %>%
  arrange(desc(n_sig), cell_type)

# Pivot to long format
heatmap_long <- heatmap_data %>%
  pivot_longer(cols = c(MILO, MASC, scCODA), names_to = "method", values_to = "value") %>%
  mutate(
    label = case_when(
      value == 1 ~ "Stroke",
      value == -1 ~ "HC",
      TRUE ~ "NS"
    )
  )

# Order cell types by n_sig then by MILO logFC
ct_order_conc <- conc %>%
  arrange(desc(n_sig), milo_logFC) %>%
  pull(cell_type)
heatmap_long$cell_type <- factor(heatmap_long$cell_type, levels = rev(ct_order_conc))
heatmap_long$method <- factor(heatmap_long$method, levels = c("MILO", "MASC", "scCODA"))

# Add n_methods annotation
n_meth_df <- heatmap_data %>%
  select(cell_type, n_sig) %>%
  mutate(cell_type = factor(cell_type, levels = rev(ct_order_conc)))

p_d <- ggplot(heatmap_long, aes(x = method, y = cell_type, fill = label)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(
    values = c("HC" = "#2166AC", "Stroke" = "#B2182B", "NS" = "grey90"),
    name = "Direction"
  ) +
  # Add concordance count as text on the right
  geom_text(
    data = n_meth_df,
    aes(x = 4, y = cell_type, label = paste0(n_sig, "/3")),
    inherit.aes = FALSE, size = 3.5, hjust = 0.5
  ) +
  scale_x_discrete(expand = c(0, 0, 0, 1.5)) +
  annotate("text", x = 4, y = length(ct_order_conc) + 0.8,
           label = "Sig.\nmethods", size = 3, fontface = "bold") +
  ggtitle("Multi-method DA concordance (L1)") +
  xlab("") + ylab("") +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

ggsave(
  file.path(outdir, "fig2_D_concordance.png"),
  p_d, width = 7, height = 8, dpi = 300, bg = "white"
)
cat("  Saved fig2_D_concordance.png\n")

# ============================================================
# Panel E: MELD P(Stroke) by cell type
# ============================================================
cat("Generating Panel E: MELD P(Stroke) by cell type...\n")

meld <- read.csv(
  "/data/user3/sobj/stroke_hc_v8_2/frequency/meld/meld_L1_summary_by_celltype.csv",
  stringsAsFactors = FALSE
)

# Clean names
meld$cell_type <- meld$anno1

# Sort by mean P(Stroke)
meld <- meld %>% arrange(desc(mean_p_stroke))
meld$cell_type <- factor(meld$cell_type, levels = rev(meld$cell_type))

# Color by whether P(Stroke) > 0.5 (Stroke-enriched) or < 0.5 (HC-enriched)
meld$direction <- ifelse(meld$mean_p_stroke > 0.5, "Stroke-enriched", "HC-enriched")

p_e <- ggplot(meld, aes(x = mean_p_stroke, y = cell_type, fill = direction)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_segment(
    aes(x = 0.5, xend = mean_p_stroke, y = cell_type, yend = cell_type),
    color = "grey60", linewidth = 0.6
  ) +
  geom_point(shape = 21, size = 4, stroke = 0.5) +
  scale_fill_manual(
    values = c("HC-enriched" = "#2166AC", "Stroke-enriched" = "#B2182B"),
    name = ""
  ) +
  scale_x_continuous(
    limits = c(0.3, 1.0),
    breaks = seq(0.3, 1.0, 0.1),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  xlab("Mean P(Stroke) per cell") +
  ylab("") +
  ggtitle("MELD: Stroke likelihood by cell type (L1)") +
  # Add text labels for top 3 and bottom 2
  geom_text(
    data = meld %>% slice(c(1:3, (n()-1):n())),
    aes(label = sprintf("%.2f", mean_p_stroke)),
    hjust = ifelse(
      meld %>% slice(c(1:3, (n()-1):n())) %>% pull(mean_p_stroke) > 0.5,
      -0.3, 1.3
    ),
    size = 3, color = "black"
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(outdir, "fig2_E_meld.png"),
  p_e, width = 9, height = 8, dpi = 300, bg = "white"
)
cat("  Saved fig2_E_meld.png\n")

cat("Figure 2 complete!\n")
