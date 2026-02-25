#!/usr/bin/env Rscript
# Figure 6: Monocyte Deep-Dive
# Panels A-D for publication

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

outdir <- "/data/user3/sobj/stroke_hc_v8_2/figures/v3_paper"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("Loading Monocyte annotated object...\n")
mono <- qread("/data/user3/sobj/stroke_hc_v8_2/subclustering/Monocyte/Monocyte_annotated.qs")
cat("Loaded:", ncol(mono), "cells\n")

# Identify biological (non-contamination) clusters
contam_clusters <- c("Platelet contamination", "T/NK contamination")
bio_clusters <- setdiff(unique(mono$subclust_identity), contam_clusters)
cat("Biological clusters (n=", length(bio_clusters), "):", paste(bio_clusters, collapse=", "), "\n")

# Subset to biological clusters for all panels
mono_bio <- subset(mono, subclust_identity %in% bio_clusters)
cat("Biological cells:", ncol(mono_bio), "\n")

# ============================================================================
# Panel A: Monocyte subcluster UMAP (annotated)
# ============================================================================
cat("\n=== Panel A: UMAP ===\n")

# Create abbreviated labels for cleaner plotting
label_map <- c(
  "Classical Mono (activated)" = "Classical\n(activated)",
  "Classical Mono (lipid-handling)" = "Classical\n(lipid-handling)",
  "Classical Mono (quiescent)" = "Classical\n(quiescent)",
  "Classical Mono (adhesion)" = "Classical\n(adhesion)",
  "S100A-high Mono (immature)" = "S100A-high\n(immature)",
  "Inflammatory Mono (IL1B+)" = "Inflammatory\n(IL1B+)",
  "Inflammatory Mono (chemokine-high)" = "Inflammatory\n(chemokine-high)",
  "Non-classical Mono (CDKN1C-hi)" = "Non-classical\n(CDKN1C-hi)",
  "Non-classical Mono (C1Q+)" = "Non-classical\n(C1Q+)",
  "Non-classical Mono (patrolling)" = "Non-classical\n(patrolling)"
)

# Palette for 10 monocyte subtypes
mono_colors <- c(
  "Classical Mono (activated)" = "#E41A1C",
  "Classical Mono (lipid-handling)" = "#FF7F00",
  "Classical Mono (quiescent)" = "#FDBF6F",
  "Classical Mono (adhesion)" = "#FB9A99",
  "S100A-high Mono (immature)" = "#A6CEE3",
  "Inflammatory Mono (IL1B+)" = "#B2DF8A",
  "Inflammatory Mono (chemokine-high)" = "#33A02C",
  "Non-classical Mono (CDKN1C-hi)" = "#6A3D9A",
  "Non-classical Mono (C1Q+)" = "#CAB2D6",
  "Non-classical Mono (patrolling)" = "#1F78B4"
)

# Get UMAP coordinates
umap_coords <- as.data.frame(Embeddings(mono_bio, "umap.scvi"))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$identity <- mono_bio$subclust_identity

# Compute cluster centroids for labels
centroids <- umap_coords %>%
  group_by(identity) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop") %>%
  mutate(label = label_map[identity])

# Randomize point order to avoid overplotting bias
set.seed(42)
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

p_umap <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = identity)) +
  geom_point(size = 0.1, alpha = 0.4, stroke = 0) +
  geom_label_repel(
    data = centroids,
    aes(label = label),
    size = 2.8,
    fontface = "bold",
    box.padding = 0.8,
    point.padding = 0.3,
    max.overlaps = 20,
    min.segment.length = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    fill = alpha("white", 0.85),
    label.size = 0.2,
    seed = 42
  ) +
  scale_color_manual(values = mono_colors, guide = "none") +
  labs(
    title = "Monocyte Subclusters",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_pub +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  coord_fixed()

png(file.path(outdir, "fig6_A_mono_umap.png"), width = 8, height = 7, units = "in", res = 300)
print(p_umap)
dev.off()
cat("Saved: fig6_A_mono_umap.png\n")


# ============================================================================
# Panel B: Monocyte marker DotPlot
# ============================================================================
cat("\n=== Panel B: DotPlot ===\n")

markers <- list(
  "Classical" = c("CD14", "S100A8", "S100A9", "VCAN", "LYZ"),
  "Non-classical" = c("FCGR3A", "CDKN1C", "CX3CR1", "MS4A7"),
  "Inflammatory" = c("IL1B", "CCL3", "CCL4", "CXCL8", "NLRP3", "TNF"),
  "Immature" = c("S100A12", "RETN", "CLU", "RNASE2"),
  "Special" = c("ABCG1", "OLR1", "C1QA", "ISG15", "THBS1", "PTGS2")
)

# Check which genes are present
all_markers <- unlist(markers)
present <- all_markers[all_markers %in% rownames(mono_bio)]
missing <- setdiff(all_markers, present)
if(length(missing) > 0) cat("Missing genes:", paste(missing, collapse=", "), "\n")

# Order identities logically
ident_order <- c(
  "Classical Mono (activated)",
  "Classical Mono (lipid-handling)",
  "Classical Mono (quiescent)",
  "Classical Mono (adhesion)",
  "S100A-high Mono (immature)",
  "Inflammatory Mono (IL1B+)",
  "Inflammatory Mono (chemokine-high)",
  "Non-classical Mono (CDKN1C-hi)",
  "Non-classical Mono (C1Q+)",
  "Non-classical Mono (patrolling)"
)

mono_bio$subclust_identity <- factor(mono_bio$subclust_identity, levels = rev(ident_order))
Idents(mono_bio) <- "subclust_identity"

p_dot <- DotPlot(mono_bio, features = present, cols = c("lightgrey", "#B2182B"),
                 dot.scale = 5, col.min = -1, col.max = 2) +
  RotatedAxis() +
  labs(title = "Monocyte Subtype Markers", x = "", y = "") +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  )

# Add group separators via vertical lines
# Find positions of group boundaries
group_sizes <- sapply(markers, function(x) sum(x %in% rownames(mono_bio)))
group_cumsum <- cumsum(group_sizes)
vlines <- group_cumsum[-length(group_cumsum)] + 0.5

p_dot <- p_dot +
  geom_vline(xintercept = vlines, linetype = "dashed", color = "grey60", linewidth = 0.3)

png(file.path(outdir, "fig6_B_mono_dotplot.png"), width = 14, height = 6, units = "in", res = 300)
print(p_dot)
dev.off()
cat("Saved: fig6_B_mono_dotplot.png\n")


# ============================================================================
# Panel C: Subcluster frequency HC vs Stroke
# ============================================================================
cat("\n=== Panel C: Frequency ===\n")

# Compute per-patient proportions from the object
mono_bio_meta <- mono_bio@meta.data
mono_bio_meta$subclust_identity <- as.character(mono_bio_meta$subclust_identity)

# Patient-level proportions
patient_freq <- mono_bio_meta %>%
  group_by(patient_name, cohort) %>%
  mutate(patient_total = n()) %>%
  group_by(patient_name, cohort, subclust_identity, patient_total) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  mutate(proportion = n_cells / patient_total * 100)

# Order by median proportion difference
median_diff <- patient_freq %>%
  group_by(subclust_identity, cohort) %>%
  summarise(med = median(proportion), .groups = "drop") %>%
  pivot_wider(names_from = cohort, values_from = med, values_fill = 0) %>%
  mutate(diff = Stroke - HC) %>%
  arrange(diff)

patient_freq$subclust_identity <- factor(
  patient_freq$subclust_identity,
  levels = median_diff$subclust_identity
)

# Wilcoxon test for each subcluster
pvals <- patient_freq %>%
  group_by(subclust_identity) %>%
  summarise(
    pval = tryCatch(
      wilcox.test(proportion[cohort == "HC"], proportion[cohort == "Stroke"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    stars = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Y position for stars
y_positions <- patient_freq %>%
  group_by(subclust_identity) %>%
  summarise(ymax = max(proportion, na.rm = TRUE), .groups = "drop")

pvals <- pvals %>% left_join(y_positions, by = "subclust_identity")

p_freq <- ggplot(patient_freq, aes(x = subclust_identity, y = proportion, fill = cohort)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8, position = position_dodge(0.75), width = 0.6) +
  geom_text(
    data = pvals,
    aes(x = subclust_identity, y = ymax * 1.1, label = stars),
    inherit.aes = FALSE,
    size = 4, vjust = 0
  ) +
  scale_fill_manual(values = c("HC" = color_hc, "Stroke" = color_stroke), name = "Cohort") +
  labs(
    title = "Monocyte Subcluster Proportions",
    x = "",
    y = "% of Monocytes per Patient"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.position = c(0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, max(patient_freq$proportion, na.rm = TRUE) * 1.2))

png(file.path(outdir, "fig6_C_mono_freq.png"), width = 12, height = 6, units = "in", res = 300)
print(p_freq)
dev.off()
cat("Saved: fig6_C_mono_freq.png\n")


# ============================================================================
# Panel D: Monocyte cNMF GEPs
# ============================================================================
cat("\n=== Panel D: cNMF GEPs ===\n")

gep_data <- read.csv("/data/user3/sobj/stroke_hc_v8_2/unbiased/cnmf/downstream/gep_named.csv",
                      stringsAsFactors = FALSE)

# Filter for monocyte-related cell types
mono_types <- c("CD14plus_Monocyte", "Inflammatory_Monocyte", "CD16plus_Monocyte")
gep_mono <- gep_data %>% filter(cell_type %in% mono_types)

cat("Mono GEPs:", nrow(gep_mono), "\n")
cat("Unique cell types:", unique(gep_mono$cell_type), "\n")

# Focus on L1 condition associations (Stroke vs HC)
gep_plot <- gep_mono %>%
  filter(!is.na(logFC_L1), !is.na(padj_L1)) %>%
  mutate(
    neg_log10_padj = -log10(pmax(padj_L1, 1e-300)),
    direction = ifelse(logFC_L1 > 0, "Up in Stroke", "Down in Stroke"),
    sig = padj_L1 < 0.05,
    cell_type_clean = gsub("plus_", "+ ", cell_type),
    cell_type_clean = gsub("_", " ", cell_type_clean),
    gep_label = paste0(cell_type_clean, "\n", gep_name_unique)
  )

# Sort by logFC
gep_plot <- gep_plot %>%
  arrange(logFC_L1) %>%
  mutate(gep_label = factor(gep_label, levels = unique(gep_label)))

# Only show significant GEPs or all if too few
if(sum(gep_plot$sig) >= 5) {
  gep_plot_sig <- gep_plot %>% filter(sig)
} else {
  gep_plot_sig <- gep_plot
}

p_cnmf <- ggplot(gep_plot_sig, aes(x = logFC_L1, y = gep_label)) +
  geom_point(aes(size = neg_log10_padj, fill = direction),
             shape = 21, color = "black", stroke = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(
    values = c("Up in Stroke" = color_stroke, "Down in Stroke" = color_hc),
    name = "Direction"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = expression(-log[10](p[adj]))
  ) +
  labs(
    title = "Monocyte Gene Expression Programs (cNMF)",
    subtitle = "Layer 1: Stroke vs HC",
    x = expression(log[2]~FC),
    y = ""
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

png(file.path(outdir, "fig6_D_mono_cnmf.png"), width = 10, height = 10, units = "in", res = 300)
print(p_cnmf)
dev.off()
cat("Saved: fig6_D_mono_cnmf.png\n")

cat("\n=== Figure 6 complete ===\n")
cat("Output files:\n")
cat("  ", file.path(outdir, "fig6_A_mono_umap.png"), "\n")
cat("  ", file.path(outdir, "fig6_B_mono_dotplot.png"), "\n")
cat("  ", file.path(outdir, "fig6_C_mono_freq.png"), "\n")
cat("  ", file.path(outdir, "fig6_D_mono_cnmf.png"), "\n")
