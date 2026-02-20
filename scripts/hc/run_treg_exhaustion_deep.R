#!/usr/bin/env Rscript
# ==============================================================================
# Treg Exhaustion Markers × Age: Deep Analysis
# ==============================================================================
# Q1: White background PNGs
# Q2: Age group (Young/Middle/Old) stratified
# Q3: Sex-stratified
# Q4: Treg subclustering → marker expression
# Focuses on analyzable markers (>20% non-zero at pseudobulk)
# ==============================================================================

cat("
================================================================================
  Treg Exhaustion Markers × Age — Deep Analysis
  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
================================================================================
\n")

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(patchwork)
  library(data.table)
})

set.seed(42)

OUTPUT_DIR <- "/data/user3/sobj/hc_only_v1/treg_exhaustion_age"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# White background theme
theme_white <- theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA))

# All markers
all_markers <- c("PTGFRN","NR4A3","EGR2","EGR3","NR4A1","ZBED2",
                  "IL1RN","ITPKA","ITGB8","GLT1D1","PTGS1","GEM",
                  "PAPSS2","PLEKHA7","SLC16A14")

# ==============================================================================
# 1. Load Seurat object (need full single-cell for subclustering)
# ==============================================================================
message("[1/6] Loading Seurat object...")
t0 <- Sys.time()
obj <- qread("/data/user3/sobj/hc_only_v1/2_hc_annotated.qs")
message(sprintf("  %d genes × %d cells (%.0fs)", nrow(obj), ncol(obj),
                as.numeric(Sys.time() - t0, units = "secs")))

# Treg subset
treg <- subset(obj, anno1 == "CD4 Treg")
message(sprintf("  Treg: %d cells, %d patients",
                ncol(treg), length(unique(treg$patient_name))))

# Check which markers exist
found_markers <- all_markers[all_markers %in% rownames(treg)]
missing_markers <- setdiff(all_markers, found_markers)
message(sprintf("  Markers found: %d / %d", length(found_markers), length(all_markers)))
if (length(missing_markers) > 0) {
  message(sprintf("  Missing: %s", paste(missing_markers, collapse = ", ")))
}

# Categorize by sparsity
pb_data <- readRDS("/data/user3/sobj/hc_only_v1/fgs_continuous_treg/pseudobulk_data.rds")
pb_counts <- pb_data$pb_counts
pct_zero <- sapply(found_markers, function(g) 100 * sum(pb_counts[g, ] == 0) / ncol(pb_counts))
analyzable <- found_markers[pct_zero < 60]  # <60% zero at pseudobulk
sparse_markers <- found_markers[pct_zero >= 60]
message(sprintf("  Analyzable (pseudobulk <60%% zero): %s", paste(analyzable, collapse = ", ")))
message(sprintf("  Sparse (>60%% zero): %s", paste(sparse_markers, collapse = ", ")))

# ==============================================================================
# 2. Q1 + Q2 + Q3: Pseudobulk scatter with white BG, age groups, sex strata
# ==============================================================================
message("\n[2/6] Pseudobulk analysis: age group & sex stratification...")

pb_logcpm <- pb_data$pb_logcpm
patient_meta <- pb_data$patient_meta

# --- Q1: Redo scatter with white background ---
cat("  [Q1] Regenerating scatter plots with white background...\n")

# Compute age_group if not present
if (!"age_group" %in% colnames(patient_meta)) {
  patient_meta$age_group <- ifelse(patient_meta$age < 35, "Young",
                             ifelse(patient_meta$age <= 50, "Middle", "Old"))
}

plot_data <- do.call(rbind, lapply(found_markers, function(g) {
  data.frame(
    gene = g, patient = rownames(patient_meta),
    age = patient_meta$age, sex = patient_meta$sex,
    age_group = patient_meta$age_group,
    expr = pb_logcpm[g, rownames(patient_meta)],
    raw_count = pb_counts[g, rownames(patient_meta)],
    stringsAsFactors = FALSE
  )
}))

# Spearman stats per gene
sp_stats <- plot_data %>%
  group_by(gene) %>%
  summarise(
    rho = cor(expr, age, method = "spearman"),
    p = cor.test(expr, age, method = "spearman", exact = FALSE)$p.value,
    pct_zero = 100 * mean(raw_count == 0),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%s\nrho=%+.2f p=%.1e\n%.0f%% zero", gene, rho, p, pct_zero))

plot_data <- merge(plot_data, sp_stats[, c("gene", "label", "rho")], by = "gene")
gene_order <- sp_stats$gene[order(sp_stats$rho)]
plot_data$gene_f <- factor(plot_data$label,
                            levels = sp_stats$label[match(gene_order, sp_stats$gene)])

ncol_p <- 4
nrow_p <- ceiling(length(found_markers) / ncol_p)

p1 <- ggplot(plot_data, aes(x = age, y = expr)) +
  geom_point(aes(color = sex), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, alpha = 0.3) +
  facet_wrap(~ gene_f, scales = "free_y", ncol = ncol_p) +
  scale_color_manual(values = c("F" = "#D6604D", "M" = "#4393C3")) +
  theme_white +
  labs(title = "Treg Exhaustion Markers × Age (pseudobulk, n=96 HC)",
       x = "Age (years)", y = "logCPM", color = "Sex") +
  theme(strip.text = element_text(face = "bold", size = 8), legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "scatter_all_markers_vs_age_white.png"),
       p1, width = 4 * ncol_p, height = 3.5 * nrow_p, dpi = 200,
       bg = "white", limitsize = FALSE)
message("  Saved: scatter_all_markers_vs_age_white.png")

# --- Q2: Age group boxplots ---
cat("\n  [Q2] Age group stratified analysis...\n")

# Make sure age_group is ordered
plot_data$age_group <- factor(plot_data$age_group, levels = c("Young", "Middle", "Old"))

# Statistical test: Kruskal-Wallis + pairwise Wilcoxon
age_grp_stats <- plot_data %>%
  group_by(gene) %>%
  summarise(
    kw_p = kruskal.test(expr ~ age_group)$p.value,
    yw_mean = mean(expr[age_group == "Young"]),
    mw_mean = mean(expr[age_group == "Middle"]),
    ow_mean = mean(expr[age_group == "Old"]),
    yo_wilcox_p = tryCatch(wilcox.test(expr[age_group == "Young"],
                                        expr[age_group == "Old"])$p.value,
                           error = function(e) NA),
    .groups = "drop"
  )

cat("\n  Age group comparison (Kruskal-Wallis + Young-vs-Old Wilcoxon):\n")
cat(sprintf("  %-15s %10s %10s %10s %10s %12s\n",
            "Gene", "Young_mean", "Mid_mean", "Old_mean", "KW_p", "Y-vs-O_p"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_len(nrow(age_grp_stats))) {
  r <- age_grp_stats[i, ]
  sig <- ifelse(r$kw_p < 0.05, "*", "")
  cat(sprintf("  %-15s %10.3f %10.3f %10.3f %10.3e %12.3e %s\n",
              r$gene, r$yw_mean, r$mw_mean, r$ow_mean, r$kw_p, r$yo_wilcox_p, sig))
}

p2 <- ggplot(plot_data, aes(x = age_group, y = expr, fill = age_group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_jitter(aes(color = sex), width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Young" = "#4DAF4A", "Middle" = "#FF7F00", "Old" = "#E41A1C")) +
  scale_color_manual(values = c("F" = "#D6604D", "M" = "#4393C3")) +
  theme_white +
  labs(title = "Treg Exhaustion Markers by Age Group (pseudobulk, HC n=96)",
       subtitle = "Young(<35) n=45, Middle(35-50) n=28, Old(>50) n=23",
       x = "", y = "logCPM") +
  theme(strip.text = element_text(face = "bold", size = 9), legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "boxplot_markers_by_age_group.png"),
       p2, width = 20, height = 3.5 * ceiling(length(found_markers) / 5), dpi = 200,
       bg = "white", limitsize = FALSE)
message("  Saved: boxplot_markers_by_age_group.png")

# --- Q3: Sex-stratified scatter ---
cat("\n  [Q3] Sex-stratified analysis...\n")

sex_stats <- plot_data %>%
  group_by(gene, sex) %>%
  summarise(
    rho = cor(expr, age, method = "spearman"),
    p = tryCatch(cor.test(expr, age, method = "spearman", exact = FALSE)$p.value,
                 error = function(e) NA),
    n = n(),
    .groups = "drop"
  )

cat("\n  Sex-stratified Spearman (age vs expression):\n")
cat(sprintf("  %-15s  %6s %8s %10s  | %6s %8s %10s\n",
            "Gene", "F_rho", "F_p", "F_n", "M_rho", "M_p", "M_n"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (g in found_markers) {
  f_row <- sex_stats[sex_stats$gene == g & sex_stats$sex == "F", ]
  m_row <- sex_stats[sex_stats$gene == g & sex_stats$sex == "M", ]
  cat(sprintf("  %-15s  %+6.3f %8.2e %10d  | %+6.3f %8.2e %10d\n",
              g,
              ifelse(nrow(f_row) > 0, f_row$rho, NA),
              ifelse(nrow(f_row) > 0, f_row$p, NA),
              ifelse(nrow(f_row) > 0, f_row$n, 0),
              ifelse(nrow(m_row) > 0, m_row$rho, NA),
              ifelse(nrow(m_row) > 0, m_row$p, NA),
              ifelse(nrow(m_row) > 0, m_row$n, 0)))
}

p3 <- ggplot(plot_data, aes(x = age, y = expr, color = sex)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_color_manual(values = c("F" = "#D6604D", "M" = "#4393C3")) +
  theme_white +
  labs(title = "Treg Exhaustion Markers × Age, Sex-Stratified (pseudobulk, HC n=96)",
       subtitle = "Separate regression lines per sex | F=57, M=39",
       x = "Age (years)", y = "logCPM") +
  theme(strip.text = element_text(face = "bold", size = 9), legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "scatter_markers_sex_stratified.png"),
       p3, width = 20, height = 3.5 * ceiling(length(found_markers) / 5), dpi = 200,
       bg = "white", limitsize = FALSE)
message("  Saved: scatter_markers_sex_stratified.png")

# ==============================================================================
# 3. Q4: Treg subclustering
# ==============================================================================
message("\n[3/6] Treg subclustering...")
rm(pb_data, pb_counts, pb_logcpm); gc(verbose = FALSE)

# Use scVI embeddings from original object for subclustering
# DO NOT recompute PCA — use original scVI reduction
if ("integrated.scvi" %in% Reductions(treg)) {
  message("  Using integrated.scvi reduction from original object")
  treg <- FindNeighbors(treg, reduction = "integrated.scvi", dims = 1:30)
  treg <- FindClusters(treg, resolution = c(0.3, 0.5, 0.8, 1.0))

  # Use resolution 0.5 as default
  Idents(treg) <- "RNA_snn_res.0.5"
  n_clusters <- length(levels(Idents(treg)))
  message(sprintf("  Subclusters (res=0.5): %d", n_clusters))

  # Print cluster sizes
  cluster_sizes <- sort(table(Idents(treg)), decreasing = TRUE)
  for (i in seq_along(cluster_sizes)) {
    cat(sprintf("    Cluster %s: %d cells\n", names(cluster_sizes)[i], cluster_sizes[i]))
  }

  # FindMarkers for subclusters (quick version)
  message("  Running FindAllMarkers for Treg subclusters...")
  treg_markers_all <- FindAllMarkers(treg, only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.25, max.cells.per.ident = 500)

  cat("\n  Top 5 markers per Treg subcluster:\n")
  top5 <- treg_markers_all %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5) %>%
    summarise(genes = paste(gene, collapse = ", "), .groups = "drop")
  for (i in seq_len(nrow(top5))) {
    cat(sprintf("    Cluster %s: %s\n", top5$cluster[i], top5$genes[i]))
  }

  # Save subcluster markers
  write.csv(treg_markers_all, file.path(OUTPUT_DIR, "treg_subcluster_markers.csv"),
            row.names = FALSE)

} else if ("umap.scvi" %in% Reductions(treg)) {
  message("  integrated.scvi not available, using umap.scvi for neighbors")
  treg <- FindNeighbors(treg, reduction = "umap.scvi", dims = 1:2)
  treg <- FindClusters(treg, resolution = c(0.3, 0.5, 0.8))
  Idents(treg) <- "RNA_snn_res.0.5"
  n_clusters <- length(levels(Idents(treg)))
  message(sprintf("  Subclusters (res=0.5): %d", n_clusters))
} else {
  message("  WARNING: No scVI reduction found. Using PCA (not ideal).")
  treg <- NormalizeData(treg) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  treg <- FindNeighbors(treg, dims = 1:20)
  treg <- FindClusters(treg, resolution = 0.5)
  n_clusters <- length(levels(Idents(treg)))
}

treg$treg_subcluster <- Idents(treg)

# --- Exhaustion markers across subclusters ---
message("\n  Exhaustion marker expression by subcluster...")

# Build per-subcluster pseudobulk per patient
treg_meta <- treg@meta.data
treg_meta$sub_patient <- paste0(treg_meta$treg_subcluster, "__", treg_meta$patient_name)

cat("\n  Mean expression (logCPM) per subcluster (analyzable markers):\n")
# Get single-cell expression
expr_mat <- GetAssayData(treg, layer = "data")

# Per-subcluster mean expression
sub_expr <- list()
for (cl in sort(unique(treg$treg_subcluster))) {
  cells <- colnames(treg)[treg$treg_subcluster == cl]
  sub_expr[[as.character(cl)]] <- rowMeans(expr_mat[found_markers, cells, drop = FALSE])
}
sub_expr_df <- as.data.frame(do.call(cbind, sub_expr))
colnames(sub_expr_df) <- paste0("Cl_", names(sub_expr))

cat(sprintf("  %-15s %s\n", "Gene", paste(colnames(sub_expr_df), collapse = "  ")))
cat(paste(rep("-", 15 + ncol(sub_expr_df) * 8), collapse = ""), "\n")
for (g in found_markers) {
  vals <- sprintf("%7.3f", as.numeric(sub_expr_df[g, ]))
  cat(sprintf("  %-15s %s\n", g, paste(vals, collapse = " ")))
}

# --- Per-subcluster age association (pseudobulk within subcluster) ---
message("\n  Per-subcluster age association...")

subcluster_age_results <- list()
for (cl in sort(unique(treg$treg_subcluster))) {
  cells_cl <- colnames(treg)[treg$treg_subcluster == cl]
  meta_cl <- treg@meta.data[cells_cl, ]
  counts_cl <- GetAssayData(treg, layer = "counts")[, cells_cl]

  # Patient-level pseudobulk within this subcluster
  patients_cl <- unique(meta_cl$patient_name)
  pat_info <- meta_cl %>%
    group_by(patient_name) %>%
    summarise(age = first(age), sex = first(sex), n_cells = n(), .groups = "drop") %>%
    filter(n_cells >= 3) %>%
    as.data.frame()

  if (nrow(pat_info) < 20) next  # Need enough patients

  keep_pats <- pat_info$patient_name
  cell_pat <- meta_cl$patient_name
  names(cell_pat) <- rownames(meta_cl)
  cell_pat <- cell_pat[cell_pat %in% keep_pats]

  pats_ordered <- sort(keep_pats)
  pi_mat <- Matrix::sparseMatrix(
    i = match(cell_pat, pats_ordered),
    j = seq_along(cell_pat), x = 1,
    dims = c(length(pats_ordered), length(cell_pat)),
    dimnames = list(pats_ordered, names(cell_pat))
  )
  pb_cl <- as.matrix(counts_cl[found_markers, names(cell_pat)] %*% Matrix::t(pi_mat))

  pat_info_ordered <- pat_info[match(pats_ordered, pat_info$patient_name), ]
  rownames(pat_info_ordered) <- pat_info_ordered$patient_name

  # logCPM (simple normalization)
  lib_sizes <- colSums(counts_cl[, names(cell_pat)] %*% Matrix::t(pi_mat))
  pb_logcpm_cl <- log2(t(t(pb_cl) / lib_sizes * 1e6) + 1)

  for (g in found_markers) {
    if (sum(pb_cl[g, ] > 0) < 5) next  # Too sparse
    rho <- cor(pb_logcpm_cl[g, ], pat_info_ordered$age, method = "spearman")
    p_val <- tryCatch(
      cor.test(pb_logcpm_cl[g, ], pat_info_ordered$age, method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA
    )
    subcluster_age_results[[length(subcluster_age_results) + 1]] <- data.frame(
      subcluster = cl, gene = g, rho = rho, pvalue = p_val,
      n_patients = nrow(pat_info_ordered),
      pct_nonzero = 100 * sum(pb_cl[g, ] > 0) / ncol(pb_cl)
    )
  }
}

sub_age_df <- do.call(rbind, subcluster_age_results)
if (nrow(sub_age_df) > 0) {
  sub_age_df$padj <- p.adjust(sub_age_df$pvalue, method = "BH")
  sub_age_df <- sub_age_df[order(sub_age_df$pvalue), ]

  cat("\n  === Per-subcluster age associations (Spearman) ===\n")
  cat(sprintf("  %-10s %-15s %8s %10s %10s %8s %10s\n",
              "Cluster", "Gene", "rho", "p-value", "padj", "n_pat", "pct_nz"))
  cat(paste(rep("-", 80), collapse = ""), "\n")
  sig_sub <- sub_age_df[sub_age_df$pvalue < 0.05, ]
  if (nrow(sig_sub) > 0) {
    for (i in seq_len(min(30, nrow(sig_sub)))) {
      r <- sig_sub[i, ]
      cat(sprintf("  %-10s %-15s %+8.3f %10.2e %10.2e %8d %9.1f%%\n",
                  r$subcluster, r$gene, r$rho, r$pvalue, r$padj,
                  r$n_patients, r$pct_nonzero))
    }
  } else {
    cat("  No nominally significant associations found.\n")
  }

  write.csv(sub_age_df, file.path(OUTPUT_DIR, "subcluster_age_associations.csv"),
            row.names = FALSE)
}

# --- Subcluster visualization ---
message("\n  Generating subcluster plots...")

# DotPlot of markers across subclusters
# Exclude GEM gene since it conflicts with GEM metadata column
dotplot_markers <- setdiff(found_markers, "GEM")
p_dot_sub <- DotPlot(treg, features = dotplot_markers, group.by = "treg_subcluster") +
  coord_flip() +
  theme_white +
  labs(title = "Exhaustion Markers by Treg Subcluster (excl. GEM gene)")
ggsave(file.path(OUTPUT_DIR, "dotplot_markers_by_subcluster.png"),
       p_dot_sub, width = max(8, n_clusters * 0.8), height = max(6, length(dotplot_markers) * 0.35),
       dpi = 200, bg = "white")
message("  Saved: dotplot_markers_by_subcluster.png")

# Subcluster composition by age_group
comp_df <- treg@meta.data %>%
  group_by(treg_subcluster, age_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_group) %>%
  mutate(pct = 100 * n / sum(n))

comp_df$age_group <- factor(comp_df$age_group, levels = c("Young", "Middle", "Old"))

p_comp <- ggplot(comp_df, aes(x = age_group, y = pct, fill = treg_subcluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_white +
  labs(title = "Treg Subcluster Composition by Age Group",
       x = "", y = "% of Treg cells", fill = "Subcluster") +
  theme(legend.position = "right")
ggsave(file.path(OUTPUT_DIR, "barplot_subcluster_by_age_group.png"),
       p_comp, width = 8, height = 6, dpi = 200, bg = "white")
message("  Saved: barplot_subcluster_by_age_group.png")

# Subcluster composition by sex
comp_sex <- treg@meta.data %>%
  group_by(treg_subcluster, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sex) %>%
  mutate(pct = 100 * n / sum(n))

p_comp_sex <- ggplot(comp_sex, aes(x = sex, y = pct, fill = treg_subcluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_white +
  labs(title = "Treg Subcluster Composition by Sex",
       x = "", y = "% of Treg cells", fill = "Subcluster") +
  theme(legend.position = "right")
ggsave(file.path(OUTPUT_DIR, "barplot_subcluster_by_sex.png"),
       p_comp_sex, width = 6, height = 6, dpi = 200, bg = "white")

# --- Also redo heatmap and dotplot from Q1 with white bg ---
message("\n  Regenerating original plots with white background...")

# Reload pseudobulk for heatmap
pb_data <- readRDS("/data/user3/sobj/hc_only_v1/fgs_continuous_treg/pseudobulk_data.rds")
pb_logcpm <- pb_data$pb_logcpm
patient_meta <- pb_data$patient_meta
pb_counts <- pb_data$pb_counts
# Recompute age_group
if (!"age_group" %in% colnames(patient_meta)) {
  patient_meta$age_group <- ifelse(patient_meta$age < 35, "Young",
                             ifelse(patient_meta$age <= 50, "Middle", "Old"))
}

mat_markers <- pb_logcpm[found_markers, ]
mat_z <- t(scale(t(mat_markers)))
age_order <- order(patient_meta$age)
mat_z <- mat_z[, age_order]
mat_z[mat_z > 3] <- 3; mat_z[mat_z < -3] <- -3

# Order rows by Spearman rho
sp <- sp_stats[match(found_markers, sp_stats$gene), ]
row_order <- sp$gene[order(sp$rho)]
mat_z <- mat_z[row_order, ]

anno_col <- data.frame(
  Age = patient_meta$age[age_order],
  Sex = patient_meta$sex[age_order],
  Age_group = patient_meta$age_group[age_order],
  row.names = colnames(mat_z)
)
ann_colors <- list(
  Sex = c("F" = "#D6604D", "M" = "#4393C3"),
  Age = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  Age_group = c("Young" = "#4DAF4A", "Middle" = "#FF7F00", "Old" = "#E41A1C")
)

png(file.path(OUTPUT_DIR, "heatmap_exhaustion_markers_age_white.png"),
    width = 1600, height = max(500, length(found_markers) * 40), res = 150)
par(bg = "white")
pheatmap(mat_z, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE,
         annotation_col = anno_col, annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         main = "Treg Exhaustion Markers (z-scored logCPM, patients ordered by age)",
         fontsize_row = 10, gaps_col = NULL)
dev.off()
message("  Saved: heatmap_exhaustion_markers_age_white.png")

# Dotplot with white bg
dot_df <- data.frame(
  gene = sp$gene,
  rho = sp$rho,
  neg_log10_p = -log10(sp$p),
  pct_zero = sp$pct_zero,
  sig = ifelse(sp$p < 0.05, "p<0.05", "NS"),
  stringsAsFactors = FALSE
)
dot_df$gene <- factor(dot_df$gene, levels = dot_df$gene[order(dot_df$rho)])

p_dot <- ggplot(dot_df, aes(x = rho, y = gene)) +
  geom_point(aes(size = neg_log10_p, color = pct_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_color_gradient(low = "steelblue", high = "red3", name = "% zero\n(pseudobulk)") +
  scale_size_continuous(range = c(2, 8), name = "-log10(p)") +
  theme_white +
  labs(title = "Treg Exhaustion Markers: Age Association (Spearman)",
       subtitle = "HC-only, n=96 patients | size=-log10(p) | color=% zero at pseudobulk",
       x = "Spearman rho (+ = up with age)", y = "")
ggsave(file.path(OUTPUT_DIR, "dotplot_markers_rho_white.png"),
       p_dot, width = 10, height = max(5, length(found_markers) * 0.35), dpi = 200, bg = "white")
message("  Saved: dotplot_markers_rho_white.png")

# ==============================================================================
# 4. Summary
# ==============================================================================
message("\n[6/6] Summary")

cat("\n========================================================================\n")
cat("  EXHAUSTION MARKERS × AGE DEEP ANALYSIS SUMMARY\n")
cat("========================================================================\n\n")

cat("  Q1: White-background plots generated ✓\n\n")

cat("  Q2: Age group stratification:\n")
sig_kw <- age_grp_stats[age_grp_stats$kw_p < 0.05, ]
if (nrow(sig_kw) > 0) {
  for (i in seq_len(nrow(sig_kw))) {
    cat(sprintf("    %s: KW p=%.3e (Y=%.3f, M=%.3f, O=%.3f)\n",
                sig_kw$gene[i], sig_kw$kw_p[i],
                sig_kw$yw_mean[i], sig_kw$mw_mean[i], sig_kw$ow_mean[i]))
  }
} else {
  cat("    No markers significantly different across age groups (KW p<0.05)\n")
}

cat("\n  Q3: Sex-stratified:\n")
sex_sig <- sex_stats[sex_stats$p < 0.05 & !is.na(sex_stats$p), ]
if (nrow(sex_sig) > 0) {
  for (i in seq_len(nrow(sex_sig))) {
    cat(sprintf("    %s in %s: rho=%+.3f, p=%.3e\n",
                sex_sig$gene[i], sex_sig$sex[i], sex_sig$rho[i], sex_sig$p[i]))
  }
} else {
  cat("    No markers with significant age-association within either sex\n")
}

cat("\n  Q4: Treg subclustering:\n")
cat(sprintf("    %d subclusters identified (scVI, res=0.5)\n", n_clusters))
if (exists("sub_age_df") && nrow(sub_age_df) > 0) {
  sub_sig <- sub_age_df[sub_age_df$pvalue < 0.05, ]
  cat(sprintf("    %d subcluster×gene pairs with p<0.05 (of %d tested)\n",
              nrow(sub_sig), nrow(sub_age_df)))
  if (nrow(sub_sig) > 0) {
    for (i in seq_len(min(10, nrow(sub_sig)))) {
      r <- sub_sig[i, ]
      cat(sprintf("      Cl%s × %s: rho=%+.3f, p=%.3e\n",
                  r$subcluster, r$gene, r$rho, r$pvalue))
    }
  }
} else {
  cat("    No subcluster associations computed\n")
}

cat(sprintf("\n  Q5: Zero-expression issue:\n"))
cat(sprintf("    10/15 markers have >70%% zero patients at pseudobulk level\n"))
cat(sprintf("    These are NOT typically expressed in Treg (Treg-irrelevant markers)\n"))
cat(sprintf("    Analyzable: NR4A3, NR4A1, GEM, EGR2, PLEKHA7 (<60%% zero)\n"))

# Save Treg object with subclusters
qsave(treg, file.path(OUTPUT_DIR, "treg_subclustered.qs"))
message(sprintf("  Saved: treg_subclustered.qs (%d cells)", ncol(treg)))

cat(sprintf("\nDone! Time: %s\n", Sys.time()))
cat(sprintf("Output: %s\n", OUTPUT_DIR))
