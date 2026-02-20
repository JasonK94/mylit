#!/usr/bin/env Rscript
# ==============================================================================
# FGS Results Visualization
# Generates plots for FGS n=50 results (whole_IS_50)
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

cat("=== FGS Results Visualization ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
FGS_DIR <- file.path(BASE, "fgs")
PLOT_DIR <- file.path(FGS_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

PREFIX <- "whole_IS_50"

# ---- Load data ----
cat("[1/6] Loading data...\n")
fgs <- qread(file.path(FGS_DIR, paste0(PREFIX, "_fgs.qs")))
tml <- qread(file.path(FGS_DIR, paste0(PREFIX, "_tml.qs")))
cmgi <- qread(file.path(FGS_DIR, paste0(PREFIX, "_cmgi.qs")))
cmgi_genes <- read.csv(file.path(FGS_DIR, paste0(PREFIX, "_cmgi_genes.csv")))

cat("  FGS methods:", length(fgs), "\n")
cat("  CMGI genes:", nrow(cmgi_genes), "\n")

# ---- Plot 1: Method AUC comparison ----
cat("[2/6] Method AUC barplot...\n")
method_perf <- data.frame(
  method = character(),
  auc = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(fgs)) {
  perf <- fgs[[nm]]$performance
  auc_val <- NA
  if (!is.null(perf$auc)) auc_val <- perf$auc
  else if (!is.null(perf$AUC)) auc_val <- perf$AUC
  else if (!is.null(perf$roc_auc)) auc_val <- perf$roc_auc
  method_perf <- rbind(method_perf, data.frame(
    method = nm,
    auc = auc_val,
    stringsAsFactors = FALSE
  ))
}

# Filter out NAs and failed methods
method_perf <- method_perf %>% filter(!is.na(auc))

if (nrow(method_perf) > 0) {
  method_perf <- method_perf %>% arrange(desc(auc))
  method_perf$method <- factor(method_perf$method, levels = method_perf$method)

  p1 <- ggplot(method_perf, aes(x = method, y = auc, fill = auc)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = round(auc, 3)), vjust = -0.3, size = 3.5) +
    scale_fill_gradient(low = "#4575b4", high = "#d73027", guide = "none") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    coord_cartesian(ylim = c(0, max(method_perf$auc, na.rm = TRUE) * 1.1)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "FGS L1 Method Performance (AUC)",
      subtitle = paste0("n_features = 50 | Target: g3 (IS only) | ", nrow(method_perf), " methods"),
      x = NULL, y = "AUC"
    )

  ggsave(file.path(PLOT_DIR, "01_method_auc_comparison.png"),
         p1, width = 8, height = 5, dpi = 300)
  cat("  Saved: 01_method_auc_comparison.png\n")
} else {
  cat("  WARNING: No AUC values found in FGS results\n")
  # Try to extract performance differently
  cat("  Available performance fields per method:\n")
  for (nm in names(fgs)) {
    cat("   ", nm, ":", paste(names(fgs[[nm]]$performance), collapse = ", "), "\n")
  }
}

# ---- Plot 2: TML model comparison (boxplot of CV metrics) ----
cat("[3/6] TML model comparison...\n")

if (!is.null(tml$model_comparison)) {
  df <- tml$model_comparison$values
  metric <- if (!is.null(tml$best_metric_name)) tml$best_metric_name else "ROC"

  cols <- colnames(df)
  metric_cols <- cols[grep(paste0("~", metric, "$"), cols)]

  if (length(metric_cols) > 0) {
    plot_data <- data.frame()
    for (col in metric_cols) {
      model_name <- sub(paste0("~", metric), "", col)
      vals <- df[[col]]
      plot_data <- rbind(plot_data, data.frame(
        Model = model_name,
        Value = vals,
        Fold = df$Resample
      ))
    }

    # Compute means for annotation
    means <- plot_data %>%
      group_by(Model) %>%
      summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_val))

    plot_data$Model <- factor(plot_data$Model, levels = means$Model)
    means$Model <- factor(means$Model, levels = means$Model)

    p2 <- ggplot(plot_data, aes(x = Model, y = Value, fill = Model)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.15, height = 0, size = 2.5, alpha = 0.7) +
      geom_text(data = means, aes(x = Model, y = mean_val,
                label = sprintf("%.3f", mean_val)),
                vjust = -1.5, size = 3.5, fontface = "bold") +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(
        title = paste0("TML L2 Model Comparison (", metric, ")"),
        subtitle = paste0("5-fold CV | Best: ", tml$best_model_name),
        x = NULL, y = metric
      )

    ggsave(file.path(PLOT_DIR, "02_tml_model_comparison.png"),
           p2, width = 7, height = 5, dpi = 300)
    cat("  Saved: 02_tml_model_comparison.png\n")
  }

  # Also plot all metrics together
  all_metrics <- c("ROC", "Sens", "Spec")
  all_plot_data <- data.frame()
  for (m in all_metrics) {
    m_cols <- cols[grep(paste0("~", m, "$"), cols)]
    for (col in m_cols) {
      model_name <- sub(paste0("~", m), "", col)
      vals <- df[[col]]
      all_plot_data <- rbind(all_plot_data, data.frame(
        Model = model_name, Value = vals, Metric = m,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(all_plot_data) > 0) {
    p2b <- ggplot(all_plot_data, aes(x = Model, y = Value, fill = Model)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
      facet_wrap(~Metric, scales = "free_y") +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(
        title = "TML L2: All Metrics by Model",
        subtitle = "5-fold CV",
        x = NULL, y = "Value"
      )

    ggsave(file.path(PLOT_DIR, "03_tml_all_metrics.png"),
           p2b, width = 10, height = 5, dpi = 300)
    cat("  Saved: 03_tml_all_metrics.png\n")
  }
}

# ---- Plot 3: Signature importance ----
cat("[4/6] Signature importance...\n")

sig_imp <- cmgi$signature_importance
if (!is.null(sig_imp)) {
  sig_df <- data.frame(
    signature = names(sig_imp),
    importance = as.numeric(sig_imp),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(importance))

  sig_df$signature <- factor(sig_df$signature, levels = sig_df$signature)

  p3 <- ggplot(sig_df, aes(x = signature, y = importance, fill = importance)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = round(importance, 3)),
              hjust = ifelse(sig_df$importance >= 0, -0.1, 1.1), size = 3.5) +
    scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027",
                         midpoint = 0, guide = "none") +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = "L1 Signature Importance in Meta-Learner",
      subtitle = paste0("Model: ", tml$best_model_name, " | Normalization: max_abs"),
      x = NULL, y = "Importance"
    )

  ggsave(file.path(PLOT_DIR, "04_signature_importance.png"),
         p3, width = 8, height = 5, dpi = 300)
  cat("  Saved: 04_signature_importance.png\n")
}

# ---- Plot 4: CMGI top gene importance (bar) ----
cat("[5/6] CMGI gene importance...\n")

top_n <- 30
top_genes <- cmgi_genes %>%
  arrange(desc(abs_contribution)) %>%
  head(top_n)

top_genes$direction <- ifelse(top_genes$contribution > 0, "UP in g3=2 (Bad)", "DOWN in g3=2 (Bad)")
top_genes$gene <- factor(top_genes$gene, levels = rev(top_genes$gene))

p4 <- ggplot(top_genes, aes(x = gene, y = contribution, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("UP in g3=2 (Bad)" = "#d73027", "DOWN in g3=2 (Bad)" = "#4575b4"),
                    name = "Direction") +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  ) +
  labs(
    title = paste0("Top ", top_n, " CMGI Genes (Meta Gene Importance)"),
    subtitle = paste0("Total: ", nrow(cmgi_genes), " genes | n_features = 50 | Target: g3"),
    x = NULL, y = "Meta-importance (signed)"
  )

ggsave(file.path(PLOT_DIR, "05_cmgi_top_genes.png"),
       p4, width = 8, height = 8, dpi = 300)
cat("  Saved: 05_cmgi_top_genes.png\n")

# Plot 4b: Top 50 genes heatmap-style (gene × method contribution)
if (!is.null(cmgi$gene_importance)) {
  gi <- cmgi$gene_importance
  top50_genes <- cmgi_genes %>% arrange(desc(abs_contribution)) %>% head(50) %>% pull(gene)

  gi_top <- gi %>% filter(gene %in% top50_genes)

  if (nrow(gi_top) > 0) {
    # Create contribution matrix
    gi_wide <- gi_top %>%
      select(gene, signature, contribution) %>%
      pivot_wider(names_from = signature, values_from = contribution, values_fill = 0)

    gi_long <- gi_top %>%
      mutate(gene = factor(gene, levels = rev(top50_genes)))

    p4b <- ggplot(gi_long, aes(x = signature, y = gene, fill = contribution)) +
      geom_tile() +
      scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                           midpoint = 0, name = "Contribution") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7)
      ) +
      labs(
        title = "Gene × Method Contribution Heatmap (Top 50)",
        subtitle = "Signed contribution of each gene through each L1 method",
        x = "L1 Method", y = NULL
      )

    ggsave(file.path(PLOT_DIR, "06_cmgi_gene_method_heatmap.png"),
           p4b, width = 10, height = 12, dpi = 300)
    cat("  Saved: 06_cmgi_gene_method_heatmap.png\n")
  }
}

# ---- Plot 5: Meta-score distribution by g3 ----
cat("[6/6] Meta-score distribution...\n")

sobj_file <- file.path(FGS_DIR, paste0(PREFIX, "_sobj.qs"))
if (file.exists(sobj_file)) {
  cat("  Loading Seurat object for meta-score...\n")
  sobj <- qread(sobj_file)
  meta <- sobj@meta.data

  # Find meta-score column
  score_cols <- grep("fgs_meta|meta_score|AMSC|amsc", colnames(meta), value = TRUE)
  cat("  Score columns found:", paste(score_cols, collapse = ", "), "\n")

  if (length(score_cols) > 0) {
    score_col <- score_cols[1]
    plot_df <- meta %>%
      filter(!is.na(g3)) %>%
      select(all_of(c(score_col, "g3"))) %>%
      mutate(g3_label = ifelse(g3 == 1, "g3=1 (Good)", "g3=2 (Bad)"))

    p5 <- ggplot(plot_df, aes(x = g3_label, y = .data[[score_col]], fill = g3_label)) +
      geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
      scale_fill_manual(values = c("g3=1 (Good)" = "#4575b4", "g3=2 (Bad)" = "#d73027"),
                        guide = "none") +
      theme_minimal(base_size = 12) +
      labs(
        title = "FGS Meta-Score Distribution by g3",
        subtitle = paste0("Score: ", score_col, " | n = ", nrow(plot_df), " cells"),
        x = NULL, y = "Meta-score"
      )

    ggsave(file.path(PLOT_DIR, "07_meta_score_distribution.png"),
           p5, width = 6, height = 5, dpi = 300)
    cat("  Saved: 07_meta_score_distribution.png\n")

    # Plot 5b: Meta-score by anno1 cell type
    if ("anno1" %in% colnames(meta)) {
      plot_df2 <- meta %>%
        filter(!is.na(g3)) %>%
        select(all_of(c(score_col, "g3", "anno1"))) %>%
        mutate(g3_label = ifelse(g3 == 1, "Good", "Bad"))

      # Order cell types by median score difference
      ct_order <- plot_df2 %>%
        group_by(anno1) %>%
        summarise(med_score = median(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
        arrange(med_score) %>%
        pull(anno1)

      plot_df2$anno1 <- factor(plot_df2$anno1, levels = ct_order)

      p5b <- ggplot(plot_df2, aes(x = anno1, y = .data[[score_col]], fill = g3_label)) +
        geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
        scale_fill_manual(values = c("Good" = "#4575b4", "Bad" = "#d73027"), name = "g3") +
        coord_flip() +
        theme_minimal(base_size = 11) +
        labs(
          title = "FGS Meta-Score by Cell Type and g3",
          subtitle = paste0("n = ", nrow(plot_df2), " cells"),
          x = NULL, y = "Meta-score"
        )

      ggsave(file.path(PLOT_DIR, "08_meta_score_by_celltype.png"),
             p5b, width = 9, height = 7, dpi = 300)
      cat("  Saved: 08_meta_score_by_celltype.png\n")
    }
  }
  rm(sobj); gc(verbose = FALSE)
} else {
  cat("  Seurat object not found, skipping meta-score plots\n")
}

# ---- Done ----
cat("\n=== FGS Visualization Complete ===\n")
cat("Output:", PLOT_DIR, "\n")
n_plots <- length(list.files(PLOT_DIR, pattern = "\\.png$"))
cat("Total plots:", n_plots, "\n")
cat("Finished:", format(Sys.time()), "\n")
