# ============================================================================
# DEG Consensus - Full Pipeline & Visualisation (Phase 6)
# ============================================================================

#' Run full DEG consensus pipeline with optional plots
#'
#' @description
#' Convenience wrapper that:
#' \enumerate{
#'   \item runs \code{run_deg_consensus()} for multiple methods,
#'   \item standardises results (\code{standardize_deg_results()}),
#'   \item builds gene × method matrices (\code{build_deg_matrices()}),
#'   \item computes agreement and meta-analysis consensus scores, and
#'   \item (optionally) generates summary plots (distributions, volcano,
#'         heatmap, PCA/UMAP) and saves intermediate objects as \code{.qs}.
#' }
#'
#' By default, this function uses a very loose thresholding strategy:
#' per-method significance is based on raw p-values
#' (\code{significance_mode = "pvalue"}) and no filtering is applied at the
#' consensus level (\code{fdr_threshold = NULL}, \code{agreement_threshold = NULL},
#' \code{min_methods = NULL}), so that all genes can be ranked/evaluated
#' (e.g., for gene set analysis).
#'
#' @param sobj Seurat object
#' @param contrast Contrast string (e.g., "2 - 1")
#' @param methods Character vector of methods to run
#' @param cluster_id Column name for cluster/cell-type (default: "anno3.scvi")
#' @param sample_id Column name for sample ID (default: "hos_no")
#' @param group_id Column name for group/condition (default: "g3")
#' @param batch_id Optional batch column (default: "GEM")
#' @param significance_mode One of "pvalue" (default) or "fdr".
#'   Controls how per-method significance is defined in \code{build_deg_matrices()}.
#' @param pvalue_threshold Raw p-value threshold when
#'   \code{significance_mode = "pvalue"} (default: 0.05).
#' @param fdr_threshold FDR threshold applied at the meta-analysis level
#'   (meta_p_adj) in \code{generate_consensus_deg_list()}. If \code{NULL}
#'   (default), no FDR filtering is applied and all genes are kept.
#' @param agreement_threshold Minimum agreement score (0–1). If \code{NULL}
#'   (default), no agreement filtering is applied.
#' @param min_methods Minimum number of methods showing significance. If
#'   \code{NULL} (default), no filter is applied.
#' @param output_dir Directory to save intermediate .qs objects and plots
#'   (default: current working directory).
#' @param prefix Optional file prefix. If \code{NULL}, a prefix of the form
#'   \code{"deg_consensus_<ds|full>_YYYYMMDD_HHMMSS"} is generated based on
#'   the number of cells and current time.
#' @param save_intermediate Logical; if TRUE (default), saves intermediate
#'   objects (raw results, standardised results, matrices, consensus scores)
#'   as \code{.qs} files.
#' @param make_plots Logical; if TRUE (default), generates and saves plots
#'   (distributions, volcano, heatmap, PCA, UMAP) using ggplot2/uwot.
#'
#' @return A list with:
#'   \describe{
#'     \item{raw_results}{Output of \code{run_deg_consensus()}}
#'     \item{standardized_results}{List of standardised DEG tables}
#'     \item{deg_matrices}{Gene × method matrices from \code{build_deg_matrices()}}
#'     \item{agreement_scores}{Agreement scores (per gene)}
#'     \item{pca_result}{Method-level PCA result (from \code{perform_deg_pca()})}
#'     \item{consensus_scores}{Per-gene consensus statistics}
#'     \item{consensus_deg_list}{Filtered/ranked consensus DEG table}
#'     \item{output_prefix}{File prefix used for saved outputs}
#'   }
#'
#' @export
run_deg_consensus_with_plots <- function(
  sobj,
  contrast = "2 - 1",
  methods = c(
    "muscat-edgeR",
    "muscat-DESeq2",
    "muscat-limma-voom",
    "muscat-limma-trend",
    "limma-voom",
    "edgeR-LRT"
  ),
  cluster_id = "anno3.scvi",
  sample_id = "hos_no",
  group_id = "g3",
  batch_id = "GEM",
  significance_mode = c("pvalue", "fdr"),
  pvalue_threshold = 0.05,
  fdr_threshold = NULL,
  agreement_threshold = NULL,
  min_methods = NULL,
  output_dir = NULL,
  prefix = NULL,
  save_intermediate = TRUE,
  make_plots = TRUE
) {
  significance_mode <- match.arg(significance_mode)

  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Build a default prefix if not provided
  if (is.null(prefix)) {
    n_cells <- tryCatch(ncol(sobj), error = function(e) NA_integer_)
    dataset_tag <- if (!is.na(n_cells) && n_cells <= 3000) "ds" else "full"
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    prefix <- file.path(output_dir, paste0("deg_consensus_", dataset_tag, "_", timestamp))
  }

  if (!requireNamespace("qs", quietly = TRUE)) {
    warning("Package 'qs' is not installed; intermediate objects will not be saved.", call. = FALSE)
    save_intermediate <- FALSE
  }
  # Ensure plotting dependencies if needed
  if (make_plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is not installed; plots will be skipped.", call. = FALSE)
      make_plots <- FALSE
    } else {
      # Attach ggplot2 so that internal helpers (e.g., gg_par) are available
      if (!"ggplot2" %in% loadedNamespaces()) {
        suppressPackageStartupMessages(library(ggplot2))
      }
      # Work around environments where gg_par is not found by legend helpers
      if (!exists("gg_par", envir = .GlobalEnv) &&
        exists("gg_par", envir = asNamespace("ggplot2"), inherits = FALSE)) {
        assign("gg_par",
          get("gg_par", envir = asNamespace("ggplot2"), inherits = FALSE),
          envir = .GlobalEnv
        )
      }
    }
  }
  # Helper functions are part of the myR package; no additional sourcing required.

  # --- 1. Run multiple DEG methods ---
  res <- run_deg_consensus(
    sobj = sobj,
    contrast = contrast,
    methods = methods,
    cluster_id = cluster_id,
    sample_id = sample_id,
    group_id = group_id,
    batch_id = batch_id,
    remove_na_groups = TRUE,
    verbose = TRUE
  )
  if (save_intermediate) {
    qs::qsave(res, paste0(prefix, "_raw_results.qs"))
  }

  # --- 2. Standardise results ---
  std <- lapply(res$methods_run, function(m) {
    standardize_deg_results(res$results[[m]], m)
  })
  names(std) <- res$methods_run
  if (save_intermediate) {
    qs::qsave(std, paste0(prefix, "_standardized_results.qs"))
  }

  # --- 3. Build gene × method matrices ---
  # build_deg_matrices() requires a valid fdr_threshold even if significance_mode="pvalue"
  fdr_for_build <- if (is.null(fdr_threshold)) 0.1 else fdr_threshold

  matrices <- build_deg_matrices(
    standardized_results_list = std,
    genes = NULL,
    fdr_threshold = fdr_for_build,
    significance_mode = significance_mode,
    pvalue_threshold = pvalue_threshold
  )
  if (save_intermediate) {
    qs::qsave(matrices, paste0(prefix, "_deg_matrices.qs"))
  }

  # --- 4. Consensus analysis ---
  agreement_scores <- compute_agreement_scores(matrices$significance)
  consensus_scores <- compute_consensus_scores(matrices, agreement_scores)
  if (save_intermediate) {
    qs::qsave(agreement_scores, paste0(prefix, "_agreement_scores.qs"))
    qs::qsave(consensus_scores, paste0(prefix, "_consensus_scores.qs"))
  }

  # Method-level PCA (optional)
  pca_result <- tryCatch(
    {
      perform_deg_pca(matrices)
    },
    error = function(e) {
      warning("perform_deg_pca() failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  # Loose filtering at consensus level by default (keep all genes unless thresholds are provided)
  consensus_deg_list <- generate_consensus_deg_list(
    consensus_scores = consensus_scores,
    fdr_threshold = fdr_threshold,
    agreement_threshold = agreement_threshold,
    min_methods = min_methods
  )
  if (save_intermediate) {
    qs::qsave(consensus_deg_list, paste0(prefix, "_consensus_deg_list.qs"))
  }

  # --- 5. Visualisation ---
  if (make_plots) {
    plot_consensus_meta_distributions(consensus_scores, file_prefix = prefix)
    plot_consensus_volcano(consensus_scores, file_prefix = prefix, fdr_threshold = fdr_threshold)
    plot_consensus_heatmap(matrices, consensus_scores, file_prefix = prefix)
    if (!is.null(pca_result)) {
      plot_method_pca(pca_result, file_prefix = prefix)
    }
    plot_gene_umap(matrices, consensus_scores, file_prefix = prefix)
  }

  final <- list(
    raw_results = res,
    standardized_results = std,
    deg_matrices = matrices,
    agreement_scores = agreement_scores,
    pca_result = pca_result,
    consensus_scores = consensus_scores,
    consensus_deg_list = consensus_deg_list,
    output_prefix = prefix
  )

  if (save_intermediate) {
    qs::qsave(final, paste0(prefix, "_final_result.qs"))
  }

  return(final)
}


# ----------------------------------------------------------------------------
# Plotting helpers
# ----------------------------------------------------------------------------

#' Plot meta-analysis p-value distributions and effect size relationships
#'
#' @param consensus_scores Data frame from \code{compute_consensus_scores()}
#' @param file_prefix File prefix for saving plots (without extension)
#'
#' @export
plot_consensus_meta_distributions <- function(
  consensus_scores,
  file_prefix,
  caption = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping meta-distribution plots.", call. = FALSE)
    return(invisible(NULL))
  }

  cs <- consensus_scores
  if (!("meta_p" %in% colnames(cs))) {
    warning("meta_p column not found in consensus_scores; skipping meta-distribution plots.", call. = FALSE)
    return(invisible(NULL))
  }

  cs$neg_log10_meta_p <- -log10(pmax(cs$meta_p, 1e-300))

  # Capping 5% outliers
  # Handle cases with few points
  valid_p <- cs$neg_log10_meta_p[is.finite(cs$neg_log10_meta_p)]
  cap_threshold <- if (length(valid_p) > 0) quantile(valid_p, 0.95, na.rm = TRUE) else Inf

  # If the distribution is not skewed enough (e.g. all non-signif), capping might be weird.
  # But assume standard DEG output.
  cs$is_capped <- cs$neg_log10_meta_p > cap_threshold
  cs$neg_log10_meta_p_capped <- ifelse(cs$is_capped, cap_threshold, cs$neg_log10_meta_p)

  cs$n_significant_factor <- factor(cs$n_significant)

  df <- cs[!is.na(cs$neg_log10_meta_p_capped) & !is.na(cs$n_significant_factor), ]
  if (nrow(df) == 0) {
    warning("No valid rows for meta-distribution plots.", call. = FALSE)
    return(invisible(NULL))
  }

  # 1) meta p-value distribution by n_significant (box + violin)
  # Capping for Visualization Stability (e.g. 95% or 99%)
  # Only cap the plotting variable, not the original data for stats if needed (though stats on capped usually better for visuals)
  pvals <- df$neg_log10_meta_p[is.finite(df$neg_log10_meta_p)]
  cap_val <- quantile(pvals, 0.95, na.rm = TRUE)
  df$neg_log10_meta_p_capped <- ifelse(df$neg_log10_meta_p > cap_val, cap_val, df$neg_log10_meta_p)

  # 1) Meta p-value distribution vs n_methods
  # Violin + Boxplot overlay
  # Use CAPPED values for visualization to prevent flattening
  p_violin <- ggplot2::ggplot(df, ggplot2::aes(x = factor(n_significant), y = neg_log10_meta_p_capped)) +
    ggplot2::geom_violin(scale = "width", trim = FALSE, fill = "grey90", colour = "grey50") +
    ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
    ggplot2::labs(
      x = "Number of significant methods",
      y = "-log10(meta p-value) [Capped 95%]",
      title = "Distribution of meta p-values by consensus strength",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_meta_p_by_n_methods_violin.png"),
    plot = p_violin,
    width = 6, height = 5, dpi = 300
  )

  # 1-B) Ridge Plot implementation for Meta P-values
  if (requireNamespace("ggridges", quietly = TRUE)) {
    stats_meta <- do.call(rbind, lapply(split(df, df$n_significant), function(sub) {
      val <- sub$neg_log10_meta_p_capped # Use capped for stats label consistency? Or raw?
      # User wants "mean, sd" displayed implies descriptive stats of distribution shown.
      m_val <- mean(val, na.rm = TRUE)
      data.frame(
        n_significant = unique(sub$n_significant),
        label = sprintf("Mean: %.2f", m_val),
        stringsAsFactors = FALSE
      )
    }))
    stats_meta$n_significant <- factor(stats_meta$n_significant)

    p_ridge_meta <- ggplot2::ggplot(df, ggplot2::aes(x = neg_log10_meta_p_capped, y = factor(n_significant), fill = factor(n_significant))) +
      ggridges::geom_density_ridges(alpha = 0.6) +
      ggplot2::geom_text(
        data = stats_meta,
        ggplot2::aes(x = Inf, y = n_significant, label = label),
        inherit.aes = FALSE,
        hjust = 1.1, size = 3, vjust = 0.5
      ) +
      ggplot2::labs(
        x = "-log10(meta p-value) [Capped 95%]",
        y = "Number of significant methods",
        title = "Ridge distribution of meta p-values",
        caption = caption
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", plot.caption = ggplot2::element_text(hjust = 0))

    ggplot2::ggsave(
      filename = paste0(file_prefix, "_meta_p_ridges.png"),
      plot = p_ridge_meta,
      width = 6, height = 6, dpi = 300
    )
  }

  # 2) Relationship between effect size and meta p-value
  # Color needs to handle outliers
  df$ColorGroup <- as.character(df$n_significant)
  df$ColorGroup[df$is_capped] <- "Capped (>95%)"

  # Ensure "Capped" is last level or distinct
  # Sort levels numerically + Capped
  lvls <- sort(unique(df$n_significant))
  df$ColorGroup <- factor(df$ColorGroup, levels = c(as.character(lvls), "Capped (>95%)"))

  # Define colors: Gradient for numbers, Red/Grey for Capped
  n_colors <- length(lvls)
  # viridis or blues
  base_pal <- grDevices::colorRampPalette(c("#deebf7", "#3182bd", "#08306b"))(n_colors)
  # Assign specific color for Capped (e.g., #d73027 Red or #525252 Dark Grey)
  # User asked for "separate color".
  pal_map <- stats::setNames(c(base_pal, "#d95f02"), levels(df$ColorGroup)) # Orange/Red for capped

  p_scatter <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = weighted_beta,
      y = neg_log10_meta_p_capped,
      colour = ColorGroup
      # shape = is_capped ?
    )
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_colour_manual(values = pal_map) +
    ggplot2::labs(
      x = "Weighted log2 fold-change (weighted_beta)",
      y = "-log10(meta p-value) [Capped]",
      colour = "n_sig / Capped",
      title = "Effect size vs meta-analysis p-value",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right") # Enable legend

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_effect_vs_meta_p_scatter.png"),
    plot = p_scatter,
    width = 8, height = 6, dpi = 300
  )

  # 3) Log-Scale X Axis version
  p_scatter_log <- p_scatter +
    ggplot2::scale_x_continuous(
      trans = "pseudo_log",
      breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
    ) +
    ggplot2::labs(title = "Effect size (Pseudo-Log Scale) vs meta-analysis p-value", caption = caption)

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_effect_vs_meta_p_scatter_logX.png"),
    plot = p_scatter_log,
    width = 8, height = 6, dpi = 300
  )

  invisible(list(violin = p_violin, scatter = p_scatter))
}


#' Volcano plot for consensus scores
#'
#' @param consensus_scores Data frame from \code{compute_consensus_scores()}
#' @param file_prefix File prefix for saving plot
#' @param fdr_threshold Optional FDR threshold on \code{meta_p_adj} to highlight
#'   significant genes (default: NULL, no highlighting by FDR)
#' @param top_n_label Number of top genes (by consensus_score) to label
#'
#' @export
plot_consensus_volcano <- function(
  consensus_scores,
  file_prefix,
  fdr_threshold = NULL,
  top_n_label = 30,
  caption = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping volcano plot.", call. = FALSE)
    return(invisible(NULL))
  }

  cs <- consensus_scores
  if (!("meta_p" %in% colnames(cs))) {
    warning("meta_p column not found in consensus_scores; skipping volcano plot.", call. = FALSE)
    return(invisible(NULL))
  }

  cs$neg_log10_meta_p <- -log10(pmax(cs$meta_p, 1e-300))

  # Determine significance for coloring
  # If fdr_threshold is provided, use meta_p_adj. Otherwise, use meta_p < 0.05 (default hline).
  cs$significant <- FALSE
  if (!is.null(fdr_threshold) && "meta_p_adj" %in% colnames(cs)) {
    cs$significant <- !is.na(cs$meta_p_adj) & cs$meta_p_adj <= fdr_threshold
  } else {
    # If no FDR threshold, use a p-value threshold of 0.05 for highlighting, matching the hline
    cs$significant <- !is.na(cs$meta_p) & cs$meta_p <= 0.05
  }

  # Select genes to label
  cs_ordered <- cs[order(-cs$consensus_score), ]
  label_genes <- head(cs_ordered$gene[!is.na(cs_ordered$gene)], top_n_label)
  cs$label <- ifelse(cs$gene %in% label_genes, cs$gene, NA_character_)

  # Capping 5% outliers for Volcano
  valid_p <- cs$neg_log10_meta_p[is.finite(cs$neg_log10_meta_p)]
  cap_threshold <- if (length(valid_p) > 0) quantile(valid_p, 0.95, na.rm = TRUE) else Inf

  cs$is_capped <- cs$neg_log10_meta_p > cap_threshold
  cs$neg_log10_meta_p_capped <- ifelse(cs$is_capped, cap_threshold, cs$neg_log10_meta_p)

  # Color logic: Significant (Red), Capped (Orange), Other (Grey)
  cs$ColorGroup <- factor(
    ifelse(cs$is_capped, "Capped (>95%)",
      ifelse(cs$significant, "Significant", "NS")
    ),
    levels = c("Significant", "NS", "Capped (>95%)")
  )

  pal_map <- c("Significant" = "#d73027", "NS" = "#636363", "Capped (>95%)" = "#fdae61")

  p_volcano <- ggplot2::ggplot(
    cs,
    ggplot2::aes(
      x = weighted_beta,
      y = neg_log10_meta_p_capped,
      colour = ColorGroup
    )
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 0.6) +
    ggplot2::scale_colour_manual(values = pal_map) +
    # Reference lines
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", alpha = 0.5) +
    ggplot2::labs(
      x = "Weighted log2 fold-change (weighted_beta)",
      y = "-log10(meta p-value) [Capped]",
      colour = "Status",
      title = "Consensus volcano plot",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")

  # Add labels for top genes (if ggrepel is available)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_volcano <- p_volcano +
      ggrepel::geom_text_repel(
        data = cs[!is.na(cs$label), ],
        ggplot2::aes(label = label),
        size = 2.5,
        max.overlaps = Inf
      )
  } else {
    p_volcano <- p_volcano +
      ggplot2::geom_text(
        data = cs[!is.na(cs$label), ],
        ggplot2::aes(label = label),
        size = 2.2,
        vjust = -0.5
      )
  }

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_volcano.png"),
    plot = p_volcano,
    width = 7, height = 5, dpi = 300
  )

  # Variant: Log Scale X
  p_volcano_log <- p_volcano +
    ggplot2::scale_x_continuous(
      trans = "pseudo_log",
      breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
    ) +
    ggplot2::labs(title = "Consensus volcano plot (Log X)", caption = caption)

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_volcano_logX.png"),
    plot = p_volcano_log,
    width = 7, height = 5, dpi = 300
  )

  invisible(p_volcano)
}


#' Heatmap of top consensus genes across methods (beta matrix)
#'
#' @param deg_matrices List from \code{build_deg_matrices()}
#' @param consensus_scores Data frame from \code{compute_consensus_scores()}
#' @param file_prefix File prefix for saving plot
#' @param top_n Number of top genes (by consensus_score) to include
#'
#' @export
plot_consensus_heatmap <- function(
  deg_matrices,
  consensus_scores,
  file_prefix,
  method_top_n = 10,
  caption = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping heatmap.", call. = FALSE)
    return(invisible(NULL))
  }

  beta <- deg_matrices$beta
  cs <- consensus_scores

  # Select Top N genes from EACH method by |beta| and take UNION
  methods <- colnames(beta)
  selected_genes <- character(0)

  for (m in methods) {
    if (m %in% colnames(beta)) {
      # Sort by abs(beta) desc
      vals <- abs(beta[, m])
      # Handle NAs
      vals[is.na(vals)] <- -1

      top_m <- head(rownames(beta)[order(vals, decreasing = TRUE)], method_top_n)
      selected_genes <- c(selected_genes, top_m)
    }
  }
  top_genes <- unique(selected_genes)

  # Filter to those in consensus scores
  top_genes <- top_genes[top_genes %in% cs$gene]
  if (length(top_genes) == 0) {
    warning("No overlapping genes between consensus_scores and beta matrix; skipping heatmap.", call. = FALSE)
    return(invisible(NULL))
  }

  beta_sub <- beta[top_genes, , drop = FALSE]

  # Z-score standardization per Method (Column) to handle scale differences (e.g. limma-trend outliers)
  beta_sub_scaled <- scale(beta_sub)

  # Long-format for ggplot heatmap
  df_long <- data.frame(
    gene = rep(rownames(beta_sub_scaled), times = ncol(beta_sub_scaled)),
    method = rep(colnames(beta_sub_scaled), each = nrow(beta_sub_scaled)),
    beta = as.vector(beta_sub_scaled),
    stringsAsFactors = FALSE
  )

  p_hm <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = method, y = gene, fill = beta)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      name = "Z-Score(Beta)\n(Per Method)"
    ) +
    ggplot2::labs(
      x = "Method",
      y = "Gene (Union of Top Genes)",
      title = "Top genes across methods (Column Z-score scaled)",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 6),
      legend.position = "right",
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_top_genes_heatmap.png"),
    plot = p_hm,
    width = 7, height = 7, dpi = 300
  )

  invisible(p_hm)
}


#' PCA plot for DEG methods
#'
#' @param pca_result Result from \code{perform_deg_pca()}
#' @param file_prefix File prefix for saving plot
#'
#' @export
plot_method_pca <- function(
  pca_result,
  file_prefix,
  caption = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping method PCA plot.", call. = FALSE)
    return(invisible(NULL))
  }

  coords <- as.data.frame(pca_result$coordinates)
  if (ncol(coords) < 2) {
    warning("PCA coordinates have fewer than 2 components; skipping PCA plot.", call. = FALSE)
    return(invisible(NULL))
  }
  coords$method <- rownames(coords)
  colnames(coords)[1:2] <- c("PC1", "PC2")

  p_pca <- ggplot2::ggplot(
    coords,
    ggplot2::aes(x = PC1, y = PC2, label = method)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "PC1",
      y = "PC2",
      title = "PCA of DEG methods",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_pca <- p_pca +
      ggrepel::geom_text_repel(size = 3)
  } else {
    p_pca <- p_pca +
      ggplot2::geom_text(vjust = -0.5, size = 3)
  }

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_methods_pca.png"),
    plot = p_pca,
    width = 6, height = 5, dpi = 300
  )

  invisible(p_pca)
}


#' UMAP of genes based on method-level beta profiles
#'
#' @param deg_matrices List from \code{build_deg_matrices()}
#' @param consensus_scores Data frame from \code{compute_consensus_scores()}
#' @param file_prefix File prefix for saving plot
#' @param top_n Number of genes (by consensus_score) to include (default: 1000)
#' @param min_methods Minimum number of methods with non-NA beta to keep a gene
#'
#' @export
plot_gene_umap <- function(
  deg_matrices,
  consensus_scores,
  file_prefix,
  top_n = 1000,
  min_methods = 1,
  caption = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping gene UMAP plot.", call. = FALSE)
    return(invisible(NULL))
  }
  if (!requireNamespace("uwot", quietly = TRUE)) {
    warning("uwot not available; skipping gene UMAP plot.", call. = FALSE)
    return(invisible(NULL))
  }

  beta <- deg_matrices$beta
  cs <- consensus_scores

  # Select genes with enough non-NA betas
  n_non_na <- rowSums(!is.na(beta))
  keep_genes <- names(n_non_na)[n_non_na >= min_methods]
  if (length(keep_genes) == 0) {
    warning("No genes with sufficient non-NA betas; skipping gene UMAP.", call. = FALSE)
    return(invisible(NULL))
  }

  # Select top genes by consensus_score among the kept genes
  cs_ordered <- cs[order(-cs$consensus_score), ]
  top_genes <- cs_ordered$gene[cs_ordered$gene %in% keep_genes]
  top_genes <- head(unique(stats::na.omit(top_genes)), top_n)

  if (length(top_genes) < 5) {
    warning("Fewer than 5 genes available for UMAP; skipping.", call. = FALSE)
    return(invisible(NULL))
  }

  mat <- beta[top_genes, , drop = FALSE]
  # Replace NA with 0 for UMAP
  mat[is.na(mat)] <- 0

  set.seed(1)
  emb <- uwot::umap(
    mat,
    n_neighbors = 15,
    min_dist = 0.3,
    metric = "cosine"
  )

  # Gene Clustering on UMAP coordinates
  set.seed(42)
  km <- kmeans(emb, centers = min(8, nrow(emb)))
  gene_clusters <- factor(km$cluster)

  df_umap <- data.frame(
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    gene = top_genes,
    consensus_score = cs$consensus_score[match(top_genes, cs$gene)],
    n_significant = cs$n_significant[match(top_genes, cs$gene)],
    GeneCluster = gene_clusters,
    stringsAsFactors = FALSE
  )
  df_umap$n_significant_factor <- factor(df_umap$n_significant)

  # Plot 1: Standard
  p_umap <- ggplot2::ggplot(
    df_umap,
    ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      colour = GeneCluster, # Group by GeneCluster as requested
      shape = n_significant_factor
    )
  ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::labs(
      x = "UMAP1",
      y = "UMAP2",
      colour = "Gene Cluster (K=8)",
      shape = "n_significant",
      title = "UMAP of genes: Grouped by Profile Similarity",
      caption = caption
    ) +
    ggplot2::scale_colour_manual(values = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
      "#A65628", "#F781BF", "#999999", "#66C2A5"
    ), name = "Gene Cluster (K=8)") +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12, 13, 14)) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 3), ncol = 1, order = 1),
      shape = ggplot2::guide_legend(override.aes = list(size = 3), ncol = 1, order = 2)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      plot.caption = ggplot2::element_text(hjust = 0),
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm") # Top margin to avoid legend cutoff
    )

  # New: Save UMAP Clusters to disk
  cluster_file_base <- paste0(file_prefix, "_gene_clusters")

  # Save simple CSV
  utils::write.csv(
    df_umap[, c("gene", "GeneCluster", "consensus_score", "n_significant")],
    paste0(cluster_file_base, ".csv"),
    row.names = FALSE
  )

  # Save QS object
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(df_umap, paste0(cluster_file_base, ".qs"))
  }

  ggplot2::ggsave(
    filename = paste0(file_prefix, "_genes_umap_clusters.png"),
    plot = p_umap,
    width = 8, height = 6, dpi = 300
  )

  invisible(p_umap)
}


#' Export Method Summary CSV
#' @export
export_method_summary_csv <- function(res, file_path) {
  if (is.null(res$raw_results)) {
    return()
  }

  # Extract info
  methods <- res$raw_results$methods_run

  # Try to find formula/args from attributes of standardized results or raw results
  info_list <- lapply(methods, function(m) {
    raw <- res$raw_results$results[[m]]

    # Attribute extraction (dependent on method implementation)
    fmla <- attr(raw, "formula")
    if (is.null(fmla)) fmla <- "Internal/Default"

    # Args?
    args <- attr(raw, "run_info")
    arg_str <- if (!is.null(args)) paste(names(args), as.character(args), sep = "=", collapse = "; ") else ""

    # Test Stat?
    # Check first row
    stat_col <- "Unknown"
    if (is.data.frame(raw) && "statistic" %in% colnames(raw)) {
      stat_col <- "statistic column"
    } else if (grepl("LRT", m)) {
      stat_col <- "LRT Statistic (Signed)"
    } else if (grepl("Wald", m)) {
      stat_col <- "Wald Statistic"
    } else {
      stat_col <- "t-statistic / F-statistic"
    }

    data.frame(
      Method = m,
      Formula = as.character(fmla),
      Arguments = arg_str,
      TestStatistic = stat_col,
      stringsAsFactors = FALSE
    )
  })

  df_sum <- do.call(rbind, info_list)
  utils::write.csv(df_sum, file_path, row.names = FALSE)
}

#' Plot P-value Ridge Plots
#' @export
plot_pvalue_ridges <- function(res, file_prefix, caption = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggridges", quietly = TRUE)) {
    return(invisible(NULL))
  }

  std <- res$standardized_results
  if (is.null(std)) {
    return()
  }

  df_list <- lapply(names(std), function(m) {
    d <- std[[m]]
    # Check standard columns
    pval <- if ("pvalue" %in% names(d)) d$pvalue else NULL
    fdr <- if ("pvalue_adj" %in% names(d)) d$pvalue_adj else pval

    if (is.null(pval)) {
      return(NULL)
    }

    data.frame(
      pvalue = pval,
      FDR = fdr,
      method = m,
      stringsAsFactors = FALSE
    )
  })

  df_all <- do.call(rbind, df_list)
  if (is.null(df_all) || nrow(df_all) == 0) {
    return()
  }

  # 1. Raw P-value
  # Use -log10? Or just distribution? Ridge plots usually nice for 0-1 uniform-ish check or skew.
  # Let's do raw 0-1 for pvalue.
  # Calculate Statistics for Labels
  stats_df <- do.call(rbind, lapply(split(df_all, df_all$method), function(sub) {
    m_val <- mean(sub$pvalue, na.rm = TRUE)
    sd_val <- sd(sub$pvalue, na.rm = TRUE)
    prop_sig <- mean(sub$pvalue < 0.05, na.rm = TRUE) * 100
    data.frame(
      method = unique(sub$method),
      label = sprintf("Mean: %.2f\nSD: %.2f\n<0.05: %.1f%%", m_val, sd_val, prop_sig),
      stringsAsFactors = FALSE
    )
  }))

  # Merge stats into labels
  # Or just use text geom

  max_y <- length(unique(df_all$method))

  p1 <- ggplot2::ggplot(df_all, ggplot2::aes(x = pvalue, y = method, fill = method)) +
    ggridges::geom_density_ridges(alpha = 0.6, scale = 1.5) +
    ggplot2::geom_text(
      data = stats_df,
      ggplot2::aes(x = 1.1, y = method, label = label),
      inherit.aes = FALSE,
      hjust = 0, size = 3, vjust = 0.5
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1.9), breaks = seq(0, 1, 0.2)) + # Wider space for text
    ggplot2::labs(
      title = "P-value Distribution by Method",
      subtitle = "Stats: Mean, SD, % Significant (<0.05)",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      plot.caption = ggplot2::element_text(hjust = 0),
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )

  ggplot2::ggsave(paste0(file_prefix, "_pvalue_ridges.png"), p1, width = 10, height = 8)

  # 2. FDR distribution
  # Similar stats for FDR?
  stats_fdr <- do.call(rbind, lapply(split(df_all, df_all$method), function(sub) {
    prop_sig <- mean(sub$FDR < 0.05, na.rm = TRUE) * 100
    data.frame(
      method = unique(sub$method),
      label = sprintf("FDR<0.05: %.1f%%", prop_sig),
      stringsAsFactors = FALSE
    )
  }))

  p2 <- ggplot2::ggplot(df_all, ggplot2::aes(x = FDR, y = method, fill = method)) +
    ggridges::geom_density_ridges(alpha = 0.6) +
    ggplot2::geom_text(
      data = stats_fdr,
      ggplot2::aes(x = 1.05, y = method, label = label),
      inherit.aes = FALSE,
      hjust = 0, size = 3
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1.3)) +
    ggplot2::labs(title = "FDR Distribution by Method", caption = caption) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  ggplot2::ggsave(paste0(file_prefix, "_fdr_ridges.png"), p2, width = 8, height = 6)
}

#' Generate Summary Table Image
#' @export
generate_summary_table_image <- function(res, file_path) {
  if (!requireNamespace("gridExtra", quietly = TRUE) || !requireNamespace("grid", quietly = TRUE)) {
    return()
  }

  # Generate DF using existing helper
  # Generate DF using existing helper
  temp_csv <- tempfile(fileext = ".csv")
  export_method_summary_csv(res, temp_csv)

  if (!file.exists(temp_csv)) {
    return()
  }
  df <- utils::read.csv(temp_csv)

  # Clean up text for plot
  # Truncate very long arguments
  df$Arguments <- substr(df$Arguments, 1, 50)

  # Create Table Grob
  t_grob <- gridExtra::tableGrob(df, rows = NULL, theme = gridExtra::ttheme_default(base_size = 9))

  h <- 0.5 * nrow(df) + 2
  w <- 12

  ggplot2::ggsave(file_path, t_grob, width = w, height = h)
}

#' Export Method Formulas Text
#' @export
export_method_formulas_text <- function(res, file_path) {
  if (is.null(res$raw_results)) {
    return()
  }
  methods <- res$raw_results$methods_run

  lines <- c("Method Formulas Summary", "=======================", "")

  for (m in methods) {
    raw <- res$raw_results$results[[m]]
    fmla <- attr(raw, "formula")
    if (is.null(fmla)) fmla <- "Internal/Default"

    # Args?
    args <- attr(raw, "run_info")

    # Construct line
    line <- paste0(m, ": ", as.character(fmla))
    lines <- c(lines, line)
  }

  writeLines(lines, file_path)
}
