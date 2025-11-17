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
  sample_id  = "hos_no",
  group_id   = "g3",
  batch_id   = "GEM",
  significance_mode = c("pvalue", "fdr"),
  pvalue_threshold  = 0.05,
  fdr_threshold     = NULL,
  agreement_threshold = NULL,
  min_methods       = NULL,
  output_dir        = NULL,
  prefix            = NULL,
  save_intermediate = TRUE,
  make_plots        = TRUE
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
               envir = .GlobalEnv)
      }
    }
  }
  # Ensure core consensus functions are available (if this file is sourced standalone)
  if (!exists("run_deg_consensus")) {
    source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")
  }
  if (!exists("standardize_deg_results")) {
    source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
  }
  if (!exists("compute_agreement_scores") || !exists("compute_consensus_scores") ||
      !exists("perform_deg_pca") || !exists("generate_consensus_deg_list")) {
    source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
  }
  # Ensure method-specific implementations are available
  if (!exists("runLIMMA_voom_v1") || !exists("runLIMMA_trend_v1")) {
    source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
  }
  if (!exists("runEDGER_LRT_v1") || !exists("runEDGER_QLF_v1")) {
    source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
  }
  if (!exists("runDESEQ2_Wald_v1") || !exists("runDESEQ2_LRT_v1")) {
    source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
  }
  
  # --- 1. Run multiple DEG methods ---
  res <- run_deg_consensus(
    sobj = sobj,
    contrast = contrast,
    methods = methods,
    cluster_id = cluster_id,
    sample_id  = sample_id,
    group_id   = group_id,
    batch_id   = batch_id,
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
  pca_result <- tryCatch({
    perform_deg_pca(matrices)
  }, error = function(e) {
    warning("perform_deg_pca() failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
  
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
  file_prefix
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
  cs$n_significant_factor <- factor(cs$n_significant)
  
  df <- cs[!is.na(cs$neg_log10_meta_p) & !is.na(cs$n_significant_factor), ]
  if (nrow(df) == 0) {
    warning("No valid rows for meta-distribution plots.", call. = FALSE)
    return(invisible(NULL))
  }
  
  # 1) meta p-value distribution by n_significant (box + violin)
  p_violin <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = n_significant_factor, y = neg_log10_meta_p)
  ) +
    ggplot2::geom_violin(fill = "#3182bd", alpha = 0.4, colour = NA) +
    ggplot2::geom_boxplot(width = 0.2, outlier.size = 0.3, outlier.alpha = 0.5) +
    ggplot2::labs(
      x = "Number of methods with significance (n_significant)",
      y = "-log10(meta p-value)",
      title = "Meta-analysis p-value distribution by n_significant"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    filename = paste0(file_prefix, "_meta_p_by_n_methods_violin.png"),
    plot = p_violin,
    width = 7, height = 5, dpi = 300
  )
  
  # 2) Relationship between effect size and meta p-value
  df2 <- df
  df2$neg_log10_meta_p <- -log10(pmax(df2$meta_p, 1e-300))
  
  p_scatter <- ggplot2::ggplot(
    df2,
    ggplot2::aes(
      x = weighted_beta,
      y = neg_log10_meta_p,
      colour = n_significant_factor
    )
  ) +
    ggplot2::geom_point(alpha = 0.4, size = 0.6) +
    ggplot2::labs(
      x = "Weighted log2 fold-change (weighted_beta)",
      y = "-log10(meta p-value)",
      colour = "n_significant",
      title = "Effect size vs meta-analysis p-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  ggplot2::ggsave(
    filename = paste0(file_prefix, "_effect_vs_meta_p_scatter.png"),
    plot = p_scatter,
    width = 7, height = 5, dpi = 300
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
  top_n_label = 30
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
  cs$significant <- FALSE
  if (!is.null(fdr_threshold) && "meta_p_adj" %in% colnames(cs)) {
    cs$significant <- !is.na(cs$meta_p_adj) & cs$meta_p_adj <= fdr_threshold
  }
  
  # Select genes to label
  cs_ordered <- cs[order(-cs$consensus_score), ]
  label_genes <- head(cs_ordered$gene[!is.na(cs_ordered$gene)], top_n_label)
  cs$label <- ifelse(cs$gene %in% label_genes, cs$gene, NA_character_)
  
  p_volcano <- ggplot2::ggplot(
    cs,
    ggplot2::aes(
      x = weighted_beta,
      y = neg_log10_meta_p,
      colour = significant
    )
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 0.6) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#d73027", `FALSE` = "#636363"),
      breaks = c(TRUE, FALSE),
      labels = c("Significant", "Other")
    ) +
    ggplot2::labs(
      x = "Weighted log2 fold-change (weighted_beta)",
      y = "-log10(meta p-value)",
      colour = if (!is.null(fdr_threshold)) {
        paste0("meta FDR <=", fdr_threshold)
      } else {
        "Significant (meta FDR)"
      },
      title = "Consensus volcano plot"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
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
  top_n = 50
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping heatmap.", call. = FALSE)
    return(invisible(NULL))
  }
  
  beta <- deg_matrices$beta
  cs <- consensus_scores
  
  # Select top genes by consensus_score
  cs_ordered <- cs[order(-cs$consensus_score), ]
  top_genes <- head(cs_ordered$gene[cs_ordered$gene %in% rownames(beta)], top_n)
  top_genes <- unique(stats::na.omit(top_genes))
  if (length(top_genes) == 0) {
    warning("No overlapping genes between consensus_scores and beta matrix; skipping heatmap.", call. = FALSE)
    return(invisible(NULL))
  }
  
  beta_sub <- beta[top_genes, , drop = FALSE]
  
  # Long-format for ggplot heatmap
  df_long <- data.frame(
    gene = rep(rownames(beta_sub), times = ncol(beta_sub)),
    method = rep(colnames(beta_sub), each = nrow(beta_sub)),
    beta = as.vector(beta_sub),
    stringsAsFactors = FALSE
  )
  
  p_hm <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = method, y = gene, fill = beta)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027") +
    ggplot2::labs(
      x = "Method",
      y = "Gene (top consensus)",
      fill = "log2FC",
      title = "Top consensus genes across methods (beta matrix)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 6),
      legend.position = "none"
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
  file_prefix
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
      title = "PCA of DEG methods"
    ) +
    ggplot2::theme_bw()
  
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
  min_methods = 1
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
  
  df_umap <- data.frame(
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    gene  = top_genes,
    consensus_score = cs$consensus_score[match(top_genes, cs$gene)],
    n_significant   = cs$n_significant[match(top_genes, cs$gene)],
    stringsAsFactors = FALSE
  )
  df_umap$n_significant_factor <- factor(df_umap$n_significant)
  
  p_umap <- ggplot2::ggplot(
    df_umap,
    ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      colour = consensus_score,
      shape = n_significant_factor
    )
  ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::scale_colour_viridis_c(option = "C") +
    ggplot2::labs(
      x = "UMAP1",
      y = "UMAP2",
      colour = "Consensus score",
      shape = "n_significant",
      title = "UMAP of genes based on method-level beta profiles"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  ggplot2::ggsave(
    filename = paste0(file_prefix, "_genes_umap.png"),
    plot = p_umap,
    width = 7, height = 5, dpi = 300
  )
  
  invisible(p_umap)
}


