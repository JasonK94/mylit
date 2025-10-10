#' Differential Expression Analysis
#'
#' This module provides comprehensive differential expression analysis methods
#' including pseudobulk (edgeR), linear regression, and linear mixed models (LMM).
#'
#' @name differential_expression
NULL

# =============================================================================
# PSEUDOBULK DIFFERENTIAL EXPRESSION (edgeR)
# =============================================================================

#' Prepare Pseudobulk Data for edgeR
#'
#' Aggregates single-cell counts to pseudo-bulk level for differential expression analysis.
#'
#' @param object Seurat object
#' @param sample_col Metadata column identifying samples/replicates
#' @param group_col Metadata column for grouping (e.g., condition, treatment)
#' @param cluster_col Optional cluster column for per-cluster analysis (default: NULL)
#' @param target_cluster If cluster_col is specified, analyze only this cluster (default: NULL)
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#' @param min_cells Minimum cells per sample to include (default: 10)
#' @param min_counts Minimum total counts per gene (default: 10)
#'
#' @return List containing counts, metadata, and cluster info
#'
#' @export
prepare_pseudobulk_edgeR <- function(object, 
                                     sample_col, 
                                     group_col,
                                     cluster_col = NULL,
                                     target_cluster = NULL,
                                     assay = "RNA",
                                     slot = "counts",
                                     min_cells = 10,
                                     min_counts = 10) {
  
  # Validate inputs
  required_cols <- c(sample_col, group_col)
  if (!is.null(cluster_col)) {
    required_cols <- c(required_cols, cluster_col)
  }
  
  missing_cols <- setdiff(required_cols, colnames(object@meta.data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  # Subset to target cluster if specified
  if (!is.null(cluster_col) && !is.null(target_cluster)) {
    cells_to_use <- colnames(object)[object@meta.data[[cluster_col]] == target_cluster]
    
    if (length(cells_to_use) == 0) {
      stop("No cells found for cluster: ", target_cluster)
    }
    
    object <- object[, cells_to_use]
    message("Analyzing cluster ", target_cluster, " (", length(cells_to_use), " cells)")
  }
  
  # Get counts
  counts <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  
  # Get metadata
  metadata <- object@meta.data[, required_cols, drop = FALSE]
  
  # Create sample identifiers
  if (!is.null(cluster_col) && is.null(target_cluster)) {
    metadata$pb_sample <- paste0(metadata[[sample_col]], "_", metadata[[cluster_col]])
  } else {
    metadata$pb_sample <- as.character(metadata[[sample_col]])
  }
  
  # Check for minimum cells per sample
  cells_per_sample <- table(metadata$pb_sample)
  valid_samples <- names(cells_per_sample)[cells_per_sample >= min_cells]
  
  if (length(valid_samples) == 0) {
    stop("No samples have >= ", min_cells, " cells")
  }
  
  if (length(valid_samples) < length(cells_per_sample)) {
    n_removed <- length(cells_per_sample) - length(valid_samples)
    message("Removing ", n_removed, " samples with < ", min_cells, " cells")
    
    keep_cells <- metadata$pb_sample %in% valid_samples
    counts <- counts[, keep_cells, drop = FALSE]
    metadata <- metadata[keep_cells, , drop = FALSE]
  }
  
  # Aggregate counts by sample
  unique_samples <- unique(metadata$pb_sample)
  pb_counts <- matrix(
    0, 
    nrow = nrow(counts), 
    ncol = length(unique_samples),
    dimnames = list(rownames(counts), unique_samples)
  )
  
  for (samp in unique_samples) {
    cells_in_sample <- metadata$pb_sample == samp
    pb_counts[, samp] <- Matrix::rowSums(counts[, cells_in_sample, drop = FALSE])
  }
  
  # Filter lowly expressed genes
  keep_genes <- Matrix::rowSums(pb_counts) >= min_counts
  pb_counts <- pb_counts[keep_genes, , drop = FALSE]
  
  message("Pseudobulk matrix: ", nrow(pb_counts), " genes x ", ncol(pb_counts), " samples")
  
  # Create sample-level metadata
  sample_metadata <- metadata[!duplicated(metadata$pb_sample), , drop = FALSE]
  rownames(sample_metadata) <- sample_metadata$pb_sample
  sample_metadata <- sample_metadata[colnames(pb_counts), , drop = FALSE]
  
  return(list(
    counts = pb_counts,
    metadata = sample_metadata,
    cluster = target_cluster
  ))
}

#' Run Pseudobulk Differential Expression with edgeR
#'
#' Performs edgeR-based pseudo-bulk differential gene expression analysis.
#'
#' @param object Seurat object or prepared pseudobulk list
#' @param sample_col Sample identifier column
#' @param group_col Group comparison column
#' @param comparison Two-element vector for comparison (e.g., c("Treatment", "Control"))
#' @param cluster_col Optional cluster column
#' @param target_cluster Specific cluster to analyze
#' @param mode Analysis mode: "overall", "per_cluster", "specific_cluster"
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#' @param min_cells Minimum cells per sample (default: 10)
#' @param min_counts Minimum total counts per gene (default: 10)
#' @param fdr_threshold FDR threshold (default: 0.05)
#' @param logfc_threshold Log fold change threshold (default: 0)
#'
#' @return Data frame with DE results or list of data frames if per_cluster
#'
#' @export
run_pseudobulk_deg <- function(object,
                               sample_col = NULL,
                               group_col = NULL,
                               comparison = NULL,
                               cluster_col = NULL,
                               target_cluster = NULL,
                               mode = c("overall", "per_cluster", "specific_cluster"),
                               assay = "RNA",
                               slot = "counts",
                               min_cells = 10,
                               min_counts = 10,
                               fdr_threshold = 0.05,
                               logfc_threshold = 0) {
  
  mode <- match.arg(mode)
  
  # Check if already prepared
  is_prepared <- is.list(object) && all(c("counts", "metadata") %in% names(object))
  
  if (!is_prepared && (is.null(sample_col) || is.null(group_col))) {
    stop("sample_col and group_col required when object is Seurat")
  }
  
  if (is.null(comparison) || length(comparison) != 2) {
    stop("comparison must be length 2 (e.g., c('Treatment', 'Control'))")
  }
  
  # Dispatch based on mode
  if (mode == "overall") {
    if (is_prepared) {
      pb_data <- object
    } else {
      pb_data <- prepare_pseudobulk_edgeR(
        object, sample_col, group_col,
        cluster_col = NULL, target_cluster = NULL,
        assay = assay, slot = slot,
        min_cells = min_cells, min_counts = min_counts
      )
    }
    
    return(.run_edger_analysis(
      pb_data$counts, pb_data$metadata, group_col,
      comparison, fdr_threshold, logfc_threshold
    ))
    
  } else if (mode == "per_cluster") {
    if (is.null(cluster_col)) {
      stop("cluster_col required for per_cluster mode")
    }
    
    clusters <- unique(object@meta.data[[cluster_col]])
    message("Running pseudobulk DEG for ", length(clusters), " clusters")
    
    results_list <- list()
    for (clust in clusters) {
      message("\n--- Cluster ", clust, " ---")
      tryCatch({
        pb_data <- prepare_pseudobulk_edgeR(
          object, sample_col, group_col,
          cluster_col = cluster_col, target_cluster = clust,
          assay = assay, slot = slot,
          min_cells = min_cells, min_counts = min_counts
        )
        results_list[[as.character(clust)]] <- .run_edger_analysis(
          pb_data$counts, pb_data$metadata, group_col,
          comparison, fdr_threshold, logfc_threshold
        )
      }, error = function(e) {
        message("Failed for cluster ", clust, ": ", e$message)
      })
    }
    return(results_list)
    
  } else if (mode == "specific_cluster") {
    if (is.null(cluster_col) || is.null(target_cluster)) {
      stop("cluster_col and target_cluster required for specific_cluster mode")
    }
    
    if (is_prepared) {
      pb_data <- object
    } else {
      pb_data <- prepare_pseudobulk_edgeR(
        object, sample_col, group_col,
        cluster_col = cluster_col, target_cluster = target_cluster,
        assay = assay, slot = slot,
        min_cells = min_cells, min_counts = min_counts
      )
    }
    
    return(.run_edger_analysis(
      pb_data$counts, pb_data$metadata, group_col,
      comparison, fdr_threshold, logfc_threshold
    ))
  }
}

#' Internal edgeR Analysis
#' @keywords internal
.run_edger_analysis <- function(counts, metadata, group_col,
                                comparison, fdr_threshold, logfc_threshold) {
  
  groups <- metadata[[group_col]]
  
  if (!all(comparison %in% groups)) {
    stop("Comparison groups not found: ", 
         paste(setdiff(comparison, groups), collapse = ", "))
  }
  
  # Create DGEList
  dge <- edgeR::DGEList(counts = counts, group = groups)
  dge <- edgeR::calcNormFactors(dge)
  
  # Design matrix
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(factor(groups))
  
  # Estimate dispersion
  dge <- edgeR::estimateDisp(dge, design)
  
  # Fit model
  fit <- edgeR::glmQLFit(dge, design)
  
  # Contrast
  contrast_formula <- paste0(comparison[1], "-", comparison[2])
  contrast <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
  
  # Test
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  results <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
  
  # Format results
  results$gene <- rownames(results)
  results$significant <- (results$FDR < fdr_threshold) & 
                         (abs(results$logFC) > logfc_threshold)
  
  col_order <- c("gene", "logFC", "logCPM", "F", "PValue", "FDR", "significant")
  results <- results[, col_order]
  
  n_sig <- sum(results$significant)
  n_up <- sum(results$significant & results$logFC > 0)
  n_down <- sum(results$significant & results$logFC < 0)
  
  message("DEG results: ", n_sig, " significant (", n_up, " up, ", n_down, " down)")
  
  return(results)
}

# =============================================================================
# LINEAR REGRESSION
# =============================================================================

#' Linear Regression for Gene Expression
#'
#' Performs linear regression analysis for gene expression against a regressor
#' with support for continuous, categorical, and ordinal predictors.
#'
#' @param sobj Seurat object
#' @param layer Expression layer: "counts", "data", "scale.data"
#' @param features Features to test (default: "all")
#' @param regressor Regressor variable name in metadata
#' @param regressor.type Type: "continuous", "categorical", "ordinal"
#' @param reference.level Reference level for categorical
#' @param ordinal.method Method for ordinal: "linear", "polynomial", "spline"
#' @param link.function Link function: "linear", "poisson", "negative.binomial"
#' @param effect Effect type: "fixed", "random"
#' @param covariates Covariate column names
#' @param min.cells Minimum cells expressing gene (default: 10)
#' @param return.full Return full results including Seurat object
#'
#' @return Data frame with regression results
#'
#' @examples
#' \dontrun{
#' # Continuous regressor
#' results <- linear_seurat(sobj, regressor = "age", regressor.type = "continuous")
#' 
#' # Categorical with covariates
#' results <- linear_seurat(
#'   sobj,
#'   regressor = "cell_type",
#'   regressor.type = "categorical",
#'   covariates = c("batch", "sex")
#' )
#' }
#'
#' @export
linear_seurat <- function(sobj, 
                          layer = c("counts", "data", "scale.data"), 
                          features = "all", 
                          regressor = "val1", 
                          regressor.type = c("continuous", "categorical", "ordinal"),
                          reference.level = NULL,
                          ordinal.method = c("linear", "polynomial", "spline"),
                          link.function = "linear",
                          effect = c("fixed", "random"), 
                          covariates = NULL,
                          min.cells = 10,
                          return.full = FALSE,
                          ...) {
  
  if (!inherits(sobj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  layer <- match.arg(layer)
  effect <- match.arg(effect)
  regressor.type <- match.arg(regressor.type)
  ordinal.method <- match.arg(ordinal.method)
  
  if (!regressor %in% colnames(sobj@meta.data)) {
    stop("Regressor '", regressor, "' not found in metadata")
  }
  
  # Get expression data
  expr_data <- Seurat::GetAssayData(sobj, slot = layer)
  
  # Select features
  if (features[1] == "all") {
    features <- rownames(expr_data)
  } else {
    missing <- setdiff(features, rownames(expr_data))
    if (length(missing) > 0) {
      warning("Features not found: ", paste(missing, collapse = ", "))
      features <- intersect(features, rownames(expr_data))
    }
  }
  
  # Filter by minimum cells
  if (min.cells > 0) {
    feature_counts <- Matrix::rowSums(expr_data[features, ] > 0)
    features <- features[feature_counts >= min.cells]
  }
  
  if (length(features) == 0) {
    stop("No features passed filtering")
  }
  
  # Prepare metadata
  meta_data <- sobj@meta.data
  regressor_values <- meta_data[[regressor]]
  
  # Process regressor based on type
  if (regressor.type == "categorical") {
    regressor_values <- as.factor(regressor_values)
    if (!is.null(reference.level)) {
      if (!reference.level %in% levels(regressor_values)) {
        stop("Reference level '", reference.level, "' not found")
      }
      regressor_values <- relevel(regressor_values, ref = reference.level)
    }
    cat("Categorical regressor levels:", paste(levels(regressor_values), collapse = ", "), "\n")
    cat("Reference level:", levels(regressor_values)[1], "\n")
    
  } else if (regressor.type == "ordinal") {
    if (!is.ordered(regressor_values)) {
      if (is.factor(regressor_values)) {
        regressor_values <- ordered(regressor_values)
      } else {
        unique_vals <- sort(unique(regressor_values[!is.na(regressor_values)]))
        regressor_values <- ordered(regressor_values, levels = unique_vals)
      }
    }
    cat("Ordinal regressor levels:", paste(levels(regressor_values), collapse = " < "), "\n")
    
  } else if (regressor.type == "continuous") {
    regressor_values <- as.numeric(regressor_values)
    if (all(is.na(regressor_values))) {
      stop("Regressor cannot be converted to numeric")
    }
  }
  
  # Handle missing values
  if (any(is.na(regressor_values))) {
    warning("Missing values in regressor variable")
    valid_cells <- !is.na(regressor_values)
    expr_data <- expr_data[, valid_cells]
    meta_data <- meta_data[valid_cells, ]
    regressor_values <- regressor_values[valid_cells]
  }
  
  # Prepare covariate formula
  covariate_formula <- ""
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta_data))
    if (length(missing_cov) > 0) {
      warning("Covariates not found: ", paste(missing_cov, collapse = ", "))
      covariates <- intersect(covariates, colnames(meta_data))
    }
    if (length(covariates) > 0) {
      covariate_formula <- paste("+", paste(covariates, collapse = " + "))
    }
  }
  
  # Fit models for each gene
  cat("Fitting models for", length(features), "features...\n")
  
  fit_gene_model <- function(gene) {
    expression <- as.numeric(expr_data[gene, ])
    
    model_data <- data.frame(
      expression = expression,
      regressor = regressor_values
    )
    
    if (!is.null(covariates) && length(covariates) > 0) {
      for (cov in covariates) {
        model_data[[cov]] <- meta_data[[cov]]
      }
    }
    
    model_data <- model_data[complete.cases(model_data), ]
    
    if (nrow(model_data) < 10) {
      return(data.frame(
        gene = gene, estimate = NA, std.error = NA,
        statistic = NA, p.value = NA, n_cells = nrow(model_data),
        model_type = paste(effect, link.function, regressor.type, sep = "_")
      ))
    }
    
    # Build formula based on regressor type
    if (regressor.type == "continuous") {
      regressor_formula <- "regressor"
    } else if (regressor.type == "categorical") {
      regressor_formula <- "regressor"
    } else if (regressor.type == "ordinal") {
      if (ordinal.method == "linear") {
        model_data$regressor_numeric <- as.numeric(model_data$regressor)
        regressor_formula <- "regressor_numeric"
      } else if (ordinal.method == "polynomial") {
        regressor_formula <- "poly(as.numeric(regressor), degree = min(3, length(levels(regressor))-1))"
      } else if (ordinal.method == "spline") {
        if (requireNamespace("splines", quietly = TRUE)) {
          regressor_formula <- "splines::ns(as.numeric(regressor), df = min(3, length(levels(regressor))-1))"
        } else {
          model_data$regressor_numeric <- as.numeric(model_data$regressor)
          regressor_formula <- "regressor_numeric"
        }
      }
    }
    
    # Fit model
    tryCatch({
      if (effect == "fixed") {
        formula_str <- paste("expression ~", regressor_formula, covariate_formula)
        
        if (link.function == "linear") {
          model <- lm(as.formula(formula_str), data = model_data, ...)
        } else if (link.function == "poisson") {
          model <- glm(as.formula(formula_str), data = model_data, family = poisson(), ...)
        } else if (link.function == "negative.binomial") {
          model <- MASS::glm.nb(as.formula(formula_str), data = model_data, ...)
        }
        
        model_summary <- broom::tidy(model)
        
        # Extract relevant rows
        if (regressor.type == "categorical") {
          regressor_rows <- model_summary[grepl("^regressor", model_summary$term), ]
          if (nrow(regressor_rows) > 1) {
            model_anova <- anova(model)
            pval_col <- which(grepl("Pr\\(>F\\)", colnames(model_anova)))
            overall_p <- if (length(pval_col) > 0) model_anova[1, pval_col] else regressor_rows$p.value[1]
          } else {
            overall_p <- regressor_rows$p.value[1]
          }
          
          result <- data.frame(
            gene = gene,
            estimate = ifelse(nrow(regressor_rows) > 0, regressor_rows$estimate[1], NA),
            std.error = ifelse(nrow(regressor_rows) > 0, regressor_rows$std.error[1], NA),
            statistic = ifelse(nrow(regressor_rows) > 0, regressor_rows$statistic[1], NA),
            p.value = overall_p,
            n_cells = nrow(model_data),
            model_type = paste(effect, link.function, regressor.type, sep = "_")
          )
        } else {
          if (regressor.type == "ordinal" && ordinal.method == "linear") {
            regressor_rows <- model_summary[model_summary$term == "regressor_numeric", ]
          } else {
            regressor_rows <- model_summary[grepl("regressor", model_summary$term), ]
          }
          
          result <- data.frame(
            gene = gene,
            estimate = ifelse(nrow(regressor_rows) > 0, regressor_rows$estimate[1], NA),
            std.error = ifelse(nrow(regressor_rows) > 0, regressor_rows$std.error[1], NA),
            statistic = ifelse(nrow(regressor_rows) > 0, regressor_rows$statistic[1], NA),
            p.value = ifelse(nrow(regressor_rows) > 0, regressor_rows$p.value[1], NA),
            n_cells = nrow(model_data),
            model_type = paste(effect, link.function, regressor.type, sep = "_")
          )
        }
      } else if (effect == "random") {
        if (is.null(covariates) || length(covariates) == 0) {
          stop("Random effects require at least one grouping variable in covariates")
        }
        
        group_var <- covariates[1]
        formula_str <- paste("expression ~", regressor_formula, "+ (1|", group_var, ")")
        if (length(covariates) > 1) {
          formula_str <- paste(formula_str, "+", paste(covariates[-1], collapse = " + "))
        }
        
        model <- lme4::lmer(as.formula(formula_str), data = model_data, ...)
        model_summary <- broom.mixed::tidy(model, effects = "fixed")
        
        regressor_rows <- model_summary[grepl("regressor", model_summary$term), ]
        
        result <- data.frame(
          gene = gene,
          estimate = ifelse(nrow(regressor_rows) > 0, regressor_rows$estimate[1], NA),
          std.error = ifelse(nrow(regressor_rows) > 0, regressor_rows$std.error[1], NA),
          statistic = ifelse(nrow(regressor_rows) > 0, regressor_rows$statistic[1], NA),
          p.value = ifelse(nrow(regressor_rows) > 0, regressor_rows$p.value[1], NA),
          n_cells = nrow(model_data),
          model_type = paste(effect, link.function, regressor.type, sep = "_")
        )
      }
      
      return(result)
      
    }, error = function(e) {
      data.frame(
        gene = gene, estimate = NA, std.error = NA,
        statistic = NA, p.value = NA,
        n_cells = nrow(model_data),
        model_type = paste(effect, link.function, regressor.type, sep = "_"),
        error = as.character(e)
      )
    })
  }
  
  # Apply to all features
  if (requireNamespace("parallel", quietly = TRUE) && length(features) > 100) {
    results <- parallel::mclapply(features, fit_gene_model, 
                                  mc.cores = parallel::detectCores() - 1)
  } else {
    results <- lapply(features, fit_gene_model)
  }
  
  results_df <- do.call(rbind, results)
  
  # Adjust p-values
  valid_pvals <- !is.na(results_df$p.value)
  results_df$adj.p.value <- NA
  if (sum(valid_pvals) > 0) {
    results_df$adj.p.value[valid_pvals] <- p.adjust(results_df$p.value[valid_pvals], method = "BH")
  }
  
  results_df <- results_df[order(results_df$p.value, na.last = TRUE), ]
  
  attr(results_df, "analysis_info") <- list(
    regressor = regressor,
    regressor.type = regressor.type,
    layer = layer,
    effect = effect,
    link.function = link.function,
    covariates = covariates,
    n_features_tested = length(features),
    n_cells = ncol(expr_data)
  )
  
  if (return.full) {
    return(list(results = results_df, seurat_object = sobj, features_tested = features))
  } else {
    return(results_df)
  }
}

# NOTE: Linear Mixed Model (LMM) functions from test.R should be added here
# fit_lmm_single_gene() and run_lmm_multiple_genes() will be integrated in the next iteration

