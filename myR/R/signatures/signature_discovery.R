#' Gene Signature Discovery Functions
#'
#' This module provides functions for discovering gene signatures that best
#' distinguish between groups using various machine learning and statistical methods.
#'
#' @name signature_discovery
NULL

#' Find Gene Signatures that Best Separate Groups
#'
#' Discovers gene signatures using multiple methods including Random Forest,
#' LASSO, limma, NMF, Wilcoxon, GAM, and PCA.
#'
#' @param data Seurat object, count matrix, or data.frame
#' @param meta.data Optional metadata data.frame (required if data is not Seurat)
#' @param target_var Column name in metadata representing the target variable
#' @param target_group For numeric: quantile cutoff or list(low, high). For factor: levels to compare
#' @param method One of: "tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings"
#' @param n_features Number of top features to return (default: 50)
#' @param preprocess Whether to normalize/scale data (default: TRUE)
#' @param min_cells Minimum cells expressing gene (default: 10)
#' @param min_pct Minimum percentage of cells expressing gene (default: 0.01)
#' @param return_model Return full model object (default: FALSE)
#' @param seed Random seed for reproducibility (default: 42)
#' @param ... Additional method-specific parameters
#'
#' @return List containing:
#'   \item{genes}{Character vector of selected genes}
#'   \item{weights}{Named numeric vector of gene weights/importance}
#'   \item{scores}{Per-cell signature scores}
#'   \item{performance}{Separation metrics (AUC, accuracy, etc.)}
#'   \item{method}{Method used}
#'   \item{model}{Full model object (if return_model=TRUE)}
#'
#' @examples
#' \dontrun{
#' # Random Forest for binary classification
#' result <- find_gene_signature(
#'   seurat_obj,
#'   target_var = "cell_type",
#'   target_group = c("TypeA", "TypeB"),
#'   method = "tree_based"
#' )
#' 
#' # LASSO for continuous variable (top vs bottom quartile)
#' result <- find_gene_signature(
#'   seurat_obj,
#'   target_var = "pseudotime",
#'   target_group = 0.25,
#'   method = "lasso"
#' )
#' }
#'
#' @export
find_gene_signature <- function(data,
                                meta.data = NULL,
                                target_var,
                                target_group = NULL,
                                method = c("tree_based", "lasso", "limma",
                                          "nmf", "wilcoxon", "gam", "pca_loadings"),
                                n_features = 50,
                                preprocess = TRUE,
                                min_cells = 10,
                                min_pct = 0.01,
                                return_model = FALSE,
                                seed = 42,
                                ...) {
  
  set.seed(seed)
  method <- match.arg(method)
  
  # ============================================================================
  # 1. Input validation and data extraction
  # ============================================================================
  
  is_seurat <- inherits(data, "Seurat")
  
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required but not installed")
    }
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    # Extract normalized data or raw counts
    if ("data" %in% names(data@assays[[1]])) {
      expr_mat <- as.matrix(Seurat::GetAssayData(data, slot = "data"))
    } else {
      expr_mat <- as.matrix(Seurat::GetAssayData(data, slot = "counts"))
    }
  } else {
    if (is.null(meta.data)) {
      stop("meta.data must be provided when data is not a Seurat object")
    }
    expr_mat <- as.matrix(data)
    # Transpose if cells are rows
    if (nrow(expr_mat) == nrow(meta.data)) {
      expr_mat <- t(expr_mat)
    }
  }
  
  # Validate target_var exists
  if (!target_var %in% colnames(meta.data)) {
    stop(sprintf("target_var '%s' not found in metadata columns: %s",
                 target_var, paste(colnames(meta.data), collapse = ", ")))
  }
  
  # Align cells
  common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
  if (length(common_cells) == 0) {
    stop("No common cells found between expression matrix and metadata")
  }
  expr_mat <- expr_mat[, common_cells]
  meta.data <- meta.data[common_cells, ]
  
  # ============================================================================
  # 2. Process target variable
  # ============================================================================
  
  target_values <- meta.data[[target_var]]
  
  if (is.numeric(target_values)) {
    # Continuous variable - create binary groups
    if (is.null(target_group)) {
      target_group <- 0.25  # default: bottom 25% vs top 25%
    }
    
    if (is.list(target_group)) {
      low_cutoff <- quantile(target_values, target_group$low, na.rm = TRUE)
      high_cutoff <- quantile(target_values, target_group$high, na.rm = TRUE)
      group_labels <- ifelse(target_values <= low_cutoff, "Low",
                             ifelse(target_values >= high_cutoff, "High", NA))
    } else if (length(target_group) == 1 && target_group < 1) {
      low_cutoff <- quantile(target_values, target_group, na.rm = TRUE)
      high_cutoff <- quantile(target_values, 1 - target_group, na.rm = TRUE)
      group_labels <- ifelse(target_values <= low_cutoff, "Low",
                             ifelse(target_values >= high_cutoff, "High", NA))
    } else {
      group_labels <- ifelse(target_values < target_group, "Low", "High")
    }
    
    # Remove NAs (middle group)
    keep_cells <- !is.na(group_labels)
    expr_mat <- expr_mat[, keep_cells]
    meta.data <- meta.data[keep_cells, ]
    target_binary <- factor(group_labels[keep_cells])
    
  } else {
    # Categorical variable
    if (!is.null(target_group)) {
      keep_cells <- target_values %in% target_group
      expr_mat <- expr_mat[, keep_cells]
      meta.data <- meta.data[keep_cells, ]
      target_binary <- factor(target_values[keep_cells])
    } else {
      target_binary <- factor(target_values)
    }
  }
  
  if (length(unique(target_binary)) < 2) {
    stop("Target variable must have at least 2 groups after processing")
  }
  
  n_groups <- length(unique(target_binary))
  
  # ============================================================================
  # 3. Filter and preprocess genes
  # ============================================================================
  
  n_cells_expr <- rowSums(expr_mat > 0)
  pct_cells_expr <- n_cells_expr / ncol(expr_mat)
  keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  
  expr_mat <- expr_mat[keep_genes, ]
  
  if (nrow(expr_mat) == 0) {
    stop("No genes pass filtering criteria")
  }
  
  # Preprocessing
  if (preprocess) {
    # Log-normalize if data appears to be raw counts
    if (max(expr_mat) > 100) {
      expr_mat <- log1p(expr_mat)
    }
    # Scale genes (z-score) for certain methods
    if (method %in% c("lasso", "gam", "pca_loadings")) {
      gene_means <- rowMeans(expr_mat)
      gene_sds <- apply(expr_mat, 1, sd)
      gene_sds[gene_sds == 0] <- 1
      expr_mat <- (expr_mat - gene_means) / gene_sds
    }
  }
  
  # ============================================================================
  # 4. Dispatch to method-specific implementation
  # ============================================================================
  
  result <- switch(
    method,
    tree_based = .find_signature_rf(expr_mat, target_binary, n_features, n_groups, return_model, ...),
    lasso = .find_signature_lasso(expr_mat, target_binary, n_features, n_groups, return_model, ...),
    limma = .find_signature_limma(expr_mat, target_binary, n_features, n_groups, return_model, ...),
    wilcoxon = .find_signature_wilcoxon(expr_mat, target_binary, n_features, n_groups, ...),
    nmf = .find_signature_nmf(expr_mat, target_binary, n_features, n_groups, return_model, ...),
    gam = .find_signature_gam(expr_mat, target_binary, n_features, n_groups, ...),
    pca_loadings = .find_signature_pca(expr_mat, target_binary, n_features, n_groups, return_model, ...)
  )
  
  # ============================================================================
  # 5. Return results
  # ============================================================================
  
  result$method <- method
  result$target_var <- target_var
  result$n_groups <- n_groups
  result$n_cells <- ncol(expr_mat)
  
  class(result) <- c("gene_signature", "list")
  return(result)
}

#' Print Method for gene_signature Objects
#'
#' @param x gene_signature object
#' @param ... Additional arguments (unused)
#'
#' @export
print.gene_signature <- function(x, ...) {
  cat("Gene Signature Object\n")
  cat("====================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Target variable: %s (%d groups)\n", x$target_var, x$n_groups))
  cat(sprintf("Number of cells: %d\n", x$n_cells))
  cat(sprintf("Number of genes in signature: %d\n", length(x$genes)))
  cat(sprintf("\nTop 10 genes:\n"))
  print(head(data.frame(gene = x$genes, weight = x$weights[x$genes]), 10))
  
  if (!is.null(x$performance)) {
    cat("\nPerformance:\n")
    print(x$performance)
  }
}

# Internal helper functions for each method will be defined below
# These are not exported but used internally by find_gene_signature

#' @keywords internal
.find_signature_rf <- function(expr_mat, target_binary, n_features, n_groups, return_model, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest package required. Install with: install.packages('randomForest')")
  }
  
  X <- t(expr_mat)
  y <- target_binary
  
  # Subsample genes if too many
  if (nrow(expr_mat) > 2000) {
    gene_vars <- apply(expr_mat, 1, var)
    top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:2000])
    X <- X[, top_var_genes]
  }
  
  # Train Random Forest
  rf_model <- randomForest::randomForest(
    x = X, y = y,
    ntree = 500,
    importance = TRUE,
    ...
  )
  
  # Extract feature importance
  importance_scores <- randomForest::importance(rf_model)
  if (n_groups == 2) {
    weights <- importance_scores[, "MeanDecreaseGini"]
  } else {
    weights <- rowMeans(importance_scores[, grep("MeanDecreaseGini",
                                                 colnames(importance_scores))])
  }
  
  # Select top features
  top_genes <- names(sort(weights, decreasing = TRUE)[1:min(n_features, length(weights))])
  weights <- weights[top_genes]
  
  # Calculate signature scores
  scores <- as.numeric(X[, top_genes] %*% weights)
  names(scores) <- rownames(X)
  
  # Performance metrics
  pred <- rf_model$predicted
  if (n_groups == 2) {
    if (requireNamespace("pROC", quietly = TRUE)) {
      roc_obj <- pROC::roc(y, scores, quiet = TRUE)
      auc <- as.numeric(pROC::auc(roc_obj))
    } else {
      auc <- NA
    }
    acc <- mean(pred == y)
    perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
  } else {
    acc <- mean(pred == y)
    perf <- list(accuracy = acc, confusion = table(pred, y))
  }
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = if (return_model) rf_model else NULL
  )
}

#' @keywords internal
.find_signature_lasso <- function(expr_mat, target_binary, n_features, n_groups, return_model, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet package required. Install with: install.packages('glmnet')")
  }
  
  X <- t(expr_mat)
  y <- target_binary
  
  # Fit LASSO with cross-validation
  if (n_groups == 2) {
    cv_fit <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 1, ...)
  } else {
    cv_fit <- glmnet::cv.glmnet(X, y, family = "multinomial", alpha = 1, ...)
  }
  
  # Extract coefficients
  coefs <- coef(cv_fit, s = "lambda.1se")
  
  if (n_groups == 2) {
    weights <- as.numeric(coefs[-1])
    names(weights) <- rownames(coefs)[-1]
  } else {
    coef_list <- lapply(coefs, function(x) as.numeric(x[-1]))
    weights <- rowMeans(do.call(cbind, coef_list))
    names(weights) <- rownames(coefs[[1]])[-1]
  }
  
  # Select non-zero coefficients
  nonzero_genes <- names(weights)[weights != 0]
  if (length(nonzero_genes) == 0) {
    warning("No genes selected by LASSO. Returning top genes by coefficient magnitude.")
    top_genes <- names(sort(abs(weights), decreasing = TRUE)[1:min(n_features, length(weights))])
  } else {
    top_genes <- nonzero_genes[1:min(n_features, length(nonzero_genes))]
  }
  weights <- weights[top_genes]
  
  # Signature scores
  scores <- as.numeric(X[, top_genes] %*% weights)
  names(scores) <- rownames(X)
  
  # Performance
  pred_probs <- predict(cv_fit, newx = X, s = "lambda.1se", type = "response")
  if (n_groups == 2) {
    pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels = levels(y))
    if (requireNamespace("pROC", quietly = TRUE)) {
      roc_obj <- pROC::roc(y, as.numeric(pred_probs), quiet = TRUE)
      auc <- as.numeric(pROC::auc(roc_obj))
    } else {
      auc <- NA
    }
    acc <- mean(pred == y)
    perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
  } else {
    pred <- levels(y)[apply(pred_probs, 1, which.max)]
    acc <- mean(pred == y)
    perf <- list(accuracy = acc, confusion = table(pred, y))
  }
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = if (return_model) cv_fit else NULL
  )
}

#' @keywords internal
.find_signature_limma <- function(expr_mat, target_binary, n_features, n_groups, return_model, ...) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package required. Install with: BiocManager::install('limma')")
  }
  
  design <- model.matrix(~0 + target_binary)
  colnames(design) <- levels(target_binary)
  
  # Fit linear model
  fit <- limma::lmFit(expr_mat, design)
  
  # Define contrasts
  if (n_groups == 2) {
    contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep = "-")
  } else {
    contrast_str <- limma::makeContrasts(
      contrasts = combn(levels(target_binary), 2, function(x) paste(x, collapse = "-")),
      levels = design
    )
  }
  
  contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2)
  
  # Extract results
  top_table <- limma::topTable(fit2, number = Inf, sort.by = "B")
  
  # Weights as moderated t-statistics
  weights <- top_table$t
  names(weights) <- rownames(top_table)
  
  top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
  weights <- weights[top_genes]
  
  # Signature scores
  scores <- colSums(expr_mat[top_genes, ] * weights)
  
  # Performance
  score_cutoff <- median(scores)
  pred <- factor(ifelse(scores > score_cutoff, levels(target_binary)[2],
                       levels(target_binary)[1]), levels = levels(target_binary))
  acc <- mean(pred == target_binary)
  perf <- list(accuracy = acc, top_table = top_table[top_genes, ])
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = if (return_model) fit2 else NULL
  )
}

#' @keywords internal
.find_signature_wilcoxon <- function(expr_mat, target_binary, n_features, n_groups, ...) {
  pvals <- numeric(nrow(expr_mat))
  effect_sizes <- numeric(nrow(expr_mat))
  
  for (i in 1:nrow(expr_mat)) {
    if (n_groups == 2) {
      group1 <- expr_mat[i, target_binary == levels(target_binary)[1]]
      group2 <- expr_mat[i, target_binary == levels(target_binary)[2]]
      test <- wilcox.test(group1, group2)
      pvals[i] <- test$p.value
      effect_sizes[i] <- median(group2) - median(group1)
    } else {
      test <- kruskal.test(expr_mat[i, ] ~ target_binary)
      pvals[i] <- test$p.value
      effect_sizes[i] <- var(tapply(expr_mat[i, ], target_binary, median))
    }
  }
  
  names(pvals) <- rownames(expr_mat)
  names(effect_sizes) <- rownames(expr_mat)
  
  # Adjust p-values
  padj <- p.adjust(pvals, method = "BH")
  
  # Select by p-value and effect size
  ranks <- rank(pvals) + rank(-abs(effect_sizes))
  top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
  
  weights <- effect_sizes[top_genes]
  scores <- colSums(expr_mat[top_genes, ] * weights)
  
  perf <- list(
    pvalues = pvals[top_genes],
    padj = padj[top_genes],
    effect_sizes = effect_sizes[top_genes]
  )
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = NULL
  )
}

#' @keywords internal
.find_signature_nmf <- function(expr_mat, target_binary, n_features, n_groups, return_model, ...) {
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("NMF package required. Install with: install.packages('NMF')")
  }
  
  # Ensure non-negative data
  expr_mat_pos <- expr_mat - min(expr_mat) + 0.01
  
  # Run NMF
  rank <- min(n_groups + 2, 10)
  nmf_res <- NMF::nmf(expr_mat_pos, rank = rank, ...)
  
  W <- NMF::basis(nmf_res)  # gene loadings
  H <- NMF::coef(nmf_res)   # cell scores
  
  # Find component most associated with target
  component_cors <- numeric(rank)
  for (k in 1:rank) {
    if (n_groups == 2) {
      component_cors[k] <- abs(cor(H[k, ], as.numeric(target_binary)))
    } else {
      component_cors[k] <- summary(aov(H[k, ] ~ target_binary))[[1]][1, "F value"]
    }
  }
  
  best_component <- which.max(component_cors)
  weights <- W[, best_component]
  names(weights) <- rownames(expr_mat)
  
  top_genes <- names(sort(weights, decreasing = TRUE)[1:min(n_features, length(weights))])
  weights <- weights[top_genes]
  
  scores <- H[best_component, ]
  names(scores) <- colnames(expr_mat)
  
  perf <- list(
    component = best_component,
    correlation = component_cors[best_component]
  )
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = if (return_model) nmf_res else NULL
  )
}

#' @keywords internal
.find_signature_gam <- function(expr_mat, target_binary, n_features, n_groups, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("mgcv package required. Install with: install.packages('mgcv')")
  }
  
  deviance_explained <- numeric(nrow(expr_mat))
  
  for (i in 1:min(500, nrow(expr_mat))) {
    gene_expr <- expr_mat[i, ]
    if (n_groups == 2) {
      gam_fit <- mgcv::gam(as.numeric(target_binary) ~ s(gene_expr), family = "binomial")
    } else {
      gam_fit <- mgcv::gam(as.numeric(target_binary) ~ s(gene_expr))
    }
    deviance_explained[i] <- summary(gam_fit)$dev.expl
  }
  
  names(deviance_explained) <- rownames(expr_mat)[1:length(deviance_explained)]
  
  top_genes <- names(sort(deviance_explained, decreasing = TRUE)[1:min(n_features,
                                                                        sum(deviance_explained > 0))])
  weights <- deviance_explained[top_genes]
  
  scores <- colSums(expr_mat[top_genes, ] * weights)
  
  perf <- list(deviance_explained = deviance_explained[top_genes])
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = NULL
  )
}

#' @keywords internal
.find_signature_pca <- function(expr_mat, target_binary, n_features, n_groups, return_model, ...) {
  # Run PCA
  pca_res <- prcomp(t(expr_mat), center = FALSE, scale. = FALSE)
  
  # Find PC most correlated with target
  pc_cors <- numeric(min(50, ncol(pca_res$x)))
  for (k in 1:length(pc_cors)) {
    if (n_groups == 2) {
      pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
    } else {
      pc_cors[k] <- summary(aov(pca_res$x[, k] ~ target_binary))[[1]][1, "F value"]
    }
  }
  
  best_pc <- which.max(pc_cors)
  weights <- pca_res$rotation[, best_pc]
  
  top_genes <- names(sort(abs(weights), decreasing = TRUE)[1:min(n_features, length(weights))])
  weights <- weights[top_genes]
  
  scores <- pca_res$x[, best_pc]
  
  perf <- list(
    PC = best_pc,
    correlation = pc_cors[best_pc],
    variance_explained = summary(pca_res)$importance[2, best_pc]
  )
  
  list(
    genes = top_genes,
    weights = weights,
    scores = scores,
    performance = perf,
    model = if (return_model) pca_res else NULL
  )
}


