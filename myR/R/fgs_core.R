# ==============================================================================
# FGS Core Implementation
# ==============================================================================

fgs_preprocess_data_v5 <- function(data,
                                   meta.data,
                                   target_var,
                                   target_group,
                                   control_vars,
                                   test_n,
                                   preprocess,
                                   min_cells,
                                   min_pct,
                                   methods_requiring_scale,
                                   methods_requiring_correction) {
  # === 1. Input validation ===
  is_seurat <- inherits(data, "Seurat")

  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required but not installed")
    }
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    default_assay <- Seurat::DefaultAssay(data)
    # [OPTIMIZATION] Keep as sparse matrix initially to save memory
    expr_mat <- tryCatch(
      {
        Seurat::GetAssayData(data, assay = default_assay, layer = "data")
      },
      error = function(e) {
        tryCatch(
          {
            Seurat::GetAssayData(data, assay = default_assay, slot = "data")
          },
          error = function(e) {
            Seurat::GetAssayData(data, assay = default_assay, slot = "counts")
          }
        )
      }
    )
  } else {
    if (is.null(meta.data)) {
      stop("meta.data must be provided when data is not a Seurat object")
    }
    expr_mat <- data # Keep as is (could be sparse or dense)
  }

  # === 2. [V5] NA 및 팩터 레벨 완벽 제거 ===
  message("V5 Preprocessing: 1. Cleaning NA values...")

  if (is.null(control_vars)) {
    vars_to_check <- target_var
  } else {
    vars_to_check <- c(target_var, control_vars)
  }

  complete_cases_idx <- complete.cases(meta.data[, vars_to_check, drop = FALSE])

  n_removed <- sum(!complete_cases_idx)
  if (n_removed > 0) {
    warning(sprintf("V5 Preprocessing: Removed %d cells with NAs in target/control vars.", n_removed))
  }

  meta.data <- meta.data[complete_cases_idx, , drop = FALSE]

  if (is.null(rownames(meta.data))) {
    stop("meta.data must have rownames corresponding to cell IDs.")
  }
  if (is.null(colnames(expr_mat))) {
    stop("Expression matrix must have column names corresponding to cell IDs.")
  }

  common_cells <- base::intersect(colnames(expr_mat), rownames(meta.data))
  if (length(common_cells) == 0) {
    stop("No overlapping cells between expression data and metadata.")
  }
  meta.data <- meta.data[common_cells, , drop = FALSE]
  expr_mat <- expr_mat[, common_cells, drop = FALSE] # Sparse indexing is efficient

  target_values <- meta.data[[target_var]]
  if (is.numeric(target_values)) {
    if (is.null(target_group)) {
      target_group <- 0.25
    }
    if (is.list(target_group)) {
      low_cutoff <- stats::quantile(target_values, target_group$low, na.rm = TRUE)
      high_cutoff <- stats::quantile(target_values, target_group$high, na.rm = TRUE)
    } else if (length(target_group) == 1 && is.numeric(target_group) && target_group < 1) {
      low_cutoff <- stats::quantile(target_values, target_group, na.rm = TRUE)
      high_cutoff <- stats::quantile(target_values, 1 - target_group, na.rm = TRUE)
    } else {
      low_cutoff <- target_group
      high_cutoff <- target_group
    }
    group_labels <- ifelse(target_values <= low_cutoff, "Low",
      ifelse(target_values >= high_cutoff, "High", NA)
    )
    keep_cells <- !is.na(group_labels)
    if (!all(keep_cells)) {
      warning(sprintf("Discarding %d cells outside numeric target thresholds.", sum(!keep_cells)))
      expr_mat <- expr_mat[, keep_cells, drop = FALSE]
      meta.data <- meta.data[keep_cells, , drop = FALSE]
      group_labels <- group_labels[keep_cells]
    }
    target_binary <- factor(group_labels)
  } else {
    if (!is.null(target_group)) {
      keep_cells <- target_values %in% target_group
      if (!all(keep_cells)) {
        warning(sprintf("Discarding %d cells not in target groups.", sum(!keep_cells)))
        expr_mat <- expr_mat[, keep_cells, drop = FALSE]
        meta.data <- meta.data[keep_cells, , drop = FALSE]
        target_values <- target_values[keep_cells]
      }
    }
    target_binary <- droplevels(factor(target_values))
  }
  meta.data$target_binary_var <- target_binary

  if (!is.null(control_vars)) {
    for (cv in control_vars) {
      if (is.factor(meta.data[[cv]])) {
        meta.data[[cv]] <- droplevels(meta.data[[cv]])
      }
    }
  }

  if (length(unique(target_binary)) < 2) {
    stop("Target variable must have at least 2 groups after NA removal")
  }

  # === 3. Gene Filter (min_cells, min_pct) ===
  message("V5 Preprocessing: 2. Filtering genes (min_cells, min_pct)...")
  # Use Matrix::rowSums for sparse support
  n_cells_expr <- Matrix::rowSums(expr_mat > 0)
  pct_cells_expr <- n_cells_expr / ncol(expr_mat)
  keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  expr_mat <- expr_mat[keep_genes, , drop = FALSE]

  if (nrow(expr_mat) == 0) stop("No genes pass filtering criteria")

  # === [OPTIMIZATION] Convert to dense NOW ===
  # Only after filtering genes and cells, we convert to dense matrix.
  message(sprintf("V5 Preprocessing: Converting to dense matrix (%d genes x %d cells)...", nrow(expr_mat), ncol(expr_mat)))
  expr_mat <- as.matrix(expr_mat)

  # === 4. [V5] test_n 필터링 (단 한 번 실행) ===
  if (!is.null(test_n) && nrow(expr_mat) > test_n) {
    message(sprintf("V5 Preprocessing: 3. Pre-filtering to top %d genes (limma)...", test_n))
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for 'test_n' pre-filtering")
    }

    if (is.null(control_vars)) {
      design_test <- model.matrix(~target_binary_var, data = meta.data)
    } else {
      formula_test <- as.formula(paste("~ target_binary_var +", paste(control_vars, collapse = "+")))
      design_test <- model.matrix(formula_test, data = meta.data)
    }

    fit_test <- limma::lmFit(expr_mat, design_test)
    fit_test <- limma::eBayes(fit_test)
    coef_indices <- grep("target_binary_var", colnames(design_test))

    top_table_test <- limma::topTable(fit_test, coef = coef_indices, number = Inf, sort.by = "P")

    top_gene_names <- rownames(top_table_test)[1:min(test_n, nrow(top_table_test))]
    expr_mat <- expr_mat[top_gene_names, ]
    message(sprintf("... reduced to %d genes.", nrow(expr_mat)))
  }

  # === 5. Preprocessing (log1p, scale) ===
  message("V5 Preprocessing: 4. Applying log1p and scaling...")

  if (preprocess) {
    if (max(expr_mat) > 100) {
      expr_mat <- log1p(expr_mat)
    }
  }

  expr_mat_scaled <- NULL
  if (length(methods_requiring_scale) > 0) {
    gene_means <- rowMeans(expr_mat)
    gene_sds <- apply(expr_mat, 1, sd)
    gene_sds[gene_sds == 0] <- 1
    expr_mat_scaled <- (expr_mat - gene_means) / gene_sds
  }

  # === 6. [V5] Confounder pre-correction (단 한 번 실행) ===
  message("V5 Preprocessing: 5. Applying confounder correction (if needed)...")

  expr_mat_corrected <- NULL
  if (!is.null(control_vars) && length(methods_requiring_correction) > 0) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for removeBatchEffect")
    }

    covariates_df <- meta.data[, control_vars, drop = FALSE]
    covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
    expr_mat_corrected <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
  }

  # --- 반환 객체 ---
  list(
    expr_mat = expr_mat,
    expr_mat_scaled = expr_mat_scaled,
    expr_mat_corrected = expr_mat_corrected,
    meta.data = meta.data,
    target_binary = target_binary,
    n_groups = length(unique(target_binary)),
    control_vars = control_vars
  )
}

#' Find Gene signature v5.4 (bug fixes for ranger/glmnet/NMF paths)
#'
#' @description Wrapper around [find_gene_signature_v5_impl] that enables the
#' `method_impl = "v5.4"` execution path, activating the latest stability
#' fixes for ranger/glmnet/NMF.
#'
#' @inheritParams find_gene_signature_v5.3
#' @return Same as [find_gene_signature_v5.3]
#' @seealso [find_gene_signature_v5.3]
#' @export
find_gene_signature_v5.4 <- function(...) {
  find_gene_signature_v5_impl(..., method_impl = "v5.4")
}


#' @export
find_gene_signature_v5.3 <- function(data,
                                     meta.data = NULL,
                                     target_var,
                                     target_group = NULL,
                                     control_vars = NULL,
                                     method = c(
                                       "random_forest", "random_forest_ranger",
                                       "lasso", "ridge", "elastic_net",
                                       "pca_loadings", "nmf_loadings",
                                       "gam", "limma", "wilcoxon",
                                       "xgboost"
                                     ),
                                     n_features = 50,
                                     test_n = NULL,
                                     preprocess = TRUE,
                                     min_cells = 10,
                                     min_pct = 0.01,
                                     return_model = FALSE,
                                     fgs_seed = 42,
                                     lambda_selection = "lambda.1se",
                                     enet.alpha = 0.5,
                                     pca.n_pcs = 1,
                                     gam.min_unique = 15,
                                     gam.k = NULL,
                                     gam.k_dynamic_factor = 5,
                                     ...) {
  find_gene_signature_v5_impl(
    data = data,
    meta.data = meta.data,
    target_var = target_var,
    target_group = target_group,
    control_vars = control_vars,
    method = method,
    n_features = n_features,
    test_n = test_n,
    preprocess = preprocess,
    min_cells = min_cells,
    min_pct = min_pct,
    return_model = return_model,
    fgs_seed = fgs_seed,
    lambda_selection = lambda_selection,
    enet.alpha = enet.alpha,
    pca.n_pcs = pca.n_pcs,
    gam.min_unique = gam.min_unique,
    gam.k = gam.k,
    gam.k_dynamic_factor = gam.k_dynamic_factor,
    method_impl = "v5.3",
    ...
  )
}

find_gene_signature_v5_impl <- function(data,
                                        meta.data = NULL,
                                        target_var,
                                        target_group = NULL,
                                        control_vars = NULL,
                                        method = c(
                                          "random_forest", "random_forest_ranger",
                                          "lasso", "ridge", "elastic_net",
                                          "pca_loadings", "nmf_loadings",
                                          "gam", "limma", "wilcoxon",
                                          "xgboost"
                                        ),
                                        n_features = 50,
                                        test_n = NULL,
                                        preprocess = TRUE,
                                        min_cells = 10,
                                        min_pct = 0.01,
                                        return_model = FALSE,
                                        fgs_seed = 42,
                                        # --- 신규/수정된 인자 ---
                                        lambda_selection = "lambda.1se",
                                        enet.alpha = 0.5, # (Elastic Net용)
                                        pca.n_pcs = 1, # (PCA용)
                                        gam.min_unique = 15, # (GAM용)
                                        gam.k = NULL,
                                        gam.k_dynamic_factor = 5,
                                        method_impl = c("v5.3", "v5.4"),
                                        ...) {
  method_impl <- match.arg(method_impl)

  # [Req 5] 메서드 순서 및 이름 변경
  all_methods <- c(
    # 1. Tree-based
    "random_forest", "random_forest_ranger", "xgboost",
    # 2. Regularization
    "lasso", "ridge", "elastic_net",
    # 3. Loadings / Dimensionality Reduction
    "pca_loadings", "nmf_loadings", # nmf -> nmf_loadings
    # 4. Statistical Modelling
    "gam", "limma", "wilcoxon"
  )

  if (is.null(method)) {
    method <- all_methods
  }

  use_dynamic_k <- is.null(gam.k)

  # === 1. [V5] 전처리 (단 1회 실행) ===

  # [Req 4] 신규 모델 스케일링 요구사항 업데이트
  methods_requiring_scale <- c(
    "lasso", "ridge", "elastic_net",
    "gam", "pca_loadings", "xgboost"
  )
  methods_requiring_correction <- c("wilcoxon", "pca_loadings")

  preprocessed_data <- fgs_preprocess_data_v5(
    data = data, meta.data = meta.data, target_var = target_var,
    target_group = target_group, control_vars = control_vars,
    test_n = test_n, preprocess = preprocess, min_cells = min_cells,
    min_pct = min_pct,
    methods_requiring_scale = base::intersect(method, methods_requiring_scale),
    methods_requiring_correction = base::intersect(method, methods_requiring_correction)
  )

  expr_mat_base <- preprocessed_data$expr_mat
  meta.data_clean <- preprocessed_data$meta.data
  target_binary <- preprocessed_data$target_binary
  n_groups <- preprocessed_data$n_groups

  set.seed(fgs_seed)

  covariate_mat_model <- NULL
  if (!is.null(control_vars)) {
    covariates_df_model <- meta.data_clean[, control_vars, drop = FALSE]
    covariate_mat_model <- model.matrix(~ . - 1, data = covariates_df_model)
  }

  results_list <- list()

  run_glmnet_signature <- function(X, alpha_value, apply_v54 = FALSE, ...) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("glmnet package required. Install with: install.packages('glmnet')")
    }

    if (is.null(control_vars)) {
      X_model <- X
      penalty_vec <- rep(1, ncol(X_model))
    } else {
      X_model <- cbind(X, covariate_mat_model)
      penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
    }

    X_model <- as.matrix(X_model)
    family <- if (n_groups == 2) "binomial" else "multinomial"
    y_glmnet <- if (n_groups == 2 && apply_v54) as.numeric(y) - 1 else y

    cv_fit <- glmnet::cv.glmnet(
      X_model, y_glmnet,
      family = family,
      alpha = alpha_value,
      penalty.factor = penalty_vec,
      ...
    )

    if (n_groups == 2) {
      coef_obj <- coef(cv_fit, s = lambda_selection)
      idx <- seq_len(ncol(X)) + 1
      weights_all <- as.numeric(coef_obj[idx])
      names(weights_all) <- rownames(coef_obj)[idx]
    } else {
      coef_list <- coef(cv_fit, s = lambda_selection)
      idx <- seq_len(ncol(X)) + 1
      weights_mat <- vapply(coef_list, function(mat) as.numeric(mat[idx]), numeric(length(idx)))
      weights_all <- rowMeans(weights_mat)
      names(weights_all) <- rownames(coef_list[[1]])[idx]
    }

    weights_all <- weights_all[!is.na(weights_all)]
    if (apply_v54) {
      weights_all <- weights_all[is.finite(weights_all)]
    }

    if (length(weights_all) == 0) {
      stop("glmnet returned no gene coefficients.")
    }

    top_genes <- names(sort(abs(weights_all), decreasing = TRUE)[1:min(n_features, length(weights_all))])
    weights <- weights_all[top_genes]

    scores <- as.numeric(X[, top_genes, drop = FALSE] %*% weights)
    names(scores) <- rownames(X)

    pred_probs <- predict(cv_fit, newx = X_model, s = lambda_selection, type = "response")
    if (n_groups == 2) {
      pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels = levels(y))
      if (requireNamespace("pROC", quietly = TRUE)) {
        roc_obj <- pROC::roc(y, scores, quiet = TRUE)
        auc <- as.numeric(pROC::auc(roc_obj))
      } else {
        auc <- NA
      }
      acc <- mean(pred == y)
      perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
    } else {
      pred_class_indices <- apply(pred_probs[, , 1], 1, which.max)
      pred <- levels(y)[pred_class_indices]
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


  # === 2. [V5] 메서드 순회 (for 루프) ===

  # Progress tracking: load previous timing data if available
  timing_file <- file.path(getwd(), ".fgs_timing_cache.rds")
  timing_cache <- if (file.exists(timing_file)) {
    tryCatch(readRDS(timing_file), error = function(e) list())
  } else {
    list()
  }

  total_methods <- length(method)
  completed_methods <- 0
  total_start_time <- Sys.time()

  # Calculate total estimated time if we have cache for all methods
  total_estimated_sec <- 0
  methods_with_estimate <- 0
  for (m_check in method) {
    if (m_check %in% all_methods) {
      method_key_check <- paste0("fgs_v5.4_", m_check)
      est_sec <- if (!is.null(timing_cache[[method_key_check]])) {
        mean(timing_cache[[method_key_check]], na.rm = TRUE)
      } else {
        NA_real_
      }
      if (!is.na(est_sec) && est_sec > 0) {
        total_estimated_sec <- total_estimated_sec + est_sec
        methods_with_estimate <- methods_with_estimate + 1
      }
    }
  }
  has_total_estimate <- (methods_with_estimate == total_methods && total_estimated_sec > 0)

  if (has_total_estimate) {
    message(sprintf(
      "\n=== FGS v5.4 시작: 총 %d개 메서드, 전체 예상 시간: %.1f분 ===\n",
      total_methods, total_estimated_sec / 60
    ))
  } else {
    message(sprintf("\n=== FGS v5.4 시작: 총 %d개 메서드 ===\n", total_methods))
  }

  for (m_idx in seq_along(method)) {
    m <- method[m_idx]

    if (!m %in% all_methods) {
      warning(sprintf("Invalid method '%s'. Skipping.", m))
      next
    }

    # Estimate time based on previous runs
    method_key <- paste0("fgs_v5.4_", m)
    estimated_sec <- if (!is.null(timing_cache[[method_key]])) {
      mean(timing_cache[[method_key]], na.rm = TRUE)
    } else {
      NA_real_
    }

    if (!is.na(estimated_sec) && estimated_sec > 0) {
      estimated_str <- sprintf("%.1f분", estimated_sec / 60)
      message(sprintf(
        "--- Running Method: %s [%d/%d] (예상: %s) ---",
        m, m_idx, total_methods, estimated_str
      ))
    } else {
      message(sprintf(
        "--- Running Method: %s [%d/%d] ---",
        m, m_idx, total_methods
      ))
    }

    method_start_time <- Sys.time()

    tryCatch(
      {
        # === 3. 데이터 선택 (메서드별) ===

        if (m %in% c("limma", "wilcoxon", "nmf_loadings", "random_forest", "random_forest_ranger")) {
          expr_mat_method <- expr_mat_base
        } else if (m %in% c("lasso", "ridge", "elastic_net", "gam", "pca_loadings", "xgboost")) {
          # 'expr_mat_scaled'가 NULL이면 (필요한 메서드가 없었으면) 원본 사용
          expr_mat_method <- if (is.null(preprocessed_data$expr_mat_scaled)) expr_mat_base else preprocessed_data$expr_mat_scaled
        }

        if (m %in% c("wilcoxon", "pca_loadings")) {
          if (!is.null(preprocessed_data$expr_mat_corrected)) {
            expr_mat_method <- preprocessed_data$expr_mat_corrected
          }
        }

        # (pca_loadings는 스케일링된 보정 데이터가 이상적이나, v5.2는 보정된 비-스케일링 데이터를 우선함)

        X <- t(expr_mat_method)
        y <- target_binary


        # === 4. Method-specific (v5.2 수정) ===

        result <- switch(m,

          # --- 1. Tree-based ---
          random_forest = { # [Req 3] 이름 변경
            if (!requireNamespace("randomForest", quietly = TRUE)) {
              stop("randomForest package required.")
            }

            # RF는 스케일링되지 않은 원본 X 사용
            X_rf <- t(expr_mat_base)
            if (!is.null(control_vars)) {
              X_rf <- cbind(X_rf, covariate_mat_model)
            }

            rf_model <- randomForest::randomForest(x = X_rf, y = y, ntree = 500, importance = TRUE, ...)

            importance_scores <- randomForest::importance(rf_model)
            if (n_groups == 2) {
              weights_magnitude_all <- importance_scores[, "MeanDecreaseGini"]
            } else {
              weights_magnitude_all <- rowMeans(importance_scores[, grep("MeanDecreaseGini", colnames(importance_scores))])
            }

            gene_names_in_model <- colnames(X_rf)[!colnames(X_rf) %in% colnames(covariate_mat_model)]
            weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]

            top_genes <- names(sort(weights_magnitude_genes, decreasing = TRUE)[1:min(n_features, length(weights_magnitude_genes))])
            weights_magnitude <- weights_magnitude_genes[top_genes]

            # 방향성 보정 및 performance 계산
            if (n_groups == 2) {
              g1_cells <- y == levels(y)[1]
              g2_cells <- y == levels(y)[2]
              mean_g1 <- colMeans(X_rf[g1_cells, top_genes, drop = FALSE])
              mean_g2 <- colMeans(X_rf[g2_cells, top_genes, drop = FALSE])
              effect_size <- mean_g2 - mean_g1
              weights <- weights_magnitude * sign(effect_size)
            } else {
              warning("random_forest: n_groups > 2. Score represents magnitude (importance), not direction.")
              weights <- weights_magnitude
            }

            scores <- as.numeric(X_rf[, top_genes] %*% weights)
            names(scores) <- rownames(X_rf)

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
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = if (return_model) rf_model else NULL
            )
          },
          random_forest_ranger = {
            if (!requireNamespace("ranger", quietly = TRUE)) {
              stop("ranger package required.")
            }

            X_ranger <- t(expr_mat_base)
            if (!is.null(control_vars)) {
              X_ranger <- cbind(X_ranger, covariate_mat_model)
            }

            # Store original names and create mapping for restoration
            orig_colnames <- colnames(X_ranger)
            # Mapping: make.names(orig) -> orig
            name_mapping <- setNames(orig_colnames, make.names(orig_colnames))

            ranger_data <- data.frame(y = y, X_ranger)
            importance_mode <- if (method_impl == "v5.4") "permutation" else "impurity"

            rf_model <- tryCatch(
              ranger::ranger(
                y ~ .,
                data = ranger_data,
                num.trees = 500,
                importance = importance_mode,
                ...
              ),
              error = function(e) {
                if (method_impl == "v5.4" && importance_mode == "permutation") {
                  warning("Permutation importance failed in ranger; falling back to impurity. ", e$message)
                  ranger::ranger(
                    y ~ .,
                    data = ranger_data,
                    num.trees = 500,
                    importance = "impurity",
                    ...
                  )
                } else {
                  stop(e)
                }
              }
            )

            weights_magnitude_all <- ranger::importance(rf_model)

            # Map names back to original before filtering/sorting
            # Only map if names don't match original (ranger modifies names)
            if (!all(names(weights_magnitude_all) %in% orig_colnames)) {
              mapped_names <- name_mapping[names(weights_magnitude_all)]
              # Keep only valid mappings
              valid_idx <- !is.na(mapped_names)
              weights_magnitude_all <- weights_magnitude_all[valid_idx]
              names(weights_magnitude_all) <- mapped_names[valid_idx]
            }

            if (!is.null(control_vars)) {
              weights_magnitude_all <- weights_magnitude_all[!names(weights_magnitude_all) %in% colnames(covariate_mat_model)]
            }
            weights_magnitude_all <- weights_magnitude_all[is.finite(weights_magnitude_all)]

            top_genes <- names(sort(weights_magnitude_all, decreasing = TRUE)[1:min(n_features, length(weights_magnitude_all))])
            weights_magnitude <- weights_magnitude_all[top_genes]

            if (n_groups == 2) {
              g1_cells <- y == levels(y)[1]
              g2_cells <- y == levels(y)[2]
              mean_g1 <- colMeans(X_ranger[g1_cells, top_genes, drop = FALSE])
              mean_g2 <- colMeans(X_ranger[g2_cells, top_genes, drop = FALSE])
              effect_size <- mean_g2 - mean_g1
              weights <- weights_magnitude * sign(effect_size)
            } else {
              warning("random_forest_ranger: n_groups > 2. Score represents magnitude (importance), not direction.")
              weights <- weights_magnitude
            }

            scores <- as.numeric(X_ranger[, top_genes, drop = FALSE] %*% weights)
            names(scores) <- rownames(X_ranger)

            pred <- rf_model$predictions
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
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = if (return_model) rf_model else NULL
            )
          },
          xgboost = { # [Req 4] 신규 (Beta)
            if (!requireNamespace("xgboost", quietly = TRUE)) {
              stop("xgboost package required.")
            }

            # X는 스케일링된 데이터 (expr_mat_scaled)
            X_xgb <- X
            if (!is.null(control_vars)) {
              X_xgb <- cbind(X_xgb, covariate_mat_model)
            }

            # y를 0/1 숫자로 변환
            y_numeric <- as.numeric(y) - 1

            dtrain <- xgboost::xgb.DMatrix(data = X_xgb, label = y_numeric)

            params <- list(
              objective = if (n_groups == 2) "binary:logistic" else "multi:softmax",
              eval_metric = if (n_groups == 2) "logloss" else "mlogloss",
              nthread = 4,
              eta = 0.1
            )
            if (n_groups > 2) {
              params$num_class <- n_groups
            }

            xgb_model <- xgboost::xgb.train(params, dtrain, nrounds = 100, ...)

            imp_matrix <- xgboost::xgb.importance(model = xgb_model)
            weights_magnitude_all <- imp_matrix$Gain
            names(weights_magnitude_all) <- imp_matrix$Feature

            gene_names_in_model <- names(weights_magnitude_all)[!names(weights_magnitude_all) %in% colnames(covariate_mat_model)]
            weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]

            top_genes <- names(sort(weights_magnitude_genes, decreasing = TRUE)[1:min(n_features, length(weights_magnitude_genes))])
            weights_magnitude <- weights_magnitude_genes[top_genes]

            # 방향성 보정 및 performance 계산
            if (n_groups == 2) {
              g1_cells <- y == levels(y)[1]
              g2_cells <- y == levels(y)[2]
              mean_g1 <- colMeans(X_xgb[g1_cells, top_genes, drop = FALSE])
              mean_g2 <- colMeans(X_xgb[g2_cells, top_genes, drop = FALSE])
              effect_size <- mean_g2 - mean_g1
              weights <- weights_magnitude * sign(effect_size)
            } else {
              warning("xgboost: n_groups > 2. Score represents magnitude (importance), not direction.")
              weights <- weights_magnitude
            }

            scores <- as.numeric(X_xgb[, top_genes] %*% weights)
            names(scores) <- rownames(X_xgb)

            pred_probs <- predict(xgb_model, newdata = X_xgb)
            if (n_groups == 2) {
              pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels = levels(y))
              if (requireNamespace("pROC", quietly = TRUE)) {
                roc_obj <- pROC::roc(y, scores, quiet = TRUE)
                auc <- as.numeric(pROC::auc(roc_obj))
              } else {
                auc <- NA
              }
              acc <- mean(pred == y)
              perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
            } else {
              pred_class_indices <- apply(pred_probs, 1, which.max)
              pred <- levels(y)[pred_class_indices]
              acc <- mean(pred == y)
              perf <- list(accuracy = acc, confusion = table(pred, y))
            }

            list(
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = if (return_model) xgb_model else NULL
            )
          },

          # --- 2. Regularization ---
          lasso = run_glmnet_signature(X, alpha_value = 1, apply_v54 = FALSE, ...),
          ridge = run_glmnet_signature(
            X,
            alpha_value = 0,
            apply_v54 = (method_impl == "v5.4"),
            ...
          ),
          elastic_net = run_glmnet_signature(
            X,
            alpha_value = enet.alpha,
            apply_v54 = (method_impl == "v5.4"),
            ...
          ),

          # --- 3. Loadings ---
          pca_loadings = {
            # [Req 2] n_pcs 로직 적용

            # expr_mat_method는 보정O, 스케일링X (또는 스케일링O, 보정X)
            # v5.2는 prcomp 전에 다시 스케일링을 보장
            pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

            n_pcs_to_test <- min(50, ncol(pca_res$x))
            pc_cors <- numeric(n_pcs_to_test)

            for (k in 1:n_pcs_to_test) {
              if (n_groups == 2) {
                pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
              } else {
                pc_cors[k] <- summary(aov(pca_res$x[, k] ~ target_binary))[[1]][1, "F value"]
              }
            }

            # 1. 가장 상관관계 높은 (Best) PC 1개
            best_pc_idx <- which.max(pc_cors)

            # 2. 상위 N개 (Top N) PC
            top_n_pcs_indices <- order(pc_cors, decreasing = TRUE)[1:min(pca.n_pcs, n_pcs_to_test)]

            # 3. 중요도 계산: Top N개 PC의 로딩 절대값 합
            top_pc_loadings <- pca_res$rotation[, top_n_pcs_indices, drop = FALSE]
            weights_magnitude_all <- rowSums(abs(top_pc_loadings))

            # 4. 방향성 계산: Best PC 1개의 부호있는 로딩
            weights_directional_all <- pca_res$rotation[, best_pc_idx]

            # 5. 유전자 선별: Top N개 기준
            top_genes <- names(sort(weights_magnitude_all, decreasing = TRUE)[1:min(n_features, length(weights_magnitude_all))])

            # 6. 최종 가중치: Best PC 1개의 방향성 적용
            weights <- weights_directional_all[top_genes]

            scores <- as.numeric(X[, top_genes] %*% weights)
            names(scores) <- rownames(X)

            perf <- list(
              Best_PC = best_pc_idx,
              Top_N_PCs = top_n_pcs_indices,
              Top_N_Correlations = pc_cors[top_n_pcs_indices],
              Best_PC_VarExplained = summary(pca_res)$importance[2, best_pc_idx]
            )

            list(
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = if (return_model) pca_res else NULL
            )
          },
          nmf_loadings = {
            if (!requireNamespace("NMF", quietly = TRUE)) {
              stop("NMF package required.")
            }

            expr_mat_nmf <- expr_mat_method

            if (!is.null(control_vars)) {
              warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s...", m))
              covariate_mat <- model.matrix(~ . - 1, data = meta.data_clean[, control_vars, drop = FALSE])
              expr_mat_nmf <- limma::removeBatchEffect(expr_mat_nmf, covariates = covariate_mat)
            }

            min_expr <- suppressWarnings(min(expr_mat_nmf, na.rm = TRUE))
            eps <- 1e-3
            shift <- if (method_impl == "v5.4") {
              if (is.finite(min_expr) && min_expr <= 0) abs(min_expr) + eps else eps
            } else {
              -min_expr + 0.01
            }
            expr_mat_pos <- expr_mat_nmf + shift

            rank_cap <- min(n_groups + 2, 10)
            if (method_impl == "v5.4") {
              max_rank_allowed <- min(nrow(expr_mat_pos) - 1, ncol(expr_mat_pos) - 1)
              if (!is.finite(max_rank_allowed) || max_rank_allowed < 2) {
                stop("nmf_loadings: insufficient dimensions to fit requested rank.")
              }
              rank <- max(2, min(rank_cap, max_rank_allowed))
            } else {
              rank <- rank_cap
            }
            dots <- list(...)

            # [Req 6 FIX] '...' (dots)에서 NMF 유효 인자만 필터링
            valid_nmf_args <- c("nrun", ".options", ".pbackend")
            filtered_dots <- dots[names(dots) %in% valid_nmf_args]

            nmf_args <- c(
              list(
                x = expr_mat_pos,
                rank = rank,
                seed = fgs_seed,
                method = "brunet"
              ), # 기본 알고리즘 명시
              filtered_dots
            )

            # Ensure NMF is loaded for registry access
            if (!"package:NMF" %in% search()) {
              suppressPackageStartupMessages(library(NMF))
            }
            # Use string "nmf" to allow S4 dispatch to work correctly with registry
            nmf_res <- do.call("nmf", nmf_args)

            W <- NMF::basis(nmf_res)
            H <- NMF::coef(nmf_res)

            component_cors <- numeric(rank)
            for (k in 1:rank) {
              if (n_groups == 2) {
                component_cors[k] <- abs(cor(H[k, ], as.numeric(target_binary)))
              } else {
                component_cors[k] <- summary(aov(H[k, ] ~ target_binary))[[1]][1, "F value"]
              }
            }

            best_component <- which.max(component_cors)
            weights_magnitude <- W[, best_component]
            names(weights_magnitude) <- rownames(expr_mat_pos)

            top_genes <- names(sort(weights_magnitude, decreasing = TRUE)[1:min(n_features, length(weights_magnitude))])
            weights_magnitude <- weights_magnitude[top_genes]

            if (n_groups == 2) {
              g1_cells <- y == levels(y)[1]
              g2_cells <- y == levels(y)[2]
              mean_g1 <- colMeans(X[g1_cells, top_genes, drop = FALSE])
              mean_g2 <- colMeans(X[g2_cells, top_genes, drop = FALSE])
              effect_size <- mean_g2 - mean_g1
              weights <- weights_magnitude * sign(effect_size)
            } else {
              warning("nmf_loadings: n_groups > 2. Score represents magnitude (importance), not direction.")
              weights <- weights_magnitude
            }

            scores <- as.numeric(X[, top_genes] %*% weights)
            names(scores) <- rownames(X)

            perf <- list(component = best_component, correlation = component_cors[best_component])

            list(
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = if (return_model) nmf_res else NULL
            )
          },

          # --- 4. Statistical Modelling ---
          gam = {
            if (!requireNamespace("mgcv", quietly = TRUE)) {
              stop("mgcv package required. Install with: install.packages('mgcv')")
            }

            deviance_explained <- numeric(nrow(expr_mat_method))
            names(deviance_explained) <- rownames(expr_mat_method)

            gam_data <- data.frame(y_var_numeric = as.numeric(target_binary) - 1)
            if (!is.null(control_vars)) {
              gam_data <- cbind(gam_data, covariate_mat_model)
              formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse = " + "))
            } else {
              formula_base <- ""
            }

            convergence_warnings <- 0
            genes_skipped <- 0
            genes_with_dynamic_k <- 0

            for (i in 1:nrow(expr_mat_method)) {
              gam_data$gene_expr <- expr_mat_method[i, ]

              n_unique_vals <- length(base::unique(gam_data$gene_expr))
              if (n_unique_vals < gam.min_unique) {
                deviance_explained[i] <- 0
                genes_skipped <- genes_skipped + 1
                next
              }

              k_value <- if (use_dynamic_k) {
                k_dyn <- floor(n_unique_vals / gam.k_dynamic_factor)
                k_dyn <- max(3, min(10, k_dyn))
                if (k_dyn != 10) genes_with_dynamic_k <- genes_with_dynamic_k + 1
                k_dyn
              } else {
                gam.k
              }

              formula_str <- paste("y_var_numeric ~ s(gene_expr, k=", k_value, ", bs='cr')", formula_base)

              fit_result <- tryCatch(
                {
                  if (n_groups == 2) {
                    mgcv::bam(as.formula(formula_str), data = gam_data, family = "binomial", ...)
                  } else {
                    mgcv::bam(as.formula(formula_str), data = gam_data, ...)
                  }
                },
                warning = function(w) {
                  if (grepl("did not converge", w$message)) {
                    convergence_warnings <<- convergence_warnings + 1
                  }
                  invokeRestart("muffleWarning")
                },
                error = function(e) {
                  warning(sprintf("GAM failed for gene %s: %s", rownames(expr_mat_method)[i], e$message))
                  return(NULL)
                }
              )

              if (is.null(fit_result)) {
                deviance_explained[i] <- 0
              } else {
                deviance_explained[i] <- summary(fit_result)$dev.expl
              }
            }

            if (genes_skipped > 0) {
              warning(sprintf("GAM: Skipped %d genes (unique values < gam.min_unique=%d).", genes_skipped, gam.min_unique))
            }
            if (convergence_warnings > 0) {
              warning(sprintf("GAM: %d genes failed to converge.", convergence_warnings))
            }
            if (use_dynamic_k && genes_with_dynamic_k > 0) {
              message(sprintf("GAM: Applied dynamic k for %d genes.", genes_with_dynamic_k))
            }

            weights_magnitude <- deviance_explained
            top_genes <- names(sort(weights_magnitude, decreasing = TRUE)[1:min(n_features, sum(weights_magnitude > 0, na.rm = TRUE))])

            if (length(top_genes) == 0) {
              warning("GAM method found no genes with deviance > 0.")
              return(list(genes = character(0), weights = numeric(0), scores = numeric(0), performance = list()))
            }

            weights_magnitude <- weights_magnitude[top_genes]

            if (n_groups == 2) {
              g1_cells <- y == levels(y)[1]
              g2_cells <- y == levels(y)[2]
              mean_g1 <- colMeans(X[g1_cells, top_genes, drop = FALSE])
              mean_g2 <- colMeans(X[g2_cells, top_genes, drop = FALSE])
              effect_size <- mean_g2 - mean_g1
              weights <- weights_magnitude * sign(effect_size)
            } else {
              warning("gam: n_groups > 2. Score represents magnitude (importance), not direction.")
              weights <- weights_magnitude
            }

            scores <- as.numeric(X[, top_genes, drop = FALSE] %*% weights)
            names(scores) <- rownames(X)

            perf <- list(deviance_explained = weights_magnitude)

            list(
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = NULL
            )
          },
          limma = {
            if (!requireNamespace("limma", quietly = TRUE)) {
              stop("limma package required. Install with: BiocManager::install('limma')")
            }

            # Sanitize target_binary levels for limma (must be valid R names)
            target_binary_limma <- target_binary
            orig_levels <- levels(target_binary_limma)
            safe_levels <- make.names(orig_levels)
            if (!all(orig_levels == safe_levels)) {
              levels(target_binary_limma) <- safe_levels
              meta.data_clean$target_binary_limma <- target_binary_limma
            } else {
              meta.data_clean$target_binary_limma <- target_binary_limma
            }

            # expr_mat_method는 'expr_mat_base' (보정되지 않은 원본)
            if (is.null(control_vars)) {
              design <- model.matrix(~ 0 + target_binary_limma, data = meta.data_clean)
              colnames(design)[1:n_groups] <- levels(target_binary_limma)
            } else {
              control_formula <- paste(control_vars, collapse = " + ")
              full_formula <- as.formula(paste("~0 + target_binary_limma +", control_formula))
              design <- model.matrix(full_formula, data = meta.data_clean)
              colnames(design)[1:n_groups] <- levels(target_binary_limma)
            }

            fit <- limma::lmFit(expr_mat_method, design)

            if (n_groups == 2) {
              contrast_str <- paste(levels(target_binary_limma)[2], levels(target_binary_limma)[1], sep = "-")
            } else {
              contrast_pairs <- combn(levels(target_binary_limma), 2, function(x) paste(x[2], x[1], sep = "-"))
              contrast_str <- contrast_pairs
              warning("limma: n_groups > 2. Using all pairwise contrasts. Weights based on average t-statistic.")
            }

            contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)
            fit2 <- limma::contrasts.fit(fit, contrast_mat)
            fit2 <- limma::eBayes(fit2)

            if (n_groups == 2) {
              top_table <- limma::topTable(fit2, number = Inf, sort.by = "B")
              weights_all <- top_table$t
            } else {
              top_table <- limma::topTable(fit2, number = Inf, sort.by = "F")
              weights_all <- rowMeans(abs(fit2$t[, contrast_str, drop = FALSE]))
            }

            names(weights_all) <- rownames(top_table)

            top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]

            if (n_groups == 2) {
              weights <- top_table[top_genes, "t"]
            } else {
              weights <- weights_all[top_genes]
            }

            scores <- as.numeric(X[, top_genes] %*% weights)
            names(scores) <- rownames(X)

            perf <- list(top_table = top_table[top_genes, ])

            list(
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = if (return_model) fit2 else NULL
            )
          },
          wilcoxon = {
            # expr_mat_method는 'expr_mat_corrected' (보정된 버전)
            pvals <- numeric(nrow(expr_mat_method))
            effect_sizes <- numeric(nrow(expr_mat_method))
            names(pvals) <- rownames(expr_mat_method)
            names(effect_sizes) <- rownames(expr_mat_method)

            for (i in 1:nrow(expr_mat_method)) {
              if (n_groups == 2) {
                group1 <- expr_mat_method[i, target_binary == levels(target_binary)[1]]
                group2 <- expr_mat_method[i, target_binary == levels(target_binary)[2]]

                test <- try(wilcox.test(group1, group2), silent = TRUE)

                if (inherits(test, "try-error")) {
                  pvals[i] <- 1.0
                  effect_sizes[i] <- 0
                } else {
                  pvals[i] <- test$p.value
                  effect_sizes[i] <- median(group2) - median(group1)
                }
              } else {
                test <- try(kruskal.test(expr_mat_method[i, ] ~ target_binary), silent = TRUE)

                if (inherits(test, "try-error")) {
                  pvals[i] <- 1.0
                  effect_sizes[i] <- 0
                } else {
                  pvals[i] <- test$p.value
                  effect_sizes[i] <- var(tapply(expr_mat_method[i, ], target_binary, median))
                }
              }
            }

            padj <- p.adjust(pvals, method = "BH")

            ranks <- rank(pvals) + rank(-abs(effect_sizes))
            top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])

            weights <- effect_sizes[top_genes]

            scores <- as.numeric(X[, top_genes] %*% weights)
            names(scores) <- rownames(X)

            perf <- list(
              pvalues = pvals[top_genes], padj = padj[top_genes],
              effect_sizes = effect_sizes[top_genes]
            )

            list(
              genes = top_genes, weights = weights, scores = scores,
              performance = perf, model = NULL
            )
          }
        ) # --- switch 끝 ---

        # === 5. 결과 저장 ===
        result$method <- m
        result$target_var <- target_var
        result$n_groups <- n_groups
        result$n_cells <- ncol(expr_mat_base)

        class(result) <- c("gene_signature", "list")
        results_list[[m]] <- result

        # Record execution time
        method_end_time <- Sys.time()
        elapsed_sec <- as.numeric(difftime(method_end_time, method_start_time, units = "secs"))
        completed_methods <- completed_methods + 1

        # Update timing cache
        method_key <- paste0("fgs_v5.4_", m)
        if (is.null(timing_cache[[method_key]])) {
          timing_cache[[method_key]] <- numeric()
        }
        timing_cache[[method_key]] <- c(timing_cache[[method_key]], elapsed_sec)
        # Keep only last 5 runs for averaging
        if (length(timing_cache[[method_key]]) > 5) {
          timing_cache[[method_key]] <- tail(timing_cache[[method_key]], 5)
        }

        # Save timing cache
        tryCatch(
          {
            saveRDS(timing_cache, timing_file)
          },
          error = function(e) {
            # Ignore save errors
          }
        )

        # Improved progress estimation: use actual elapsed time from first completed method
        # After first method completes, we can estimate remaining time more accurately
        elapsed_total <- as.numeric(difftime(method_end_time, total_start_time, units = "secs"))
        remaining_methods <- total_methods - completed_methods

        if (completed_methods == 1 && remaining_methods > 0) {
          # First method completed: estimate based on this method's time
          # So remaining methods will take approximately: elapsed_sec * remaining_methods
          estimated_remaining <- elapsed_sec * remaining_methods
          message(sprintf(
            "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (첫 메서드 기준)",
            m, elapsed_sec, elapsed_sec / 60, completed_methods, total_methods, estimated_remaining / 60
          ))
        } else if (completed_methods > 1 && remaining_methods > 0) {
          # Multiple methods completed: use average time
          avg_time_per_method <- elapsed_total / completed_methods
          estimated_remaining <- avg_time_per_method * remaining_methods
          message(sprintf(
            "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (평균 %.1f초/메서드)",
            m, elapsed_sec, elapsed_sec / 60, completed_methods, total_methods,
            estimated_remaining / 60, avg_time_per_method
          ))
        } else {
          # Last method or no remaining methods
          message(sprintf(
            "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 완료!",
            m, elapsed_sec, elapsed_sec / 60, completed_methods, total_methods
          ))
        }
      },
      error = function(e) {
        warning(sprintf("Method '%s' failed with error: %s", m, e$message))
        results_list[[m]] <<- list(method = m, error = e$message)

        # Record execution time even for errors
        method_end_time <- Sys.time()
        elapsed_sec <- as.numeric(difftime(method_end_time, method_start_time, units = "secs"))
        completed_methods <<- completed_methods + 1

        # Calculate remaining time estimate even for errors
        elapsed_total <- as.numeric(difftime(method_end_time, total_start_time, units = "secs"))
        remaining_methods <- total_methods - completed_methods

        if (completed_methods == 1 && remaining_methods > 0) {
          estimated_remaining <- elapsed_sec * remaining_methods
          message(sprintf(
            "✗ %s 실패: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (첫 메서드 기준)",
            m, elapsed_sec, elapsed_sec / 60, completed_methods, total_methods, estimated_remaining / 60
          ))
        } else if (completed_methods > 1 && remaining_methods > 0) {
          avg_time_per_method <- elapsed_total / completed_methods
          estimated_remaining <- avg_time_per_method * remaining_methods
          message(sprintf(
            "✗ %s 실패: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (평균 %.1f초/메서드)",
            m, elapsed_sec, elapsed_sec / 60, completed_methods, total_methods,
            estimated_remaining / 60, avg_time_per_method
          ))
        } else {
          message(sprintf(
            "✗ %s 실패: %.1f초 (%.1f분) | 진행: %d/%d",
            m, elapsed_sec, elapsed_sec / 60, completed_methods, total_methods
          ))
        }
      }
    )
  } # --- for 루프 끝 ---

  # Final summary
  total_end_time <- Sys.time()
  total_elapsed <- as.numeric(difftime(total_end_time, total_start_time, units = "secs"))
  message(sprintf(
    "\n=== FGS v5.4 완료: 총 %d개 메서드, %.1f분 소요 ===",
    completed_methods, total_elapsed / 60
  ))

  return(results_list)
}

#' Train Meta-Learner (TML7): Stacked Ensemble Model for Signature Scores
#'
#' @description
#'   Implements a two-level stacking approach for combining multiple Level-1 (L1)
#'   gene signatures into a unified prediction model. The function standardizes
#'   diverse signature formats (character vectors, named numeric weights,
#'   up/down lists, data frames) into a common representation, computes signature
#'   scores on the holdout expression data, and trains one or more Level-2 (L2)
#'   caret models on the resulting score matrix. The best-performing model,
#'   selected via cross-validation, is returned along with the processed training
#'   data and standardized signatures.
#'
#'   **Mechanism:**
#'   1. **Signature Standardization**: Each L1 signature is converted to a named
#'      numeric vector of gene weights using `as_signature()`, which handles
#'      multiple input formats (character vectors → uniform weights, named
#'      numeric → direct weights, up/down lists → +1/-1 weights, data frames →
#'      extracted gene-weight pairs).
#'   2. **Score Computation**: For each signature, a per-cell score is computed
#'      as the weighted average of expression values, normalized by the sum of
#'      absolute weights and optionally z-scored across cells.
#'   3. **L2 Feature Matrix**: The signature scores form the columns of the L2
#'      training matrix (one row per cell, one column per signature).
#'   4. **Model Training**: Multiple caret models (e.g., glm, ranger, xgbTree)
#'      are trained via k-fold cross-validation, with the best model selected
#'      based on the specified metric (AUC/ROC for binary, Accuracy/Kappa for
#'      multi-class).
#'   5. **Parallel Safety**: Parallel execution is opt-in (`allow_parallel=TRUE`)
#'      and guarded to prevent runaway worker creation on shared servers.
#'      Worker count is capped at 8 by default.
#'
#' @param l1_signatures Named list of L1 signatures. Each element can be:
#'   \itemize{
#'     \item A character vector of gene names (uniform weights = 1)
#'     \item A named numeric vector (gene names → weights)
#'     \item A list with `up` and/or `down` components (up-regulated genes get
#'           weight +1, down-regulated get -1)
#'     \item A data frame with gene identifiers and weight/score columns
#'   }
#' @param holdout_data Expression container. Either:
#'   \itemize{
#'     \item A Seurat object (expression extracted via `GetAssayData(layer=layer)`)
#'     \item A matrix or data frame with rows = genes/features, columns = cells/samples
#'   }
#' @param target_var For Seurat input: column name in `meta.data` containing the
#'   outcome to predict. For matrix input: the target vector itself (must match
#'   cell order).
#' @param l2_methods Character vector of caret model identifiers to evaluate
#'   (default: `c("glm", "ranger", "xgbTree")`). Supported values:
#'   `c("glm","ranger","xgbTree","glmnet","svmRadial","mlp",
#'     "mlpKerasDropout","nnet","earth")`. Any unsupported method names are
#'   dropped with a warning; each candidate requires its backing package
#'   (e.g. `glmnet`, `kernlab`, `RSNNS`, `keras`, `earth`).
#' @param k_folds Number of cross-validation folds (default: 5).
#' @param metric Performance metric for model selection. Options:
#'   \itemize{
#'     \item `"AUC"` or `"ROC"`: For binary classification (requires `twoClassSummary`)
#'     \item `"Accuracy"`: Classification accuracy
#'     \item `"Kappa"`: Cohen's kappa
#'   }
#'   For multi-class targets, AUC/ROC is invalid and falls back to Accuracy.
#' @param fgs_seed Random seed for reproducibility (default: 42).
#' @param layer When `holdout_data` is a Seurat object, the assay layer to extract
#'   (default: `"data"`). Options: `"counts"`, `"data"`, `"scale.data"`.
#' @param allow_parallel Logical. If `TRUE` and `future`/`doFuture` packages are
#'   available, enables parallel cross-validation via `future::multisession` plan
#'   (default: `FALSE` for safety).
#' @param parallel_workers Optional integer overriding the number of parallel
#'   workers when `allow_parallel=TRUE`. If `NULL`, defaults to
#'   `min(8, detectCores(logical=FALSE))`, falling back to 4 if detection fails.
#' @param cv_group_var Optional metadata column (Seurat input only) used to
#'   define group-wise cross-validation splits (default: `"emrid"`). When the
#'   column is missing, contains `NA`, or when `holdout_data` is not Seurat,
#'   the function automatically falls back to standard cell-wise CV with a log
#'   message describing the reason.
#'
#' @return A list with components:
#'   \describe{
#'     \item{best_model}{Caret `train` object for the selected L2 model.}
#'     \item{best_model_name}{Character string identifier of the best model.}
#'     \item{best_metric_name}{Name of the metric used for selection.}
#'     \item{model_comparison}{A `resamples` object comparing all trained models
#'           (or `NULL` if only one model succeeded).}
#'     \item{trained_models}{Named list of all successful caret `train` objects
#'           prior to best-model selection.}
#'     \item{l2_train}{Data frame with signature scores (columns) and target
#'           variable (`.target` column).}
#'     \item{l1_signatures}{Standardized signatures (named numeric weight vectors).}
#'     \item{positive_class}{Name of the positive class (tail of `.target`
#'           levels) stored for downstream importance utilities.}
#'   }
#'
#' @details
#'   **Signature Scoring Formula:**
#'   For a signature with gene weights \eqn{w_g} and expression values \eqn{x_{gc}}
#'   for gene \eqn{g} and cell \eqn{c}:
#'   \deqn{s_c = \frac{\sum_g w_g \cdot x_{gc}}{\sum_g |w_g|}}
#'   If `normalize=TRUE`, scores are further z-scored: \eqn{s_c' = (s_c - \mu_s) / \sigma_s}.
#'
#'   **Model Selection:**
#'   The best model is chosen by maximizing the cross-validated metric across all
#'   folds. For binary classification with AUC/ROC, `caret::twoClassSummary` is
#'   used; for multi-class or Accuracy/Kappa, `caret::defaultSummary` is used.
#'
#'   **Group-wise CV Logging:**
#'   When `cv_group_var` is available, folds are created at the group level and
#'   per-fold index/indexOut sizes (cells + groups) are logged to aid debugging.
#'   Missing metadata (or insufficient group counts) triggers an automatic
#'   fallback to standard cell-wise folds.
#'
#'   **L2 Registry:**
#'   Only registered caret methods are allowed. Dependencies are verified via
#'   `requireNamespace()` and missing packages cause the corresponding methods
#'   to be dropped before training.
#'
#' @examples
#' \dontrun{
#' # Example: Combine multiple signature-finding methods
#' sigs <- list(
#'   lasso = find_gene_signature(data, target_var = "g3", method = "lasso"),
#'   rf = find_gene_signature(data, target_var = "g3", method = "tree_based"),
#'   limma = find_gene_signature(data, target_var = "g3", method = "limma")
#' )
#'
#' # Train meta-learner on holdout data
#' meta_model <- TML7(
#'   l1_signatures = sigs,
#'   holdout_data = seurat_obj,
#'   target_var = "g3",
#'   l2_methods = c("glm", "ranger", "xgbTree"),
#'   metric = "AUC",
#'   cv_folds = 5,
#'   cv_method = "cv",
#'   repeats = 1,
#'   cv_group_var = NULL,
#'   fgs_seed = 42,
#'   cv_group_var = "emrid" # Group-wise CV to prevent patient-level leakage
#' )
#'
#' # Use for prediction
#' predictions <- predict(meta_model$best_model, newdata = new_scores)
#' }
#'
#' @seealso
#'   \code{\link[caret]{train}} for L2 model training,
#'   \code{\link{compute_meta_gene_importance}} for deriving gene-level
#'   importance from the trained meta-learner
#'
#' @export

#' Find Gene Signature (FGS)
#'
#' @description
#' Alias for `find_gene_signature_v5.4`. A comprehensive gene signature
#' discovery function supporting multiple methods.
#'
#' @param ... All arguments passed to `find_gene_signature_v5.4`
#' @return See `find_gene_signature_v5.4`
#' @export
FGS <- function(...) {
  find_gene_signature_v5.4(...)
}
