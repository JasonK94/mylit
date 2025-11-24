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
  
  # === 1. Input validation (v4.2의 안정적 로직) ===
  is_seurat <- inherits(data, "Seurat")
  
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required but not installed")
    }
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    default_assay <- Seurat::DefaultAssay(data)
    expr_mat <- tryCatch({
        as.matrix(Seurat::GetAssayData(data, assay = default_assay, layer = "data"))
    }, error = function(e) {
        tryCatch({
            as.matrix(Seurat::GetAssayData(data, assay = default_assay, slot = "data"))
        }, error = function(e) {
            as.matrix(Seurat::GetAssayData(data, assay = default_assay, slot = "counts"))
        })
    })
  } else {
    if (is.null(meta.data)) {
      stop("meta.data must be provided when data is not a Seurat object")
    }
    expr_mat <- as.matrix(data)
  }
  
  # === 2. [V5] NA 및 팩터 레벨 완벽 제거 ===
  message("V5 Preprocessing: 1. Cleaning NA values...")
  
  if (is.null(control_vars)) {
    vars_to_check <- target_var
  } else {
    vars_to_check <- c(target_var, control_vars)
  }
  
  complete_cases_idx <- complete.cases(meta.data[, vars_to_check, drop=FALSE])
  
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
  
  common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
  if (length(common_cells) == 0) {
    stop("No overlapping cells between expression data and metadata.")
  }
  meta.data <- meta.data[common_cells, , drop = FALSE]
  expr_mat <- expr_mat[, common_cells, drop = FALSE]

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
                           ifelse(target_values >= high_cutoff, "High", NA))
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
  n_cells_expr <- rowSums(expr_mat > 0)
  pct_cells_expr <- n_cells_expr / ncol(expr_mat)
  keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  expr_mat <- expr_mat[keep_genes, ]
  
  if (nrow(expr_mat) == 0) stop("No genes pass filtering criteria")
  
  # === 4. [V5] test_n 필터링 (단 한 번 실행) ===
  if (!is.null(test_n) && nrow(expr_mat) > test_n) {
    message(sprintf("V5 Preprocessing: 3. Pre-filtering to top %d genes (limma)...", test_n))
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for 'test_n' pre-filtering")
    }
    
    if (is.null(control_vars)) {
      design_test <- model.matrix(~ target_binary_var, data = meta.data)
    } else {
      formula_test <- as.formula(paste("~ target_binary_var +", paste(control_vars, collapse="+")))
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

