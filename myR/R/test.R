

# NOTE: Package dependencies should be declared in DESCRIPTION, not with library() calls

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
  
  # Input validation
  if (!inherits(sobj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  layer <- match.arg(layer)
  effect <- match.arg(effect)
  regressor.type <- match.arg(regressor.type)
  ordinal.method <- match.arg(ordinal.method)
  
  # Check if regressor exists in metadata
  if (!regressor %in% colnames(sobj@meta.data)) {
    stop(paste("Regressor", regressor, "not found in metadata"))
  }
  
  # Get expression data based on layer
  if (layer == "counts") {
    expr_data <- GetAssayData(sobj, layer = "counts")
  } else if (layer == "data") {
    expr_data <- GetAssayData(sobj, layer = "data")
  } else if (layer == "scale.data") {
    expr_data <- GetAssayData(sobj, layer = "scale.data")
  }
  
  # Select features
  if (features[1] == "all") {
    features <- rownames(expr_data)
  } else {
    # Check if all features exist
    missing_features <- setdiff(features, rownames(expr_data))
    if (length(missing_features) > 0) {
      warning(paste("Features not found:", paste(missing_features, collapse = ", ")))
      features <- intersect(features, rownames(expr_data))
    }
  }
  
  # Filter features by minimum cell expression
  if (min.cells > 0) {
    feature_counts <- Matrix::rowSums(expr_data[features, ] > 0)
    features <- features[feature_counts >= min.cells]
  }
  
  if (length(features) == 0) {
    stop("No features passed filtering criteria")
  }
  
  # Prepare metadata
  meta_data <- sobj@meta.data
  regressor_values <- meta_data[[regressor]]
  
  # Process regressor based on type
  if (regressor.type == "categorical") {
    regressor_values <- as.factor(regressor_values)
    
    # Set reference level if specified
    if (!is.null(reference.level)) {
      if (!reference.level %in% levels(regressor_values)) {
        stop(paste("Reference level", reference.level, "not found in regressor levels"))
      }
      regressor_values <- relevel(regressor_values, ref = reference.level)
    }
    
    cat("Categorical regressor with levels:", paste(levels(regressor_values), collapse = ", "), "\n")
    cat("Reference level:", levels(regressor_values)[1], "\n")
    
  } else if (regressor.type == "ordinal") {
    # Check if it's already ordered
    if (!is.ordered(regressor_values)) {
      if (is.factor(regressor_values)) {
        regressor_values <- ordered(regressor_values)
      } else {
        # Try to convert to ordered factor
        unique_vals <- sort(unique(regressor_values[!is.na(regressor_values)]))
        regressor_values <- ordered(regressor_values, levels = unique_vals)
      }
    }
    
    cat("Ordinal regressor with levels:", paste(levels(regressor_values), collapse = " < "), "\n")
    cat("Ordinal method:", ordinal.method, "\n")
    
  } else if (regressor.type == "continuous") {
    regressor_values <- as.numeric(regressor_values)
    if (all(is.na(regressor_values))) {
      stop("Regressor cannot be converted to numeric for continuous analysis")
    }
  }
  
  # Check for missing values in regressor
  if (any(is.na(regressor_values))) {
    warning("Missing values found in regressor variable")
    valid_cells <- !is.na(regressor_values)
    expr_data <- expr_data[, valid_cells]
    meta_data <- meta_data[valid_cells, ]
    regressor_values <- regressor_values[valid_cells]
  }
  
  # Prepare covariate formula
  covariate_formula <- ""
  if (!is.null(covariates)) {
    # Check if covariates exist in metadata
    missing_covariates <- setdiff(covariates, colnames(meta_data))
    if (length(missing_covariates) > 0) {
      warning(paste("Covariates not found:", paste(missing_covariates, collapse = ", ")))
      covariates <- intersect(covariates, colnames(meta_data))
    }
    if (length(covariates) > 0) {
      covariate_formula <- paste("+", paste(covariates, collapse = " + "))
    }
  }
  
  # Function to fit model for each gene
  fit_gene_model <- function(gene) {
    expression <- as.numeric(expr_data[gene, ])
    
    # Create data frame for modeling
    model_data <- data.frame(
      expression = expression,
      regressor = regressor_values
    )
    
    # Add covariates if specified
    if (!is.null(covariates) && length(covariates) > 0) {
      for (cov in covariates) {
        model_data[[cov]] <- meta_data[[cov]]
      }
    }
    
    # Remove rows with missing values
    model_data <- model_data[complete.cases(model_data), ]
    
    if (nrow(model_data) < 10) {
      return(data.frame(
        gene = gene,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        n_cells = nrow(model_data),
        regressor_type = regressor.type
      ))
    }
    
    # Prepare formula based on regressor type
    if (regressor.type == "continuous") {
      regressor_formula <- "regressor"
    } else if (regressor.type == "categorical") {
      regressor_formula <- "regressor"  # R automatically handles factor contrasts
    } else if (regressor.type == "ordinal") {
      if (ordinal.method == "linear") {
        # Convert to numeric for linear trend
        model_data$regressor_numeric <- as.numeric(model_data$regressor)
        regressor_formula <- "regressor_numeric"
      } else if (ordinal.method == "polynomial") {
        # Use polynomial contrasts
        regressor_formula <- "poly(as.numeric(regressor), degree = min(3, length(levels(regressor))-1))"
      } else if (ordinal.method == "spline") {
        # Use natural splines (requires splines package)
        if (requireNamespace("splines", quietly = TRUE)) {
          regressor_formula <- "splines::ns(as.numeric(regressor), df = min(3, length(levels(regressor))-1))"
        } else {
          warning("splines package not available, using linear method")
          model_data$regressor_numeric <- as.numeric(model_data$regressor)
          regressor_formula <- "regressor_numeric"
        }
      }
    }
    
    # Fit model based on effect type
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
        
        # Extract results - handle different regressor types
        model_summary <- tidy(model)
        
        if (regressor.type == "categorical") {
          # For categorical variables, extract all factor levels
          regressor_rows <- model_summary[grepl("^regressor", model_summary$term), ]
          
          if (nrow(regressor_rows) == 0) {
            # Single level case
            regressor_rows <- model_summary[model_summary$term == "regressor", ]
          }
          
          # Calculate overall F-test p-value for categorical variables
          if (nrow(regressor_rows) > 1) {
            model_anova <- anova(model)
            # More robust way to get p-value from anova table
            pval_col <- which(grepl("Pr\\(>F\\)", colnames(model_anova)))
            if (length(pval_col) > 0) {
              overall_p <- model_anova[1, pval_col]
            } else {
              overall_p <- regressor_rows$p.value[1]
            }
          } else {
            overall_p <- regressor_rows$p.value[1]
          }
          
          # Return summary statistics
          result <- data.frame(
            gene = gene,
            estimate = ifelse(nrow(regressor_rows) > 0, regressor_rows$estimate[1], NA),
            std.error = ifelse(nrow(regressor_rows) > 0, regressor_rows$std.error[1], NA),
            statistic = ifelse(nrow(regressor_rows) > 0, regressor_rows$statistic[1], NA),
            p.value = overall_p,
            n_cells = nrow(model_data),
            model_type = paste(effect, link.function, regressor.type, sep = "_"),
            n_levels = length(levels(model_data$regressor))
          )
          
        } else {
          # For continuous and ordinal (treated as continuous)
          if (regressor.type == "ordinal" && ordinal.method == "linear") {
            regressor_rows <- model_summary[model_summary$term == "regressor_numeric", ]
            overall_p <- regressor_rows$p.value[1]
          } else if (regressor.type == "ordinal" && ordinal.method %in% c("polynomial", "spline")) {
            regressor_rows <- model_summary[grepl("poly|ns", model_summary$term), ]
            # For polynomial/spline, use F-test
            if (nrow(regressor_rows) > 1) {
              model_anova <- anova(model)
              # Find the row corresponding to polynomial/spline term
              poly_spline_rows <- grepl("poly|ns", rownames(model_anova))
              pval_col <- which(grepl("Pr\\(>F\\)", colnames(model_anova)))
              
              if (any(poly_spline_rows) && length(pval_col) > 0) {
                overall_p <- model_anova[which(poly_spline_rows)[1], pval_col]
              } else {
                overall_p <- regressor_rows$p.value[1]
              }
            } else {
              overall_p <- regressor_rows$p.value[1]
            }
          } else {
            regressor_rows <- model_summary[model_summary$term == "regressor", ]
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
        }
        
      } else if (effect == "random") {
        # Random effects implementation
        if (is.null(covariates) || length(covariates) == 0) {
          stop("Random effects models require at least one grouping variable in covariates")
        }
        
        group_var <- covariates[1]
        formula_str <- paste("expression ~", regressor_formula, "+ (1|", group_var, ")")
        if (length(covariates) > 1) {
          formula_str <- paste(formula_str, "+", paste(covariates[-1], collapse = " + "))
        }
        
        model <- lmer(as.formula(formula_str), data = model_data, ...)
        model_summary <- tidy(model, effects = "fixed")
        
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
      
      return(result)
      
    }, error = function(e) {
      data.frame(
        gene = gene,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        n_cells = nrow(model_data),
        model_type = paste(effect, link.function, regressor.type, sep = "_"),
        error = as.character(e)
      )
    })
  }
  
  # Apply function to all features
  cat("Fitting models for", length(features), "features...\n")
  
  # Use parallel processing if available
  if (requireNamespace("parallel", quietly = TRUE) && length(features) > 100) {
    results <- parallel::mclapply(features, fit_gene_model, mc.cores = parallel::detectCores() - 1)
  } else {
    results <- lapply(features, fit_gene_model)
  }
  
  # Combine results
  results_df <- do.call(rbind, results)
  
  # Multiple testing correction
  valid_pvals <- !is.na(results_df$p.value)
  results_df$adj.p.value <- NA
  if (sum(valid_pvals) > 0) {
    results_df$adj.p.value[valid_pvals] <- p.adjust(results_df$p.value[valid_pvals], method = "BH")
  }
  
  # Sort by p-value
  results_df <- results_df[order(results_df$p.value, na.last = TRUE), ]
  
  # Add metadata about the analysis
  attr(results_df, "analysis_info") <- list(
    regressor = regressor,
    regressor.type = regressor.type,
    layer = layer,
    effect = effect,
    link.function = link.function,
    covariates = covariates,
    n_features_tested = length(features),
    n_cells = ncol(expr_data),
    reference.level = ifelse(regressor.type == "categorical", levels(regressor_values)[1], NA),
    ordinal.method = ifelse(regressor.type == "ordinal", ordinal.method, NA)
  )
  
  if (return.full) {
    return(list(
      results = results_df,
      seurat_object = sobj,
      features_tested = features
    ))
  } else {
    return(results_df)
  }
}

# Helper function to visualize top results (linear_seurat)
plot_top_genes <- function(results, sobj, layer = "data", top_n = 6) {
  library(ggplot2)
  library(gridExtra)
  
  # Get analysis info
  analysis_info <- attr(results, "analysis_info")
  regressor <- analysis_info$regressor
  regressor.type <- analysis_info$regressor.type
  
  # Get top significant genes
  top_genes <- head(results[!is.na(results$p.value), ], top_n)$gene
  
  # Get expression data
  if (layer == "counts") {
    expr_data <- GetAssayData(sobj, layer = "counts")
  } else if (layer == "data") {
    expr_data <- GetAssayData(sobj, layer = "data")
  } else if (layer == "scale.data") {
    expr_data <- GetAssayData(sobj, layer = "scale.data")
  }
  
  # Prepare regressor data
  regressor_data <- sobj@meta.data[[regressor]]
  if (regressor.type == "categorical") {
    regressor_data <- as.factor(regressor_data)
  } else if (regressor.type == "ordinal") {
    if (!is.ordered(regressor_data)) {
      if (is.factor(regressor_data)) {
        regressor_data <- ordered(regressor_data)
      } else {
        unique_vals <- sort(unique(regressor_data[!is.na(regressor_data)]))
        regressor_data <- ordered(regressor_data, levels = unique_vals)
      }
    }
  } else {
    regressor_data <- as.numeric(regressor_data)
  }
  
  # Create plots
  plots <- list()
  for (i in seq_along(top_genes)) {
    gene <- top_genes[i]
    gene_result <- results[results$gene == gene, ]
    
    plot_data <- data.frame(
      expression = as.numeric(expr_data[gene, ]),
      regressor = regressor_data
    )
    
    # Remove missing values
    plot_data <- plot_data[complete.cases(plot_data), ]
    
    if (regressor.type == "categorical") {
      p <- ggplot(plot_data, aes(x = regressor, y = expression)) +
        geom_boxplot(aes(fill = regressor), alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(
          title = paste0(gene, " (p=", format(gene_result$p.value, digits = 3), ")"),
          x = regressor,
          y = paste("Expression (", layer, ")", sep = "")
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = "none")
      
    } else if (regressor.type == "ordinal") {
      p <- ggplot(plot_data, aes(x = regressor, y = expression)) +
        geom_boxplot(aes(group = regressor, fill = regressor), alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red") +
        labs(
          title = paste0(gene, " (p=", format(gene_result$p.value, digits = 3), ")"),
          x = regressor,
          y = paste("Expression (", layer, ")", sep = "")
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = "none")
      
    } else {  # continuous
      p <- ggplot(plot_data, aes(x = regressor, y = expression)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(
          title = paste0(gene, " (p=", format(gene_result$p.value, digits = 3), ")"),
          x = regressor,
          y = paste("Expression (", layer, ")", sep = "")
        ) +
        theme_minimal()
    }
    
    plots[[i]] <- p
  }
  
  # Arrange plots
  do.call(grid.arrange, c(plots, ncol = 2))
}

# Example usage function
example_usage <- function() {
  # === CATEGORICAL VARIABLES ===
  # results_cat <- linear_seurat(sobj, 
  #                             regressor = "cell_type",
  #                             regressor.type = "categorical",
  #                             reference.level = "Control")
  
  # === ORDINAL VARIABLES ===
  # results_ord <- linear_seurat(sobj,
  #                             regressor = "disease_stage",
  #                             regressor.type = "ordinal",
  #                             ordinal.method = "linear")
  
  # === CONTINUOUS VARIABLES ===
  # results_cont <- linear_seurat(sobj,
  #                              regressor = "age",
  #                              regressor.type = "continuous")
}


# analysis -----
NEBULA=function(){}


#' @title MAST 파이프라인 함수
#'
#' @description [수정] Seurat -> SCE로 변환하여 MAST를 실행
#'
#' @param sobj (Seurat) Seurat 객체
#' @param formula (formula or character) lme4 문법의 포뮬러
#' @param min_cells_expr (numeric) 유전자 필터링 기준 (최소 발현 세포 수)
#' @param n_cores (numeric) 병렬 처리 코어 수
#' @param lrt_variable (character) LRT 검정을 수행할 변수명 (예: "type")
#'
#' @export
runMAST <- function(sobj,
                    formula,
                    min_cells_expr = 10,
                    n_cores = 4,
                    lrt_variable = NULL) {
  
  if (!requireNamespace("MAST", quietly = TRUE)) stop("MAST 패키지가 필요합니다.")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("SCE 패키지가 필요합니다.")
  if (is.null(lrt_variable)) stop("'lrt_variable' 인자 (예: 'g3')를 지정해야 합니다.")
  
  if (is.character(formula)) {
    formula_obj <- as.formula(formula)
  } else if (inherits(formula, "formula")) {
    formula_obj <- formula
  } else {
    stop("'formula'는 문자열 또는 formula 객체여야 합니다.")
  }
  
  # --- 1. SCA 객체 생성 및 정규화 ---
  message("1/5: Seurat -> SingleCellExperiment(SCE) 객체 변환 중...")
  
  # [수정] MAST Problem 1 해결: FromSeurat 대신 as.SingleCellExperiment 사용
  sca <- as.SingleCellExperiment(sobj)
  
  # --- 2. 유전자 필터링 ---
  message(sprintf("2/5: 유전자 필터링 (min %d cells)...", min_cells_expr))
  # freq()는 발현 비율 (0~1)을 반환
  keep_genes <- (MAST::freq(sca) * ncol(sca)) >= min_cells_expr
  sca_filtered <- sca[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(sca)))

  # --- 3. 정규화 (MAST는 log2cpm 사용) ---
  message("3/5: Log2(CPM+1) 정규화 중...")
  SummarizedExperiment::assay(sca_filtered, "logcpm") <- MAST::cpm(sca_filtered, log = TRUE)
  
  # --- 4. zlm (Hurdle LMM) 실행 ---
  message(sprintf("4/5: MAST::zlm 실행 (Cores: %d). 시간이 오래 걸릴 수 있습니다...", n_cores))
  
  zfit <- MAST::zlm(formula_obj, 
                    sca = sca_filtered, 
                    method = "glmer", 
                    parallel = TRUE,
                    nCores = n_cores)
  
  # --- 5. LRT 결과 요약 ---
  message(sprintf("5/5: LRT 검정 수행 (변수: %s)...", lrt_variable))
  summary_res <- summary(zfit, doLRT = lrt_variable)
  summary_dt <- summary_res$datatable
  
  results_df <- merge(
      summary_dt[component == 'H', .(primerid, `Pr(>Chisq)`)], # Hurdle (logistic)
      summary_dt[component == 'logcpm', .(primerid, coef, ci.hi, ci.lo)], # Continuous
      by = 'primerid'
  )
  colnames(results_df)[2] <- "p_value_hurdle"
  results_df <- results_df[order(p_value_hurdle), ]
  
  message("MAST 분석 완료.")
  return(results_df)
}

# [V4] find_gene_signature_v4
#
# 변경점:
# 1. control_vars = NULL 인자 추가 (교란 변수 보정)
# 2. test_n = NULL 인자 추가 (속도 향상을 위한 사전 필터링)
# 3. limma, lasso, tree_based, gam: Native 방식 보정
# 4. wilcoxon, nmf, pca_loadings: limma::removeBatchEffect 방식 보정

#' @export
find_gene_signature_v4 <- function(data, 
                                 meta.data = NULL,
                                 target_var,
                                 target_group = NULL,
                                 control_vars = NULL,   # <<< [V4 UPGRADE]
                                 method = c("tree_based", "lasso", "limma", 
                                            "nmf", "wilcoxon", "gam", "pca_loadings"),
                                 n_features = 50,
                                 test_n = NULL,         # <<< [V4 UPGRADE]
                                 preprocess = TRUE,
                                 min_cells = 10,
                                 min_pct = 0.01,
                                 return_model = FALSE,
                                 fgs_seed = 42,
                                 lambda_selection = "lambda.1se",
                                 ...) {
  
  # ===
  # 0. Batch Processing (v3와 동일)
  # ===
  
  all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
  if (is.null(method)) {
    method <- all_methods
  }
  
  if (length(method) > 1) {
    current_call <- match.call()
    
    results_list <- lapply(method, function(m) {
      single_call <- current_call
      single_call$method <- m
      tryCatch({
        eval(single_call)
      }, error = function(e) {
        warning(sprintf("Method '%s' failed with error: %s", m, e$message))
        return(list(method = m, error = e$message))
      })
    })
    
    names(results_list) <- method
    return(results_list)
  }
  
  if (!method %in% all_methods) {
    stop(sprintf("Invalid method '%s'. Choose from: %s", 
                 method, paste(all_methods, collapse=", ")))
  }
  
  set.seed(fgs_seed)
  
  # ===
  # 1. Input validation and data extraction (v3와 동일)
  # ===
  
  is_seurat <- inherits(data, "Seurat")
  
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required but not installed")
    }
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    if ("data" %in% names(data@assays[[1]])) {
      expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
    } else {
      expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
    }
  } else {
    if (is.null(meta.data)) {
      stop("meta.data must be provided when data is not a Seurat object")
    }
    expr_mat <- as.matrix(data)
    if (nrow(expr_mat) == nrow(meta.data)) {
      expr_mat <- t(expr_mat)
    }
  }
  
  if (!target_var %in% colnames(meta.data)) {
    stop(sprintf("target_var '%s' not found", target_var))
  }
  
  # [V4 UPGRADE] control_vars 유효성 검사
  if (!is.null(control_vars)) {
    missing_vars <- control_vars[!control_vars %in% colnames(meta.data)]
    if (length(missing_vars) > 0) {
      stop(sprintf("control_vars not found in metadata: %s", 
                   paste(missing_vars, collapse=", ")))
    }
  }
  
  common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
  if (length(common_cells) == 0) {
    stop("No common cells found between expression matrix and metadata")
  }
  expr_mat <- expr_mat[, common_cells]
  meta.data <- meta.data[common_cells, ]
  
  # ===
  # 2. Process target variable (v3와 동일)
  # ===
  
  target_values <- meta.data[[target_var]]
  target_type <- class(target_values)[1]
  
  if (is.numeric(target_values)) {
    # ... (numeric 처리 로직, v3와 동일) ...
     if (is.null(target_group)) { target_group <- 0.25 }
     if (is.list(target_group)) {
       low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
       high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
     } else if (length(target_group) == 1 && target_group < 1) {
       low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
       high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
     } else {
       low_cutoff <- high_cutoff <- target_group
     }
     group_labels <- ifelse(target_values <= low_cutoff, "Low",
                            ifelse(target_values >= high_cutoff, "High", NA))
     keep_cells <- !is.na(group_labels)
  } else {
    if (!is.null(target_group)) {
      keep_cells <- target_values %in% target_group
    } else {
      keep_cells <- TRUE # 모든 팩터 레벨 사용
    }
    group_labels <- target_values
  }
  
  expr_mat <- expr_mat[, keep_cells]
  meta.data <- meta.data[keep_cells, ]
  target_binary <- factor(group_labels[keep_cells]) # factor로 변환
  
  if (length(unique(target_binary)) < 2) {
    stop("Target variable must have at least 2 groups after processing")
  }
  
  n_groups <- length(unique(target_binary))
  
  # ===
  # 3. Filter and preprocess genes (v3와 동일)
  # ===
  
  n_cells_expr <- rowSums(expr_mat > 0)
  pct_cells_expr <- n_cells_expr / ncol(expr_mat)
  keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  expr_mat <- expr_mat[keep_genes, ]
  
  if (nrow(expr_mat) == 0) {
    stop("No genes pass filtering criteria")
  }
  
  # ===
  # 3.5 [V4 UPGRADE] Speed-up: Pre-filter genes using test_n
  # ===
  if (!is.null(test_n) && nrow(expr_mat) > test_n) {
    message(sprintf("Pre-filtering to top %d genes based on limma p-value...", test_n))
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for 'test_n' pre-filtering")
    }
    
    # limma는 교란변수를 '미리' 보정하고 p-value를 뽑을 수 있음
    if (is.null(control_vars)) {
      design_test <- model.matrix(~ target_binary, data = meta.data)
    } else {
      formula_test <- as.formula(paste("~ target_binary +", paste(control_vars, collapse="+")))
      design_test <- model.matrix(formula_test, data = meta.data)
    }
    
    fit_test <- limma::lmFit(expr_mat, design_test)
    fit_test <- limma::eBayes(fit_test)
    
    # p-value 추출 (target_binary 그룹 비교)
    top_table_test <- limma::topTable(fit_test, 
                                      coef = grep("target_binary", colnames(design_test)), 
                                      number = Inf, sort.by = "P")
    
    top_gene_names <- rownames(top_table_test)[1:min(test_n, nrow(top_table_test))]
    expr_mat <- expr_mat[top_gene_names, ]
    message(sprintf("... reduced to %d genes.", nrow(expr_mat)))
  }
  
  if (preprocess) {
    if (max(expr_mat) > 100) {
      expr_mat <- log1p(expr_mat)
    }
    if (method %in% c("lasso", "gam", "pca_loadings")) {
      gene_means <- rowMeans(expr_mat)
      gene_sds <- apply(expr_mat, 1, sd)
      gene_sds[gene_sds == 0] <- 1
      expr_mat <- (expr_mat - gene_means) / gene_sds
    }
  }
  
  # ===
  # 3.8 [V4 UPGRADE] Confounder pre-correction (for non-native methods)
  # ===
  if (!is.null(control_vars) && method %in% c("wilcoxon", "nmf", "pca_loadings")) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for removeBatchEffect")
    }
    
    message(sprintf("Method %s: Applying limma::removeBatchEffect for: %s",
                    method, paste(control_vars, collapse=", ")))
    
    # removeBatchEffect는 design matrix가 아닌 'batch' 또는 'covariates'를 받음
    covariates_df <- meta.data[, control_vars, drop = FALSE]
    
    # Factor 변수들을 model.matrix로 변환
    covariate_mat <- model.matrix(~ . - 1, data = covariates_df) 
    
    expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
  }
  
  
  # ===
  # 4. Method-specific signature discovery
  # ===
  
  # X (samples x genes)는 여러 메서드에서 공통으로 사용
  X <- t(expr_mat)
  y <- target_binary
  
  # [V4 UPGRADE] Native-supporting methods에서 사용할 교란변수 매트릭스
  if (!is.null(control_vars) && method %in% c("tree_based", "lasso", "gam")) {
     covariate_mat_model <- model.matrix(~ . - 1, 
                                data = meta.data[, control_vars, drop=FALSE])
  }

  
  result <- switch(method,
                 
                 tree_based = {
                   if (!requireNamespace("randomForest", quietly = TRUE)) { ... }
                   
                   X_rf <- X
                   
                   # [V4 UPGRADE] 교란 변수를 X 매트릭스에 추가
                   if (!is.null(control_vars)) {
                     X_rf <- cbind(X_rf, covariate_mat_model)
                   }
                   
                   if (ncol(X_rf) > 2000) {
                     # ... (v3의 top_var_genes 로직은 test_n으로 대체됨) ...
                     # ... 단, control_vars가 추가되었으므로 컬럼 이름으로 필터링
                     gene_vars <- apply(X, 2, var) # 원본 X에서만 var 계산
                     top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
                     
                     # 교란 변수는 항상 포함
                     X_rf <- X_rf[, c(top_var_genes, colnames(covariate_mat_model))]
                   }
                   
                   rf_model <- randomForest::randomForest(x = X_rf, y = y, ntree = 500, importance = TRUE, ...)
                   
                   importance_scores <- randomForest::importance(rf_model)
                   
                   # ... (v3의 MeanDecreaseGini 로직) ...
                   if (n_groups == 2) {
                     weights_magnitude_all <- importance_scores[, "MeanDecreaseGini"]
                   } else {
                     weights_magnitude_all <- rowMeans(importance_scores[, grep("MeanDecreaseGini", colnames(importance_scores))])
                   }
                   
                   # [V4 UPGRADE] 교란 변수를 제외하고 '유전자'만 선택
                   gene_names_in_model <- colnames(X_rf)[!colnames(X_rf) %in% colnames(covariate_mat_model)]
                   weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]
                   
                   top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
                   weights_magnitude <- weights_magnitude_genes[top_genes]
                   
                   # ... (v3의 방향성 보정 및 나머지 로직 동일) ...
                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]; g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1 
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     weights <- weights_magnitude
                   }
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   # ... (v3의 performance 로직) ...
                   list(...) 
                 },
                 
                 lasso = {
                   if (!requireNamespace("glmnet", quietly = TRUE)) { ... }
                   
                   # [V4 UPGRADE] X 매트릭스 구성 및 penalty.factor 설정
                   if (is.null(control_vars)) {
                     X_model <- X
                     penalty_vec <- rep(1, ncol(X_model)) # 모든 유전자에 페널티 1
                   } else {
                     X_model <- cbind(X, covariate_mat_model)
                     # 유전자는 페널티 1, 교란 변수는 페널티 0 (Lasso에서 제외)
                     penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
                   }
                   
                   # [V4 UPGRADE] cv.glmnet 호출 시 penalty.factor 전달
                   if (n_groups == 2) {
                     cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=1, 
                                                 penalty.factor = penalty_vec, ...)
                   } else {
                     cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=1, 
                                                 penalty.factor = penalty_vec, ...)
                   }
                   
                   coefs <- coef(cv_fit, s = lambda_selection) 
                   
                   # [V4 UPGRADE] weights 추출 시 유전자 부분만 선택
                   if (n_groups == 2) {
                     # 1번은 intercept, (ncol(X)+1)까지가 유전자
                     weights_all <- as.numeric(coefs[2:(ncol(X)+1)]) 
                     names(weights_all) <- rownames(coefs)[2:(ncol(X)+1)]
                   } else {
                     coef_list <- lapply(coefs, function(x) as.numeric(x[2:(ncol(X)+1)]))
                     weights_all <- rowMeans(do.call(cbind, coef_list))
                     names(weights_all) <- rownames(coefs[[1]])[2:(ncol(X)+1)]
                   }
                   
                   # ... (v3의 non-zero 유전자 선택 및 나머지 로직 동일) ...
                   nonzero_genes <- names(weights_all)[weights_all != 0]
                   if (length(nonzero_genes) == 0) {
                     top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
                   } else {
                     nonzero_weights <- weights_all[nonzero_genes]
                     top_genes <- names(sort(abs(nonzero_weights), decreasing=TRUE)[1:min(n_features, length(nonzero_weights))])
                   }
                   weights <- weights_all[top_genes]
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   # ... (v3의 performance 로직) ...
                   list(...)
                 },
                 
                 limma = {
                   if (!requireNamespace("limma", quietly = TRUE)) { ... }
                   
                   # [V4 UPGRADE] control_vars를 포함하는 formula 생성
                   if (is.null(control_vars)) {
                     design <- model.matrix(~0 + target_binary, data = meta.data)
                     colnames(design)[1:n_groups] <- levels(target_binary)
                   } else {
                     control_formula <- paste(control_vars, collapse = " + ")
                     full_formula <- as.formula(paste("~0 + target_binary +", control_formula))
                     design <- model.matrix(full_formula, data = meta.data)
                     colnames(design)[1:n_groups] <- levels(target_binary) 
                   }

                   fit <- limma::lmFit(expr_mat, design)
                   
                   # ... (v3의 contrast 생성 및 나머지 로직 동일) ...
                   if (n_groups == 2) {
                     contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
                   } else {
                      # ... (v3의 n_groups > 2 로직) ...
                   }
                   contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
                   fit2 <- limma::contrasts.fit(fit, contrast_mat)
                   fit2 <- limma::eBayes(fit2)
                   top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
                   
                   weights_all <- top_table$t
                   names(weights_all) <- rownames(top_table)
                   top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
                   weights <- weights_all[top_genes]
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = list(top_table = top_table[top_genes, ]), 
                        model = if(return_model) fit2 else NULL)
                 },
                 
                 wilcoxon = {
                   # [V4] 이 블록은 수정할 필요가 없습니다.
                   # 만약 control_vars가 있었다면, 함수 상단의
                   # '3.8 Pre-correction' 블록에서 expr_mat이 이미 보정되었습니다.
                   
                   # ... (v3의 wilcoxon 로직 동일) ...
                   pvals <- numeric(nrow(expr_mat))
                   effect_sizes <- numeric(nrow(expr_mat))
                   for (i in 1:nrow(expr_mat)) { ... }
                   # ... (v3 로직 계속) ...
                   list(...)
                 },
                 
                 nmf = {
                   # [V4] 이 블록도 수정할 필요가 없습니다.
                   # control_vars가 있었다면 expr_mat이 이미 보정되었습니다.
                   
                   # ... (v3의 nmf 로직 동일) ...
                   list(...)
                 },
                 
                 gam = {
                   if (!requireNamespace("mgcv", quietly = TRUE)) { ... }
                   
                   # ... (v3의 n_test_genes 로직 동일) ...
                   
                   deviance_explained <- numeric(nrow(expr_mat))
                   names(deviance_explained) <- rownames(expr_mat)
                   
                   # [V4 UPGRADE] 교란 변수를 포함하는 GAM formula 생성
                   if (is.null(control_vars)) {
                     formula_base <- ""
                     gam_data <- data.frame(y_var = target_binary)
                   } else {
                     formula_base <- paste(" +", paste(control_vars, collapse=" + "))
                     gam_data <- data.frame(y_var = target_binary, 
                                            meta.data[, control_vars, drop=FALSE])
                   }
                   
                   # ... (v3의 gam 루프) ...
                   for (i in test_genes_idx) {
                     gam_data$gene_expr <- expr_mat[i, ]
                     
                     if (n_groups == 2) {
                       gam_data$y_var_numeric <- as.numeric(gam_data$y_var) - 1
                       full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
                       gam_fit <- mgcv::gam(as.formula(full_formula_str), 
                                            data = gam_data, family="binomial", ...)
                     } else {
                       gam_data$y_var_numeric <- as.numeric(gam_data$y_var)
                       full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
                       gam_fit <- mgcv::gam(as.formula(full_formula_str), 
                                            data = gam_data, ...)
                     }
                     deviance_explained[i] <- summary(gam_fit)$dev.expl
                   }
                   
                   # ... (v3의 나머지 로직 동일) ...
                   weights_magnitude <- deviance_explained
                   top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0))])
                   # ... (방향성 보정) ...
                   list(...)
                 },
                 
                 pca_loadings = {
                   # [V4] 이 블록도 수정할 필요가 없습니다.
                   # control_vars가 있었다면 expr_mat이 이미 보정되었습니다.
                   
                   # ... (v3의 pca_loadings 로직 동일) ...
                   list(...)
                 }
  )
  
  # ===
  # 5. Return results (v3와 동일)
  # ===
  
  result$method <- method
  result$target_var <- target_var
  result$n_groups <- n_groups
  result$n_cells <- ncol(expr_mat)
  result$formula <- paste(deparse(match.call()), collapse = " ")
  
  class(result) <- c("gene_signature", "list")
  return(result)
}

# [V(4.1)] find_gene_signature_v4.1
#
# 변경점 (2025-11-10):
# 1. [FIX-NA] 2.Process target: meta.data 서브셋 후 'droplevels()'를 호출하여
#    control_vars/target_binary의 미사용 레벨을 제거 (model.matrix NA 에러 방지)
# 2. [FIX-Lasso] 4.Method (lasso): 'nonzero_genes' 필터링 로직 제거.
#    항상 'abs(weights_all)' 기준 상위 n_features 반환 (n_features 개수 보장)
# 3. [FIX-NMF] 3.8 Pre-correction: 'nmf'를 'removeBatchEffect' 대상에서 제외.
# 4. [FIX-NMF] 4.Method (nmf): 'nmf' 블록 *내부*에서 control_vars를
#    'removeBatchEffect'로 보정하고, 즉시 'min-shift'로 양수화.
# 5. [FIX-GAM] 4.Method (gam): 'gam' 블록 내부의 중복된 유전자 필터링
#    ('test_genes_idx') 로직 제거. 'test_n'으로 필터링된 모든 유전자 사용.
# 6. [FIX-Wilcox] 4.Method (wilcoxon): 'wilcox.test'/'kruskal.test'에서 '...' 인자 제거.
#' @export
find_gene_signature_v4.1 <- function(data, 
                                 meta.data = NULL,
                                 target_var,
                                 target_group = NULL,
                                 control_vars = NULL,   
                                 method = c("tree_based", "lasso", "limma", 
                                            "nmf", "wilcoxon", "gam", "pca_loadings"),
                                 n_features = 50,
                                 test_n = NULL,         
                                 preprocess = TRUE,
                                 min_cells = 10,
                                 min_pct = 0.01,
                                 return_model = FALSE,
                                 fgs_seed = 42,
                                 lambda_selection = "lambda.1se",
                                 ...) {
  
  # ===
  # 0. Batch Processing (v4와 동일)
  # ===
  
  all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
  if (is.null(method)) {
    method <- all_methods
  }
  
  if (length(method) > 1) {
    current_call <- match.call()
    
    results_list <- lapply(method, function(m) {
      single_call <- current_call
      single_call$method <- m
      tryCatch({
        eval(single_call)
      }, error = function(e) {
        warning(sprintf("Method '%s' failed with error: %s", m, e$message))
        return(list(method = m, error = e$message))
      })
    })
    
    names(results_list) <- method
    return(results_list)
  }
  
  if (!method %in% all_methods) {
    stop(sprintf("Invalid method '%s'. Choose from: %s", 
                 method, paste(all_methods, collapse=", ")))
  }
  
  set.seed(fgs_seed)
  
  # ===
  # 1. Input validation and data extraction (v4와 동일)
  # ===
  
  is_seurat <- inherits(data, "Seurat")
  
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required but not installed")
    }
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    # Seurat v5/v4 호환
    if ("data" %in% Seurat::Layers(data)) {
        expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
    } else if ("counts" %in% Seurat::Layers(data)) {
        expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
    } else if ("data" %in% names(data@assays[[Seurat::DefaultAssay(data)]])) {
        expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
    } else {
        expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
    }
  } else {
    if (is.null(meta.data)) {
      stop("meta.data must be provided when data is not a Seurat object")
    }
    expr_mat <- as.matrix(data)
    if (nrow(expr_mat) == nrow(meta.data)) {
      expr_mat <- t(expr_mat)
    }
  }
  
  if (!target_var %in% colnames(meta.data)) {
    stop(sprintf("target_var '%s' not found in metadata columns: %s", 
                 target_var, paste(colnames(meta.data), collapse=", ")))
  }
  
  if (!is.null(control_vars)) {
    missing_vars <- control_vars[!control_vars %in% colnames(meta.data)]
    if (length(missing_vars) > 0) {
      stop(sprintf("control_vars not found in metadata: %s", 
                   paste(missing_vars, collapse=", ")))
    }
  }
  
  common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
  if (length(common_cells) == 0) {
    stop("No common cells found between expression matrix and metadata")
  }
  expr_mat <- expr_mat[, common_cells]
  meta.data <- meta.data[common_cells, ]
  
  # ===
  # 2. Process target variable (v4.1 수정)
  # ===
  
  target_values <- meta.data[[target_var]]
  target_type <- class(target_values)[1]
  
  if (is.numeric(target_values)) {
    if (is.null(target_group)) {
      target_group <- 0.25
    }
    if (is.list(target_group)) {
      low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
      high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
      group_labels <- ifelse(target_values <= low_cutoff, "Low",
                             ifelse(target_values >= high_cutoff, "High", NA))
    } else if (length(target_group) == 1 && target_group < 1) {
      low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
      high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
      group_labels <- ifelse(target_values <= low_cutoff, "Low",
                             ifelse(target_values >= high_cutoff, "High", NA))
    } else {
      group_labels <- ifelse(target_values < target_group, "Low", "High")
    }
    keep_cells <- !is.na(group_labels)
  } else {
    if (!is.null(target_group)) {
      keep_cells <- target_values %in% target_group
    } else {
      keep_cells <- rep(TRUE, length(target_values))
    }
    group_labels <- target_values
    keep_cells <- keep_cells & !is.na(group_labels)
  }
  
  expr_mat <- expr_mat[, keep_cells]
  meta.data <- meta.data[keep_cells, ]
  
  # [FIX-NA] 사용하지 않는 팩터 레벨 제거 (model.matrix 에러 방지)
  target_binary <- factor(group_labels[keep_cells])
  
  if (!is.null(control_vars)) {
    for (cv in control_vars) {
      if (is.factor(meta.data[[cv]])) {
        meta.data[[cv]] <- droplevels(meta.data[[cv]])
      }
    }
  }

  if (length(unique(target_binary)) < 2) {
    stop("Target variable must have at least 2 groups after processing")
  }
  n_groups <- length(unique(target_binary))
  
  # ===
  # 3. Filter and preprocess genes (v4와 동일)
  # ===
  
  n_cells_expr <- rowSums(expr_mat > 0)
  pct_cells_expr <- n_cells_expr / ncol(expr_mat)
  keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  expr_mat <- expr_mat[keep_genes, ]
  
  if (nrow(expr_mat) == 0) {
    stop("No genes pass filtering criteria")
  }
  
  # ===
  # 3.5 [V4] Speed-up: Pre-filter genes using test_n
  # ===
  if (!is.null(test_n) && nrow(expr_mat) > test_n) {
    message(sprintf("Pre-filtering to top %d genes based on limma p-value...", test_n))
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for 'test_n' pre-filtering")
    }
    
    if (is.null(control_vars)) {
      design_test <- model.matrix(~ target_binary, data = meta.data)
    } else {
      formula_test <- as.formula(paste("~ target_binary +", paste(control_vars, collapse="+")))
      design_test <- model.matrix(formula_test, data = meta.data)
    }
    
    fit_test <- limma::lmFit(expr_mat, design_test)
    fit_test <- limma::eBayes(fit_test)
    
    coef_indices <- grep("target_binary", colnames(design_test))
    if (length(coef_indices) == 0) {
        warning("Could not find target_binary coefficients for test_n pre-filtering.")
        # Fallback: use all coefficients (less ideal)
        coef_indices <- 1:ncol(design_test)
    }
    
    top_table_test <- limma::topTable(fit_test, 
                                      coef = coef_indices, 
                                      number = Inf, sort.by = "P")
    
    top_gene_names <- rownames(top_table_test)[1:min(test_n, nrow(top_table_test))]
    expr_mat <- expr_mat[top_gene_names, ]
    message(sprintf("... reduced to %d genes.", nrow(expr_mat)))
  }
  
  if (preprocess) {
    if (max(expr_mat) > 100) {
      expr_mat <- log1p(expr_mat)
    }
    if (method %in% c("lasso", "gam", "pca_loadings")) {
      gene_means <- rowMeans(expr_mat)
      gene_sds <- apply(expr_mat, 1, sd)
      gene_sds[gene_sds == 0] <- 1
      expr_mat <- (expr_mat - gene_means) / gene_sds
    }
  }
  
  # ===
  # 3.8 [V4.1] Confounder pre-correction (non-native methods)
  # ===
  
  # [FIX-NMF] 'nmf'는 removeBatchEffect 대상에서 제외 (음수 값 방지)
  if (!is.null(control_vars) && method %in% c("wilcoxon", "pca_loadings")) { 
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for removeBatchEffect")
    }
    
    message(sprintf("Method %s: Applying limma::removeBatchEffect for: %s",
                    method, paste(control_vars, collapse=", ")))
    
    covariates_df <- meta.data[, control_vars, drop = FALSE]
    
    # Check for factors and create model matrix
    if(any(sapply(covariates_df, is.factor)) || any(sapply(covariates_df, is.character))) {
        covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
    } else {
        covariate_mat <- as.matrix(covariates_df)
    }
    
    expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
  }
  
  
  # ===
  # 4. Method-specific signature discovery (v4.1 수정)
  # ===
  
  X <- t(expr_mat)
  y <- target_binary
  
  covariate_mat_model <- NULL # Initialize
  if (!is.null(control_vars) && method %in% c("tree_based", "lasso", "gam")) {
     covariates_df_model <- meta.data[, control_vars, drop = FALSE]
     if(any(sapply(covariates_df_model, is.factor)) || any(sapply(covariates_df_model, is.character))) {
        covariate_mat_model <- model.matrix(~ . - 1, data = covariates_df_model)
     } else {
        covariate_mat_model <- as.matrix(covariates_df_model)
     }
  }

  
  result <- switch(method,
                 
                 tree_based = {
                   if (!requireNamespace("randomForest", quietly = TRUE)) {
                     stop("randomForest package required. Install with: install.packages('randomForest')")
                   }
                   
                   X_rf <- X
                   
                   if (!is.null(control_vars)) {
                     X_rf <- cbind(X_rf, covariate_mat_model)
                   }
                   
                   # Pre-filter genes if matrix is still too large
                   if (ncol(X_rf) > 2000) {
                     gene_vars <- apply(X, 2, var) # var from original X
                     top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
                     
                     control_var_names <- if(is.null(control_vars)) character(0) else colnames(covariate_mat_model)
                     X_rf <- X_rf[, c(top_var_genes, control_var_names)]
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
                   
                   top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
                   weights_magnitude <- weights_magnitude_genes[top_genes]
                   
                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]
                     g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1 
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     warning("tree_based: n_groups > 2. Score represents magnitude (importance), not direction.")
                     weights <- weights_magnitude
                   }
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   pred <- rf_model$predicted
                   if (n_groups == 2) {
                     if (requireNamespace("pROC", quietly = TRUE)) {
                       roc_obj <- pROC::roc(y, scores, quiet=TRUE)
                       auc <- as.numeric(pROC::auc(roc_obj))
                     } else { auc <- NA }
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
                   } else {
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, confusion = table(pred, y))
                   }
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) rf_model else NULL)
                 },
                 
                 lasso = {
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
                   
                   if (n_groups == 2) {
                     cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=1, 
                                                 penalty.factor = penalty_vec, ...)
                   } else {
                     cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=1, 
                                                 penalty.factor = penalty_vec, ...)
                   }
                   
                   coefs <- coef(cv_fit, s = lambda_selection) 
                   
                   if (n_groups == 2) {
                     weights_all <- as.numeric(coefs[2:(ncol(X)+1)]) 
                     names(weights_all) <- rownames(coefs)[2:(ncol(X)+1)]
                   } else {
                     coef_list <- lapply(coefs, function(x) as.numeric(x[2:(ncol(X)+1)]))
                     weights_all <- rowMeans(do.call(cbind, coef_list))
                     names(weights_all) <- rownames(coefs[[1]])[2:(ncol(X)+1)]
                   }
                   
                   if (length(weights_all) == 0) {
                     stop("LASSO returned no gene coefficients.")
                   }

                   # [FIX-Lasso] 'nonzero' 필터링 제거. 항상 abs(weight)로 정렬.
                   top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
                   
                   weights <- weights_all[top_genes]
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   pred_probs <- predict(cv_fit, newx=X_model, s = lambda_selection, type="response")
                   if (n_groups == 2) {
                     pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
                     if (requireNamespace("pROC", quietly = TRUE)) {
                       roc_obj <- pROC::roc(y, scores, quiet=TRUE) 
                       auc <- as.numeric(pROC::auc(roc_obj))
                     } else { auc <- NA }
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
                   } else {
                     pred_class_indices <- apply(pred_probs[,,1], 1, which.max)
                     pred <- levels(y)[pred_class_indices]
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, confusion = table(pred, y))
                   }
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) cv_fit else NULL)
                 },
                 
                 limma = {
                   if (!requireNamespace("limma", quietly = TRUE)) {
                     stop("limma package required. Install with: BiocManager::install('limma')")
                   }
                   
                   if (is.null(control_vars)) {
                     design <- model.matrix(~0 + target_binary, data = meta.data)
                     colnames(design)[1:n_groups] <- levels(target_binary)
                   } else {
                     control_formula <- paste(control_vars, collapse = " + ")
                     full_formula <- as.formula(paste("~0 + target_binary +", control_formula))
                     design <- model.matrix(full_formula, data = meta.data)
                     colnames(design)[1:n_groups] <- levels(target_binary) 
                   }

                   fit <- limma::lmFit(expr_mat, design)
                   
                   if (n_groups == 2) {
                     contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
                   } else {
                     # Create all pairwise contrasts if n_groups > 2
                     contrast_pairs <- combn(levels(target_binary), 2, function(x) paste(x[2], x[1], sep="-"))
                     contrast_str <- contrast_pairs
                     warning("limma: n_groups > 2. Using all pairwise contrasts. Weights based on average t-statistic.")
                   }
                   
                   contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
                   fit2 <- limma::contrasts.fit(fit, contrast_mat)
                   fit2 <- limma::eBayes(fit2)
                   
                   if(n_groups == 2) {
                       top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
                       weights_all <- top_table$t
                   } else {
                       # For multiple contrasts, sort by F-statistic
                       top_table <- limma::topTable(fit2, number=Inf, sort.by="F")
                       # Use average absolute t-stat as weight magnitude
                       weights_all <- rowMeans(abs(fit2$t[, contrast_str, drop=FALSE]))
                   }
                   
                   names(weights_all) <- rownames(top_table)
                   
                   top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
                   
                   # Re-fetch signed weights for top genes
                   if(n_groups == 2) {
                       weights <- top_table[top_genes, "t"]
                   } else {
                       # For n_groups > 2, weights are just magnitude (no clear direction)
                       weights <- weights_all[top_genes]
                   }

                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(top_table = top_table[top_genes, ])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) fit2 else NULL)
                 },
                 
                 wilcoxon = {
                   # [V4.1] control_vars가 있었다면 3.8에서 expr_mat이 이미 보정됨.
                   
                   pvals <- numeric(nrow(expr_mat))
                   effect_sizes <- numeric(nrow(expr_mat))
                   names(pvals) <- rownames(expr_mat)
                   names(effect_sizes) <- rownames(expr_mat)

                   for (i in 1:nrow(expr_mat)) {
                     if (n_groups == 2) {
                       group1 <- expr_mat[i, target_binary == levels(target_binary)[1]]
                       group2 <- expr_mat[i, target_binary == levels(target_binary)[2]]
                       
                       # [FIX-Wilcox] '...' 인자 제거
                       test <- try(wilcox.test(group1, group2), silent=TRUE) 
                       
                       if(inherits(test, "try-error")) {
                         pvals[i] <- 1.0
                         effect_sizes[i] <- 0
                       } else {
                         pvals[i] <- test$p.value
                         effect_sizes[i] <- median(group2) - median(group1)
                       }
                     } else {
                       # [FIX-Wilcox] '...' 인자 제거
                       test <- try(kruskal.test(expr_mat[i, ] ~ target_binary), silent=TRUE) 
                       
                       if(inherits(test, "try-error")) {
                          pvals[i] <- 1.0
                          effect_sizes[i] <- 0
                       } else {
                          pvals[i] <- test$p.value
                          # Use variance of medians as effect size (magnitude only)
                          effect_sizes[i] <- var(tapply(expr_mat[i, ], target_binary, median))
                       }
                     }
                   }
                   
                   padj <- p.adjust(pvals, method="BH")
                   
                   ranks <- rank(pvals) + rank(-abs(effect_sizes))
                   top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
                   
                   weights <- effect_sizes[top_genes]
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(pvalues = pvals[top_genes], padj = padj[top_genes],
                                effect_sizes = effect_sizes[top_genes])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = NULL)
                 },
                 
                 nmf = {
                   if (!requireNamespace("NMF", quietly = TRUE)) {
                     stop("NMF package required. Install with: install.packages('NMF')")
                   }
                   
                   # [FIX-NMF] 'nmf' 블록 내부에서 보정 실행
                   if (!is.null(control_vars)) {
                     if (!requireNamespace("limma", quietly = TRUE)) {
                       stop("limma package required for NMF confounder correction")
                     }
                     warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s. Re-shifting to non-negative.",
                                     paste(control_vars, collapse=", ")))
                     
                     covariates_df <- meta.data[, control_vars, drop = FALSE]
                     if(any(sapply(covariates_df, is.factor)) || any(sapply(covariates_df, is.character))) {
                        covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
                     } else {
                        covariate_mat <- as.matrix(covariates_df)
                     }
                     
                     expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
                   }
                   
                   # [FIX-NMF] 보정 후 (혹은 보정 없이) 무조건 양수로 이동
                   expr_mat_pos <- expr_mat - min(expr_mat) + 0.01 
                   
                   rank <- min(n_groups + 2, 10)
                   dots <- list(...)
                   dots$seed <- NULL 
                   
                   nmf_res <- do.call(NMF::nmf, c(list(x = expr_mat_pos, rank = rank, seed=fgs_seed), dots))
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
                   
                   top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, length(weights_magnitude))])
                   weights_magnitude <- weights_magnitude[top_genes]
                   
                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]
                     g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     warning("nmf: n_groups > 2. Score represents magnitude (importance), not direction.")
                     weights <- weights_magnitude
                   }
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(component = best_component, correlation = component_cors[best_component])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) nmf_res else NULL)
                 },
                 
                 gam = {
                   if (!requireNamespace("mgcv", quietly = TRUE)) {
                     stop("mgcv package required. Install with: install.packages('mgcv')")
                   }
                   
                   deviance_explained <- numeric(nrow(expr_mat))
                   names(deviance_explained) <- rownames(expr_mat)
                   
                   if (is.null(control_vars)) {
                     formula_base <- ""
                     gam_data <- data.frame(y_var = target_binary)
                   } else {
                     # Ensure control vars are in the data frame
                     gam_data <- data.frame(y_var = target_binary, 
                                            meta.data[, control_vars, drop=FALSE])
                     # Create formula string *from names in gam_data*
                     formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse=" + "))
                     # Add the covariates matrix to gam_data
                     gam_data <- cbind(gam_data, covariate_mat_model)
                   }
                   
                   # [FIX-GAM] Use all genes remaining in expr_mat (1:nrow)
                   for (i in 1:nrow(expr_mat)) {
                     gam_data$gene_expr <- expr_mat[i, ]
                     
                     if (n_groups == 2) {
                       gam_data$y_var_numeric <- as.numeric(gam_data$y_var) - 1
                       full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
                       gam_fit <- try(mgcv::gam(as.formula(full_formula_str), 
                                            data = gam_data, family="binomial", ...), silent=TRUE)
                     } else {
                       gam_data$y_var_numeric <- as.numeric(gam_data$y_var)
                       full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
                       gam_fit <- try(mgcv::gam(as.formula(full_formula_str), 
                                            data = gam_data, ...), silent=TRUE)
                     }
                     
                     if(inherits(gam_fit, "try-error")) {
                        deviance_explained[i] <- 0
                     } else {
                        deviance_explained[i] <- summary(gam_fit)$dev.expl
                     }
                   }
                   
                   weights_magnitude <- deviance_explained
                   top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0, na.rm=TRUE))])
                   
                   if (length(top_genes) == 0) {
                     warning("GAM method found no genes with deviance > 0.")
                     return(list(genes=character(0), weights=numeric(0), scores=numeric(0), performance=list()))
                   }
                   
                   weights_magnitude <- weights_magnitude[top_genes]

                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]
                     g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     warning("gam: n_groups > 2. Score represents magnitude (importance), not direction.")
                     weights <- weights_magnitude
                   }
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(deviance_explained = weights_magnitude)
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = NULL)
                 },
                 
                 pca_loadings = {
                   # [V4.1] control_vars가 있었다면 3.8에서 expr_mat이 이미 보정됨.
                   # PCA는 스케일링된 데이터를 사용해야 함 (3.x에서 처리됨)
                   
                   pca_res <- prcomp(X, center=FALSE, scale.=FALSE) 
                   
                   n_pcs_to_test <- min(50, ncol(pca_res$x))
                   pc_cors <- numeric(n_pcs_to_test)
                   
                   for (k in 1:n_pcs_to_test) {
                     if (n_groups == 2) {
                       pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
                     } else {
                       aov_res <- try(summary(aov(pca_res$x[, k] ~ target_binary)), silent=TRUE)
                       if(inherits(aov_res, "try-error")) {
                           pc_cors[k] <- 0
                       } else {
                           pc_cors[k] <- aov_res[[1]][1, "F value"]
                       }
                     }
                   }
                   
                   best_pc <- which.max(pc_cors)
                   weights_all <- pca_res$rotation[, best_pc]
                   
                   top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
                   weights <- weights_all[top_genes]
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(PC = best_pc, correlation = pc_cors[best_pc],
                                variance_explained = summary(pca_res)$importance[2, best_pc])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) pca_res else NULL)
                 }
  )
  
  # ===
  # 5. Return results (v4와 동일)
  # ===
  
  result$method <- method
  result$target_var <- target_var
  result$n_groups <- n_groups
  result$n_cells <- ncol(expr_mat)
  result$formula <- paste(deparse(match.call()), collapse = " ")
  
  class(result) <- c("gene_signature", "list")
  return(result)
}

#' @export
find_gene_signature_v4.2 <- function(data, 
                                 meta.data = NULL,
                                 target_var,
                                 target_group = NULL,
                                 control_vars = NULL,   
                                 method = c("tree_based", "lasso", "limma", 
                                            "nmf", "wilcoxon", "gam", "pca_loadings"),
                                 n_features = 50,
                                 test_n = NULL,         
                                 preprocess = TRUE,
                                 min_cells = 10,
                                 min_pct = 0.01,
                                 return_model = FALSE,
                                 fgs_seed = 42,
                                 lambda_selection = "lambda.1se",
                                 ...) {
  
  # ===
  # 0. Batch Processing (v4와 동일)
  # ===
  
  all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
  if (is.null(method)) {
    method <- all_methods
  }
  
  if (length(method) > 1) {
    current_call <- match.call()
    
    results_list <- lapply(method, function(m) {
      single_call <- current_call
      single_call$method <- m
      tryCatch({
        eval(single_call)
      }, error = function(e) {
        warning(sprintf("Method '%s' failed with error: %s", m, e$message))
        return(list(method = m, error = e$message))
      })
    })
    
    names(results_list) <- method
    return(results_list)
  }
  
  if (!method %in% all_methods) {
    stop(sprintf("Invalid method '%s'. Choose from: %s", 
                 method, paste(all_methods, collapse=", ")))
  }
  
  set.seed(fgs_seed)
  
  # ===
  # 1. Input validation and data extraction (v4와 동일)
  # ===
  
  is_seurat <- inherits(data, "Seurat")
  
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required but not installed")
    }
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    # 기본 Assay 설정 
    default_assay <- Seurat::DefaultAssay(data)
    
    # 데이터 추출 우선순위:
    # 1. 'data' 슬롯 (일반적으로 정규화된 데이터)
    # 2. 'counts' 슬롯 (로우 데이터)
    
    # tryCatch를 사용하여 V5 레이어 객체와 V4/V3 slot 객체 모두 호환되도록 처리
    expr_mat <- tryCatch({
        # V5 방식 (레이어 접근) 시도
        as.matrix(Seurat::GetAssayData(data, assay = default_assay, layer = "data"))
    }, error = function(e) {
        # V4/V3 방식 (슬롯 접근) 시도 또는 'data' 슬롯이 없을 경우
        tryCatch({
            as.matrix(Seurat::GetAssayData(data, assay = default_assay, slot = "data"))
        }, error = function(e) {
            # 'counts' 슬롯으로 최종 시도
            as.matrix(Seurat::GetAssayData(data, assay = default_assay, slot = "counts"))
        })
    })
  } else {
    if (is.null(meta.data)) {
      stop("meta.data must be provided when data is not a Seurat object")
    }
    expr_mat <- as.matrix(data)
    if (nrow(expr_mat) == nrow(meta.data)) {
      expr_mat <- t(expr_mat)
    }
  }
  
  if (!target_var %in% colnames(meta.data)) {
    stop(sprintf("target_var '%s' not found in metadata columns: %s", 
                 target_var, paste(colnames(meta.data), collapse=", ")))
  }
  
  if (!is.null(control_vars)) {
    missing_vars <- control_vars[!control_vars %in% colnames(meta.data)]
    if (length(missing_vars) > 0) {
      stop(sprintf("control_vars not found in metadata: %s", 
                   paste(missing_vars, collapse=", ")))
    }
  }
  
  common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
  if (length(common_cells) == 0) {
    stop("No common cells found between expression matrix and metadata")
  }
  expr_mat <- expr_mat[, common_cells]
  meta.data <- meta.data[common_cells, ]
  
  # ===
  # 2. Process target variable (v4.1 수정)
  # ===
  
  target_values <- meta.data[[target_var]]
  target_type <- class(target_values)[1]
  
  if (is.numeric(target_values)) {
    if (is.null(target_group)) {
      target_group <- 0.25
    }
    if (is.list(target_group)) {
      low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
      high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
      group_labels <- ifelse(target_values <= low_cutoff, "Low",
                             ifelse(target_values >= high_cutoff, "High", NA))
    } else if (length(target_group) == 1 && target_group < 1) {
      low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
      high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
      group_labels <- ifelse(target_values <= low_cutoff, "Low",
                             ifelse(target_values >= high_cutoff, "High", NA))
    } else {
      group_labels <- ifelse(target_values < target_group, "Low", "High")
    }
    keep_cells <- !is.na(group_labels)
  } else {
    if (!is.null(target_group)) {
      keep_cells <- target_values %in% target_group
    } else {
      keep_cells <- rep(TRUE, length(target_values))
    }
    group_labels <- target_values
    keep_cells <- keep_cells & !is.na(group_labels)
  }
  
  expr_mat <- expr_mat[, keep_cells]
  meta.data <- meta.data[keep_cells, ]
  
  # [FIX-NA] 사용하지 않는 팩터 레벨 제거 (model.matrix 에러 방지)
  target_binary <- factor(group_labels[keep_cells])
  
  if (!is.null(control_vars)) {
    for (cv in control_vars) {
      if (is.factor(meta.data[[cv]])) {
        meta.data[[cv]] <- droplevels(meta.data[[cv]])
      }
    }
  }

  if (length(unique(target_binary)) < 2) {
    stop("Target variable must have at least 2 groups after processing")
  }
  n_groups <- length(unique(target_binary))
  
  # ===
  # 3. Filter and preprocess genes (v4와 동일)
  # ===
  
  n_cells_expr <- rowSums(expr_mat > 0)
  pct_cells_expr <- n_cells_expr / ncol(expr_mat)
  keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  expr_mat <- expr_mat[keep_genes, ]
  
  if (nrow(expr_mat) == 0) {
    stop("No genes pass filtering criteria")
  }
  
  # ===
  # 3.5 [V4] Speed-up: Pre-filter genes using test_n
  # ===
  if (!is.null(test_n) && nrow(expr_mat) > test_n) {
    message(sprintf("Pre-filtering to top %d genes based on limma p-value...", test_n))
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for 'test_n' pre-filtering")
    }
    
    if (is.null(control_vars)) {
      design_test <- model.matrix(~ target_binary, data = meta.data)
    } else {
      formula_test <- as.formula(paste("~ target_binary +", paste(control_vars, collapse="+")))
      design_test <- model.matrix(formula_test, data = meta.data)
    }
    
    fit_test <- limma::lmFit(expr_mat, design_test)
    fit_test <- limma::eBayes(fit_test)
    
    coef_indices <- grep("target_binary", colnames(design_test))
    if (length(coef_indices) == 0) {
        warning("Could not find target_binary coefficients for test_n pre-filtering.")
        # Fallback: use all coefficients (less ideal)
        coef_indices <- 1:ncol(design_test)
    }
    
    top_table_test <- limma::topTable(fit_test, 
                                      coef = coef_indices, 
                                      number = Inf, sort.by = "P")
    
    top_gene_names <- rownames(top_table_test)[1:min(test_n, nrow(top_table_test))]
    expr_mat <- expr_mat[top_gene_names, ]
    message(sprintf("... reduced to %d genes.", nrow(expr_mat)))
  }
  
  if (preprocess) {
    if (max(expr_mat) > 100) {
      expr_mat <- log1p(expr_mat)
    }
    if (method %in% c("lasso", "gam", "pca_loadings")) {
      gene_means <- rowMeans(expr_mat)
      gene_sds <- apply(expr_mat, 1, sd)
      gene_sds[gene_sds == 0] <- 1
      expr_mat <- (expr_mat - gene_means) / gene_sds
    }
  }
  
  # ===
  # 3.8 [V4.1] Confounder pre-correction (non-native methods)
  # ===
  
  # [FIX-NMF] 'nmf'는 removeBatchEffect 대상에서 제외 (음수 값 방지)
  if (!is.null(control_vars) && method %in% c("wilcoxon", "pca_loadings")) { 
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for removeBatchEffect")
    }
    
    message(sprintf("Method %s: Applying limma::removeBatchEffect for: %s",
                    method, paste(control_vars, collapse=", ")))
    
    covariates_df <- meta.data[, control_vars, drop = FALSE]
    
    # Check for factors and create model matrix
    if(any(sapply(covariates_df, is.factor)) || any(sapply(covariates_df, is.character))) {
        covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
    } else {
        covariate_mat <- as.matrix(covariates_df)
    }
    
    expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
  }
  
  
  # ===
  # 4. Method-specific signature discovery (v4.1 수정)
  # ===
  
  X <- t(expr_mat)
  y <- target_binary
  
  covariate_mat_model <- NULL # Initialize
  if (!is.null(control_vars) && method %in% c("tree_based", "lasso", "gam")) {
     covariates_df_model <- meta.data[, control_vars, drop = FALSE]
     if(any(sapply(covariates_df_model, is.factor)) || any(sapply(covariates_df_model, is.character))) {
        covariate_mat_model <- model.matrix(~ . - 1, data = covariates_df_model)
     } else {
        covariate_mat_model <- as.matrix(covariates_df_model)
     }
  }

  
  result <- switch(method,
                 
                 tree_based = {
                   if (!requireNamespace("randomForest", quietly = TRUE)) {
                     stop("randomForest package required. Install with: install.packages('randomForest')")
                   }
                   
                   X_rf <- X
                   
                   if (!is.null(control_vars)) {
                     X_rf <- cbind(X_rf, covariate_mat_model)
                   }
                   
                   # Pre-filter genes if matrix is still too large
                   if (ncol(X_rf) > 2000) {
                     gene_vars <- apply(X, 2, var) # var from original X
                     top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
                     
                     control_var_names <- if(is.null(control_vars)) character(0) else colnames(covariate_mat_model)
                     X_rf <- X_rf[, c(top_var_genes, control_var_names)]
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
                   
                   top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
                   weights_magnitude <- weights_magnitude_genes[top_genes]
                   
                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]
                     g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1 
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     warning("tree_based: n_groups > 2. Score represents magnitude (importance), not direction.")
                     weights <- weights_magnitude
                   }
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   pred <- rf_model$predicted
                   if (n_groups == 2) {
                     if (requireNamespace("pROC", quietly = TRUE)) {
                       roc_obj <- pROC::roc(y, scores, quiet=TRUE)
                       auc <- as.numeric(pROC::auc(roc_obj))
                     } else { auc <- NA }
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
                   } else {
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, confusion = table(pred, y))
                   }
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) rf_model else NULL)
                 },
                 
                 lasso = {
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
                   
                   if (n_groups == 2) {
                     cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=1, 
                                                 penalty.factor = penalty_vec, ...)
                   } else {
                     cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=1, 
                                                 penalty.factor = penalty_vec, ...)
                   }
                   
                   coefs <- coef(cv_fit, s = lambda_selection) 
                   
                   if (n_groups == 2) {
                     weights_all <- as.numeric(coefs[2:(ncol(X)+1)]) 
                     names(weights_all) <- rownames(coefs)[2:(ncol(X)+1)]
                   } else {
                     coef_list <- lapply(coefs, function(x) as.numeric(x[2:(ncol(X)+1)]))
                     weights_all <- rowMeans(do.call(cbind, coef_list))
                     names(weights_all) <- rownames(coefs[[1]])[2:(ncol(X)+1)]
                   }
                   
                   if (length(weights_all) == 0) {
                     stop("LASSO returned no gene coefficients.")
                   }

                   # [FIX-Lasso] 'nonzero' 필터링 제거. 항상 abs(weight)로 정렬.
                   top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
                   
                   weights <- weights_all[top_genes]
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   pred_probs <- predict(cv_fit, newx=X_model, s = lambda_selection, type="response")
                   if (n_groups == 2) {
                     pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
                     if (requireNamespace("pROC", quietly = TRUE)) {
                       roc_obj <- pROC::roc(y, scores, quiet=TRUE) 
                       auc <- as.numeric(pROC::auc(roc_obj))
                     } else { auc <- NA }
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
                   } else {
                     pred_class_indices <- apply(pred_probs[,,1], 1, which.max)
                     pred <- levels(y)[pred_class_indices]
                     acc <- mean(pred == y)
                     perf <- list(accuracy = acc, confusion = table(pred, y))
                   }
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) cv_fit else NULL)
                 },
                 
                 limma = {
                   if (!requireNamespace("limma", quietly = TRUE)) {
                     stop("limma package required. Install with: BiocManager::install('limma')")
                   }
                   
                   if (is.null(control_vars)) {
                     design <- model.matrix(~0 + target_binary, data = meta.data)
                     colnames(design)[1:n_groups] <- levels(target_binary)
                   } else {
                     control_formula <- paste(control_vars, collapse = " + ")
                     full_formula <- as.formula(paste("~0 + target_binary +", control_formula))
                     design <- model.matrix(full_formula, data = meta.data)
                     colnames(design)[1:n_groups] <- levels(target_binary) 
                   }

                   fit <- limma::lmFit(expr_mat, design)
                   
                   if (n_groups == 2) {
                     contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
                   } else {
                     # Create all pairwise contrasts if n_groups > 2
                     contrast_pairs <- combn(levels(target_binary), 2, function(x) paste(x[2], x[1], sep="-"))
                     contrast_str <- contrast_pairs
                     warning("limma: n_groups > 2. Using all pairwise contrasts. Weights based on average t-statistic.")
                   }
                   
                   contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
                   fit2 <- limma::contrasts.fit(fit, contrast_mat)
                   fit2 <- limma::eBayes(fit2)
                   
                   if(n_groups == 2) {
                       top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
                       weights_all <- top_table$t
                   } else {
                       # For multiple contrasts, sort by F-statistic
                       top_table <- limma::topTable(fit2, number=Inf, sort.by="F")
                       # Use average absolute t-stat as weight magnitude
                       weights_all <- rowMeans(abs(fit2$t[, contrast_str, drop=FALSE]))
                   }
                   
                   names(weights_all) <- rownames(top_table)
                   
                   top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
                   
                   # Re-fetch signed weights for top genes
                   if(n_groups == 2) {
                       weights <- top_table[top_genes, "t"]
                   } else {
                       # For n_groups > 2, weights are just magnitude (no clear direction)
                       weights <- weights_all[top_genes]
                   }

                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(top_table = top_table[top_genes, ])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) fit2 else NULL)
                 },
                 
                 wilcoxon = {
                   # [V4.1] control_vars가 있었다면 3.8에서 expr_mat이 이미 보정됨.
                   
                   pvals <- numeric(nrow(expr_mat))
                   effect_sizes <- numeric(nrow(expr_mat))
                   names(pvals) <- rownames(expr_mat)
                   names(effect_sizes) <- rownames(expr_mat)

                   for (i in 1:nrow(expr_mat)) {
                     if (n_groups == 2) {
                       group1 <- expr_mat[i, target_binary == levels(target_binary)[1]]
                       group2 <- expr_mat[i, target_binary == levels(target_binary)[2]]
                       
                       # [FIX-Wilcox] '...' 인자 제거
                       test <- try(wilcox.test(group1, group2), silent=TRUE) 
                       
                       if(inherits(test, "try-error")) {
                         pvals[i] <- 1.0
                         effect_sizes[i] <- 0
                       } else {
                         pvals[i] <- test$p.value
                         effect_sizes[i] <- median(group2) - median(group1)
                       }
                     } else {
                       # [FIX-Wilcox] '...' 인자 제거
                       test <- try(kruskal.test(expr_mat[i, ] ~ target_binary), silent=TRUE) 
                       
                       if(inherits(test, "try-error")) {
                          pvals[i] <- 1.0
                          effect_sizes[i] <- 0
                       } else {
                          pvals[i] <- test$p.value
                          # Use variance of medians as effect size (magnitude only)
                          effect_sizes[i] <- var(tapply(expr_mat[i, ], target_binary, median))
                       }
                     }
                   }
                   
                   padj <- p.adjust(pvals, method="BH")
                   
                   ranks <- rank(pvals) + rank(-abs(effect_sizes))
                   top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
                   
                   weights <- effect_sizes[top_genes]
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(pvalues = pvals[top_genes], padj = padj[top_genes],
                                effect_sizes = effect_sizes[top_genes])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = NULL)
                 },
                 
                 nmf = {
                   if (!requireNamespace("NMF", quietly = TRUE)) {
                     stop("NMF package required. Install with: install.packages('NMF')")
                   }
                   
                   # [FIX-NMF] 'nmf' 블록 내부에서 보정 실행
                   if (!is.null(control_vars)) {
                     if (!requireNamespace("limma", quietly = TRUE)) {
                       stop("limma package required for NMF confounder correction")
                     }
                     warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s. Re-shifting to non-negative.",
                                     paste(control_vars, collapse=", ")))
                     
                     covariates_df <- meta.data[, control_vars, drop = FALSE]
                     if(any(sapply(covariates_df, is.factor)) || any(sapply(covariates_df, is.character))) {
                        covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
                     } else {
                        covariate_mat <- as.matrix(covariates_df)
                     }
                     
                     expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
                   }
                   
                   # [FIX-NMF] 보정 후 (혹은 보정 없이) 무조건 양수로 이동
                   expr_mat_pos <- expr_mat - min(expr_mat) + 0.01 
                   
                   rank <- min(n_groups + 2, 10)
                   dots <- list(...)
                   dots$seed <- NULL 
                   # [V4.2 FIX] 'method' 인자 충돌 방지
                   dots$method <- NULL

                   nmf_res <- do.call(NMF::nmf, c(list(x = expr_mat_pos, rank = rank, seed=fgs_seed), dots))
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
                   
                   top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, length(weights_magnitude))])
                   weights_magnitude <- weights_magnitude[top_genes]
                   
                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]
                     g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     warning("nmf: n_groups > 2. Score represents magnitude (importance), not direction.")
                     weights <- weights_magnitude
                   }
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(component = best_component, correlation = component_cors[best_component])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) nmf_res else NULL)
                 },
                 
                 gam = {
                   if (!requireNamespace("mgcv", quietly = TRUE)) {
                     stop("mgcv package required. Install with: install.packages('mgcv')")
                   }
                   
                   deviance_explained <- numeric(nrow(expr_mat))
                   names(deviance_explained) <- rownames(expr_mat)
                   
                   if (is.null(control_vars)) {
                     formula_base <- ""
                     gam_data <- data.frame(y_var = target_binary)
                   } else {
                     # Ensure control vars are in the data frame
                     gam_data <- data.frame(y_var = target_binary, 
                                            meta.data[, control_vars, drop=FALSE])
                     # Create formula string *from names in gam_data*
                     formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse=" + "))
                     # Add the covariates matrix to gam_data
                     gam_data <- cbind(gam_data, covariate_mat_model)
                   }
                   
                   # [FIX-GAM] Use all genes remaining in expr_mat (1:nrow)
                   for (i in 1:nrow(expr_mat)) {
                     gam_data$gene_expr <- expr_mat[i, ]
                     
                     if (n_groups == 2) {
                       gam_data$y_var_numeric <- as.numeric(gam_data$y_var) - 1
                       full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
                       gam_fit <- try(mgcv::gam(as.formula(full_formula_str), 
                                            data = gam_data, family="binomial", ...), silent=TRUE)
                     } else {
                       gam_data$y_var_numeric <- as.numeric(gam_data$y_var)
                       full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
                       gam_fit <- try(mgcv::gam(as.formula(full_formula_str), 
                                            data = gam_data, ...), silent=TRUE)
                     }
                     
                     if(inherits(gam_fit, "try-error")) {
                        deviance_explained[i] <- 0
                     } else {
                        deviance_explained[i] <- summary(gam_fit)$dev.expl
                     }
                   }
                   
                   weights_magnitude <- deviance_explained
                   top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0, na.rm=TRUE))])
                   
                   if (length(top_genes) == 0) {
                     warning("GAM method found no genes with deviance > 0.")
                     return(list(genes=character(0), weights=numeric(0), scores=numeric(0), performance=list()))
                   }
                   
                   weights_magnitude <- weights_magnitude[top_genes]

                   if (n_groups == 2) {
                     g1_cells <- y == levels(y)[1]
                     g2_cells <- y == levels(y)[2]
                     mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
                     mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
                     effect_size <- mean_g2 - mean_g1
                     weights <- weights_magnitude * sign(effect_size)
                   } else {
                     warning("gam: n_groups > 2. Score represents magnitude (importance), not direction.")
                     weights <- weights_magnitude
                   }
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(deviance_explained = weights_magnitude)
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = NULL)
                 },
                 
                 pca_loadings = {
                   # [V4.1] control_vars가 있었다면 3.8에서 expr_mat이 이미 보정됨.
                   # PCA는 스케일링된 데이터를 사용해야 함 (3.x에서 처리됨)
                   
                   pca_res <- prcomp(X, center=FALSE, scale.=FALSE) 
                   
                   n_pcs_to_test <- min(50, ncol(pca_res$x))
                   pc_cors <- numeric(n_pcs_to_test)
                   
                   for (k in 1:n_pcs_to_test) {
                     if (n_groups == 2) {
                       pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
                     } else {
                       aov_res <- try(summary(aov(pca_res$x[, k] ~ target_binary)), silent=TRUE)
                       if(inherits(aov_res, "try-error")) {
                           pc_cors[k] <- 0
                       } else {
                           pc_cors[k] <- aov_res[[1]][1, "F value"]
                       }
                     }
                   }
                   
                   best_pc <- which.max(pc_cors)
                   weights_all <- pca_res$rotation[, best_pc]
                   
                   top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
                   weights <- weights_all[top_genes]
                   
                   scores <- as.numeric(X[, top_genes] %*% weights)
                   names(scores) <- rownames(X)
                   
                   perf <- list(PC = best_pc, correlation = pc_cors[best_pc],
                                variance_explained = summary(pca_res)$importance[2, best_pc])
                   
                   list(genes = top_genes, weights = weights, scores = scores,
                        performance = perf, model = if(return_model) pca_res else NULL)
                 }
  )
  
  # ===
  # 5. Return results (v4와 동일)
  # ===
  
  result$method <- method
  result$target_var <- target_var
  result$n_groups <- n_groups
  result$n_cells <- ncol(expr_mat)
  result$formula <- paste(deparse(match.call()), collapse = " ")
  
  class(result) <- c("gene_signature", "list")
  return(result)
}



#' FGS v5 Preprocessing Helper Function
#'
#' @description
#' find_gene_signature_v5의 전처리를 담당하는 헬퍼 함수.
#' 1. NA 제거 (target_var, control_vars 기준)
#' 2. expr_mat 추출 (as.matrix)
#' 3. test_n 필터링 (limma, 단 한 번 실행)
#' 4. Preprocessing (log1p, scale)
#' 5. Confounder pre-correction (removeBatchEffect, 단 한 번 실행)
#'
#' @return 전처리가 완료된 데이터 객체(list)
#'
fgs_preprocess_data <- function(data, 
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
  
  # === 1. Input validation (v4.2와 유사) ===
  is_seurat <- inherits(data, "Seurat")
  
  if (is_seurat) {
    if (is.null(meta.data)) {
      meta.data <- data@meta.data
    }
    # V4.2의 안정적인 데이터 추출 로직
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
  
  # === 2. [V5 핵심] NA 및 팩터 레벨 완벽 제거 ===
  message("V5 Preprocessing: 1. Cleaning NA values...")
  vars_to_check <- c(target_var, control_vars)
  
  # control_vars가 NULL이 아닐 때만 vars_to_check에 포함
  if (is.null(control_vars)) {
    vars_to_check <- target_var
  }
  
  complete_cases_idx <- complete.cases(meta.data[, vars_to_check, drop=FALSE])
  
  n_removed <- sum(!complete_cases_idx)
  if (n_removed > 0) {
    warning(sprintf("V5 Preprocessing: Removed %d cells with NAs in target/control vars.", n_removed))
  }
  
  meta.data <- meta.data[complete_cases_idx, ]
  expr_mat <- expr_mat[, rownames(meta.data)]

  # 팩터 레벨 정리
  target_binary <- factor(meta.data[[target_var]])
  meta.data$target_binary_var <- target_binary # 포뮬러용 변수
  
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
  
  # === 4. [V5 핵심] test_n 필터링 (단 한 번 실행) ===
  if (!is.null(test_n) && nrow(expr_mat) > test_n) {
    message(sprintf("V5 Preprocessing: 3. Pre-filtering to top %d genes (limma)...", test_n))
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for 'test_n' pre-filtering")
    }
    
    # NA가 제거된 meta.data로 디자인 생성
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
  
  # 스케일링이 필요한 메서드를 위한 데이터 (원본 expr_mat 보존)
  expr_mat_scaled <- NULL
  if (length(methods_requiring_scale) > 0) {
    gene_means <- rowMeans(expr_mat)
    gene_sds <- apply(expr_mat, 1, sd)
    gene_sds[gene_sds == 0] <- 1
    expr_mat_scaled <- (expr_mat - gene_means) / gene_sds
  }

  # === 6. [V5 핵심] Confounder pre-correction (단 한 번 실행) ===
  message("V5 Preprocessing: 5. Applying confounder correction (if needed)...")
  
  expr_mat_corrected <- NULL
  if (!is.null(control_vars) && length(methods_requiring_correction) > 0) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package required for removeBatchEffect")
    }
    
    covariates_df <- meta.data[, control_vars, drop = FALSE]
    covariate_mat <- model.matrix(~ . - 1, data = covariates_df) 
    
    # 스케일링 *전*의 expr_mat을 보정 (pca_loadings는 스케일링된 매트릭스를 보정해야 함)
    expr_mat_corrected <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
  }

  # --- 반환 객체 ---
  list(
    expr_mat = expr_mat,                         # 원본 (log1p만 적용)
    expr_mat_scaled = expr_mat_scaled,           # 스케일링된 버전
    expr_mat_corrected = expr_mat_corrected,     # 교란변수 보정된 버전 (스케일링 X)
    meta.data = meta.data,                       # NA 제거된 메타데이터
    target_binary = target_binary,               # 타겟 변수
    n_groups = length(unique(target_binary)),
    control_vars = control_vars,
    fgs_seed = 42 # 임시 (v5에서는 fgs_seed를 메인에서 받도록)
  )
}

# [V5] find_gene_signature_v5
#
# 변경점:
# 1. 'lapply', 'match.call' 제거.
# 2. 'fgs_preprocess_data' 헬퍼 함수를 처음에 1회 호출.
# 3. 'for (m in method)' 루프로 메서드를 순회. (모든 반복 문제 해결)
# 4. 'nmf' 블록: 'do.call'의 'method' 인자 충돌을 'dots$method <- NULL'로 해결.
# 5. 'gam' 블록: 'tryCatch'를 사용하여 수렴 경고를 캡처하고 1회만 요약 리포트.

#' @export
find_gene_signature_v5 <- function(data, 
                                 meta.data = NULL,
                                 target_var,
                                 target_group = NULL,
                                 control_vars = NULL,   
                                 method = c("tree_based", "lasso", "limma", 
                                            "nmf", "wilcoxon", "gam", "pca_loadings"),
                                 n_features = 50,
                                 test_n = NULL,         
                                 preprocess = TRUE,
                                 min_cells = 10,
                                 min_pct = 0.01,
                                 return_model = FALSE,
                                 fgs_seed = 42,
                                 lambda_selection = "lambda.1se",
                                 ...) {
  
  all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
  if (is.null(method)) {
    method <- all_methods
  }
  
  # === 1. [V5] 전처리 (단 1회 실행) ===
  
  # 각 메서드가 어떤 데이터가 필요한지 정의
  methods_requiring_scale <- c("lasso", "gam", "pca_loadings")
  methods_requiring_correction <- c("wilcoxon", "pca_loadings") # nmf는 내부 처리
  
  # 헬퍼 함수 호출
  preprocessed_data <- fgs_preprocess_data(
    data = data, 
    meta.data = meta.data, 
    target_var = target_var, 
    target_group = target_group, 
    control_vars = control_vars, 
    test_n = test_n, 
    preprocess = preprocess,
    min_cells = min_cells,
    min_pct = min_pct,
    methods_requiring_scale = intersect(method, methods_requiring_scale),
    methods_requiring_correction = intersect(method, methods_requiring_correction)
  )

  # 전처리된 객체 할당
  expr_mat_base <- preprocessed_data$expr_mat
  meta.data_clean <- preprocessed_data$meta.data
  target_binary <- preprocessed_data$target_binary
  n_groups <- preprocessed_data$n_groups
  
  set.seed(preprocessed_data$fgs_seed) # 헬퍼에서 정의된 시드 사용
  
  # 교란변수 model.matrix (native 지원 메서드용)
  covariate_mat_model <- NULL
  if (!is.null(control_vars)) {
     covariates_df_model <- meta.data_clean[, control_vars, drop = FALSE]
     covariate_mat_model <- model.matrix(~ . - 1, data = covariates_df_model)
  }

  # --- [V5] 결과 저장을 위한 리스트 ---
  results_list <- list()
  
  
  # === 2. [V5] 메서드 순회 (lapply 대신 for 루프) ===
  
  for (m in method) {
    
    if (!m %in% all_methods) {
      warning(sprintf("Invalid method '%s'. Skipping.", m))
      next
    }
    
    message(sprintf("--- Running Method: %s ---", m))
    
    tryCatch({
      
      # === 3. 데이터 선택 (메서드별) ===
      
      # 1) expr_mat (유전자 x 샘플) 선택
      if (m %in% c("limma", "wilcoxon", "nmf")) {
        expr_mat_method <- expr_mat_base
      } else if (m %in% c("lasso", "gam", "pca_loadings")) {
        expr_mat_method <- preprocessed_data$expr_mat_scaled
      } else { # tree_based
        expr_mat_method <- expr_mat_base 
      }
      
      # 2) 교란변수 사전 보정된 데이터 사용
      if (m %in% c("wilcoxon", "pca_loadings")) {
        # 'pca_loadings'는 스케일링된 + 보정된 데이터가 필요 (v5.1 개선 필요)
        # 우선순위: 보정된 데이터가 있으면 사용
        if (!is.null(preprocessed_data$expr_mat_corrected)) {
          expr_mat_method <- preprocessed_data$expr_mat_corrected
        }
      }

      # 3) X (샘플 x 유전자) 생성
      X <- t(expr_mat_method)
      y <- target_binary

      
      # === 4. Method-specific (v4.1 로직과 대부분 동일) ===
      
      result <- switch(m,
        
        tree_based = {
          # ... (v4.1 tree_based 로직과 동일) ...
          # 단, X_rf 생성 시 expr_mat_method가 아닌 'X' (t(expr_mat_base)) 사용
          X_rf <- t(expr_mat_base) 
          if (!is.null(control_vars)) { X_rf <- cbind(X_rf, covariate_mat_model) }
          # ...
          list(...) # v4.1 코드 참조
        },
        
        lasso = {
          # ... (v4.1 lasso 로직과 동일) ...
          # X는 이미 스케일링된 t(expr_mat_scaled)임
          if (is.null(control_vars)) {
            X_model <- X
            penalty_vec <- rep(1, ncol(X_model))
          } else {
            X_model <- cbind(X, covariate_mat_model)
            penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
          }
          # ...
          list(...) # v4.1 코드 참조
        },
        
        limma = {
          # ... (v4.1 limma 로직과 동일) ...
          # expr_mat_method는 스케일링 안 된 'expr_mat_base'임
          if (is.null(control_vars)) {
            design <- model.matrix(~0 + target_binary_var, data = meta.data_clean)
          } else {
            control_formula <- paste(control_vars, collapse = " + ")
            full_formula <- as.formula(paste("~0 + target_binary_var +", control_formula))
            design <- model.matrix(full_formula, data = meta.data_clean)
          }
          colnames(design)[1:n_groups] <- levels(target_binary)
          # ...
          list(...) # v4.1 코드 참조
        },
        
        wilcoxon = {
          # ... (v4.1 wilcoxon 로직과 동일) ...
          # expr_mat_method는 이미 'expr_mat_corrected'임
          pvals <- numeric(nrow(expr_mat_method))
          # ...
          list(...) # v4.1 코드 참조
        },
        
        nmf = {
          # [V5 FIX] v4.1의 'nmf' 블록 로직 (내부 보정 + 인자 충돌 해결)
          if (!requireNamespace("NMF", quietly = TRUE)) { ... }
          
          # 'expr_mat_method'는 'expr_mat_base'임
          expr_mat_nmf <- expr_mat_method 
          
          if (!is.null(control_vars)) {
            warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s...", m))
            covariate_mat <- model.matrix(~ . - 1, data = meta.data_clean[, control_vars, drop=FALSE])
            expr_mat_nmf <- limma::removeBatchEffect(expr_mat_nmf, covariates = covariate_mat)
          }
          
          expr_mat_pos <- expr_mat_nmf - min(expr_mat_nmf) + 0.01 
          
          rank <- min(n_groups + 2, 10)
          dots <- list(...)
          
          # [V5 FIX] 'method' 및 'seed' 인자 충돌 방지
          dots$method <- NULL 
          dots$seed <- NULL 
          
          nmf_args <- c(list(x = expr_mat_pos, rank = rank, seed = fgs_seed), dots)
          nmf_res <- do.call(NMF::nmf, nmf_args)
          
          # ... (v4.1 nmf 로직) ...
          list(...) # v4.1 코드 참조
        },

        gam = {
          # [V5 FIX] v4.1의 'gam' 블록 + 경고 억제
          if (!requireNamespace("mgcv", quietly = TRUE)) { ... }
          
          deviance_explained <- numeric(nrow(expr_mat_method))
          names(deviance_explained) <- rownames(expr_mat_method)
          
          # GAM 데이터 준비 (X는 이미 스케일링됨)
          gam_data <- data.frame(y_var_numeric = as.numeric(target_binary) - 1)
          if (!is.null(control_vars)) {
            gam_data <- cbind(gam_data, covariate_mat_model)
            formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse=" + "))
          } else {
            formula_base <- ""
          }
          full_formula_str_base <- paste("y_var_numeric ~ s(gene_expr)", formula_base)

          # [V5 FIX] 경고 카운터
          convergence_warnings <- 0
          
          for (i in 1:nrow(expr_mat_method)) {
            gam_data$gene_expr <- expr_mat_method[i, ]
            
            # tryCatch로 경고 캡처
            fit_result <- tryCatch({
              if (n_groups == 2) {
                mgcv::gam(as.formula(full_formula_str_base), data = gam_data, family="binomial", ...)
              } else {
                mgcv::gam(as.formula(full_formula_str_base), data = gam_data, ...)
              }
            }, warning = function(w) {
              if (grepl("did not converge", w$message)) {
                convergence_warnings <<- convergence_warnings + 1
              }
              # 경고를 억제하고 결과 반환 (suppressWarning과 유사)
              invokeRestart("muffleWarning")
            }, error = function(e) {
              return(NULL) # 에러 발생 시 NULL 반환
            })
            
            if (is.null(fit_result) || inherits(fit_result, "try-error")) {
              deviance_explained[i] <- 0
            } else {
              deviance_explained[i] <- summary(fit_result)$dev.expl
            }
          }
          
          if (convergence_warnings > 0) {
            warning(sprintf("GAM: %d genes failed to converge (e.g., sparse data or NA issues).", convergence_warnings))
          }
          
          # ... (v4.1 gam 나머지 로직) ...
          list(...) # v4.1 코드 참조
        },
        
        pca_loadings = {
          # ... (v4.1 pca_loadings 로직) ...
          # X는 'expr_mat_corrected'의 전치행렬
          list(...) # v4.1 코드 참조
        }
      ) # --- switch 끝 ---

      # === 5. 결과 저장 ===
      result$method <- m
      result$target_var <- target_var
      result$n_groups <- n_groups
      result$n_cells <- ncol(expr_mat_base)
      
      class(result) <- c("gene_signature", "list")
      
      results_list[[m]] <- result
      
    }, error = function(e) {
      # --- 메서드 개별 에러 처리 ---
      warning(sprintf("Method '%s' failed with error: %s", m, e$message))
      results_list[[m]] <- list(method = m, error = e$message)
    })
    
  } # --- for 루프 끝 ---

  # === 6. [V5] 최종 반환 ===
  return(results_list)
}

#' FGS v5.2 Preprocessing Helper Function
#'
#' @description
#' find_gene_signature_v5.2의 전처리를 담당하는 헬퍼 함수.
#' 1. NA 제거 (target_var, control_vars 기준)
#' 2. expr_mat 추출 (as.matrix)
#' 3. test_n 필터링 (limma, 단 한 번 실행)
#' 4. Preprocessing (log1p, scale)
#' 5. Confounder pre-correction (removeBatchEffect, 단 한 번 실행)
#'
#' @export
fgs_preprocess_data_v5.2 <- function(data, 
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
  
  meta.data <- meta.data[complete_cases_idx, ]
  expr_mat <- expr_mat[, rownames(meta.data)]

  target_binary <- factor(meta.data[[target_var]])
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

# [V5.2] find_gene_signature_v5.2
#
# 변경점 (v5.1 대비):
# 1. [Req 2] `pca.n_pcs = 1` 인자 추가. 상위 N개 PC의 로딩 합으로 중요도 계산.
# 2. [Req 3] `tree_based` -> `random_forest`로 메서드 이름 변경.
# 3. [Req 4] 신규 베타 메서드 추가:
#    - `random_forest_ranger` (ranger)
#    - `ridge` (glmnet, alpha=0)
#    - `elastic_net` (glmnet, alpha=0.5 (기본값)), `enet.alpha` 인자 추가.
#    - `xgboost` (xgboost)
# 4. [Req 5] `all_methods` 순서 재정렬 (Tree, Regularization, Loadings, Stats)
# 5. [Req 6] `nmf` 블록: `...` 인자 필터링 로직('valid_nmf_args') 적용 (충돌 해결)
# 6. [Req 1] `gam` 블록: `gam.min_unique = 15`, `gam.k = 10` 인자 추가.
#    - 고유 값 15개 미만 유전자 스킵.
#    - `bam()` 사용하여 속도 향상.

#' @export 
find_gene_signature_v5.2 <- function(data, 
                                 meta.data = NULL,
                                 target_var,
                                 target_group = NULL,
                                 control_vars = NULL,   
                                 method = c("random_forest", "random_forest_ranger",
                                            "lasso", "ridge", "elastic_net",
                                            "pca_loadings", "nmf_loadings",
                                            "gam", "limma", "wilcoxon",
                                            "xgboost"),
                                 n_features = 50,
                                 test_n = NULL,         
                                 preprocess = TRUE,
                                 min_cells = 10,
                                 min_pct = 0.01,
#' @param fgs_seed The seed for reproducibility.
#' @param lambda_selection The lambda selection method.
#' @param enet.alpha The alpha for elastic net.
#' @param pca.n_pcs The number of principal components to use.
#' @param gam.min_unique The minimum number of unique values to keep a gene.
#' @param gam.k The number of knots to use.
#' @param ... Additional arguments to pass to the methods.
                                 return_model = FALSE,
                                 fgs_seed = 42,
                                 # --- 신규/수정된 인자 ---
                                 lambda_selection = "lambda.1se",
                                 enet.alpha = 0.5,        # (Elastic Net용)
                                 pca.n_pcs = 1,           # (PCA용)
                                 gam.min_unique = 15,     # (GAM용)
                                 gam.k = 10,              # (GAM용)
                                 ...) {
  
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
  
  # === 1. [V5] 전처리 (단 1회 실행) ===
  
  # [Req 4] 신규 모델 스케일링 요구사항 업데이트
  methods_requiring_scale <- c("lasso", "ridge", "elastic_net", 
                               "gam", "pca_loadings", "xgboost")
  methods_requiring_correction <- c("wilcoxon", "pca_loadings")
  
  preprocessed_data <- fgs_preprocess_data_v5.2(
    data = data, meta.data = meta.data, target_var = target_var, 
    target_group = target_group, control_vars = control_vars, 
    test_n = test_n, preprocess = preprocess, min_cells = min_cells,
    min_pct = min_pct,
    methods_requiring_scale = intersect(method, methods_requiring_scale),
    methods_requiring_correction = intersect(method, methods_requiring_correction)
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
  
  
  # === 2. [V5] 메서드 순회 (for 루프) ===
  
  for (m in method) {
    
    if (!m %in% all_methods) {
      warning(sprintf("Invalid method '%s'. Skipping.", m))
      next
    }
    
    message(sprintf("--- Running Method: %s ---", m))
    
    tryCatch({
      
      # === 3. 데이터 선택 (메서드별) ===
      
      if (m %in% c("limma", "wilcoxon", "nmf_loadings", "random_forest", "random_forest_ranger")) {
        expr_mat_method <- expr_mat_base
      } else if (m %in% c("lasso", "ridge", "elastic_net", "gam", "pca_loadings", "xgboost")) {
        # 'expr_mat_scaled'가 NULL이면 (필요한 메서드가 없었으면) 원본 사용
        expr_mat_method <- if(is.null(preprocessed_data$expr_mat_scaled)) expr_mat_base else preprocessed_data$expr_mat_scaled
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
          if (!is.null(control_vars)) { X_rf <- cbind(X_rf, covariate_mat_model) }
          
          rf_model <- randomForest::randomForest(x = X_rf, y = y, ntree = 500, importance = TRUE, ...)
          
          importance_scores <- randomForest::importance(rf_model)
          if (n_groups == 2) {
            weights_magnitude_all <- importance_scores[, "MeanDecreaseGini"]
          } else {
            weights_magnitude_all <- rowMeans(importance_scores[, grep("MeanDecreaseGini", colnames(importance_scores))])
          }
          
          gene_names_in_model <- colnames(X_rf)[!colnames(X_rf) %in% colnames(covariate_mat_model)]
          weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]
          
          top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
          weights_magnitude <- weights_magnitude_genes[top_genes]
          
          # 방향성 보정 및 performance 계산
          if (n_groups == 2) {
            g1_cells <- y == levels(y)[1]
            g2_cells <- y == levels(y)[2]
            mean_g1 <- colMeans(X_rf[g1_cells, top_genes, drop=FALSE])
            mean_g2 <- colMeans(X_rf[g2_cells, top_genes, drop=FALSE])
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
              roc_obj <- pROC::roc(y, scores, quiet=TRUE)
              auc <- as.numeric(pROC::auc(roc_obj))
            } else { auc <- NA }
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
          } else {
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, confusion = table(pred, y))
          }
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) rf_model else NULL)
        },
        
        random_forest_ranger = { # [Req 4] 신규
          if (!requireNamespace("ranger", quietly = TRUE)) {
            stop("ranger package required.")
          }
          
          # Ranger는 스케일링되지 않은 원본 X 사용
          X_ranger <- t(expr_mat_base) 
          if (!is.null(control_vars)) { X_ranger <- cbind(X_ranger, covariate_mat_model) }
          
          # ranger는 formula 인터페이스가 더 안정적임
          ranger_data <- data.frame(y = y, X_ranger)
          
          rf_model <- ranger::ranger(
            y ~ ., 
            data = ranger_data, 
            num.trees = 500, 
            importance = "impurity", # Gini
            ...
          )
          
          weights_magnitude_all <- ranger::importance(rf_model)
          
          gene_names_in_model <- names(weights_magnitude_all)[!names(weights_magnitude_all) %in% colnames(covariate_mat_model)]
          weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]
          
          top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
          weights_magnitude <- weights_magnitude_genes[top_genes]

          # 방향성 보정 및 performance 계산
          if (n_groups == 2) {
            g1_cells <- y == levels(y)[1]
            g2_cells <- y == levels(y)[2]
            mean_g1 <- colMeans(X_ranger[g1_cells, top_genes, drop=FALSE])
            mean_g2 <- colMeans(X_ranger[g2_cells, top_genes, drop=FALSE])
            effect_size <- mean_g2 - mean_g1 
            weights <- weights_magnitude * sign(effect_size)
          } else {
            warning("random_forest_ranger: n_groups > 2. Score represents magnitude (importance), not direction.")
            weights <- weights_magnitude
          }
          
          scores <- as.numeric(X_ranger[, top_genes] %*% weights)
          names(scores) <- rownames(X_ranger)
          
          pred <- rf_model$predictions
          if (n_groups == 2) {
            if (requireNamespace("pROC", quietly = TRUE)) {
              roc_obj <- pROC::roc(y, scores, quiet=TRUE)
              auc <- as.numeric(pROC::auc(roc_obj))
            } else { auc <- NA }
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
          } else {
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, confusion = table(pred, y))
          }
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) rf_model else NULL)
        },

        xgboost = { # [Req 4] 신규 (Beta)
          if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("xgboost package required.")
          }
          
          # X는 스케일링된 데이터 (expr_mat_scaled)
          X_xgb <- X 
          if (!is.null(control_vars)) { X_xgb <- cbind(X_xgb, covariate_mat_model) }
          
          # y를 0/1 숫자로 변환
          y_numeric <- as.numeric(y) - 1
          
          dtrain <- xgboost::xgb.DMatrix(data = X_xgb, label = y_numeric)
          
          params <- list(
            objective = if (n_groups == 2) "binary:logistic" else "multi:softmax",
            eval_metric = if (n_groups == 2) "logloss" else "mlogloss",
            nthread = 4,
            eta = 0.1
          )
          if (n_groups > 2) { params$num_class = n_groups }
          
          xgb_model <- xgboost::xgb.train(params, dtrain, nrounds = 100, ...)
          
          imp_matrix <- xgboost::xgb.importance(model = xgb_model)
          weights_magnitude_all <- imp_matrix$Gain
          names(weights_magnitude_all) <- imp_matrix$Feature
          
          gene_names_in_model <- names(weights_magnitude_all)[!names(weights_magnitude_all) %in% colnames(covariate_mat_model)]
          weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]
          
          top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
          weights_magnitude <- weights_magnitude_genes[top_genes]

          # 방향성 보정 및 performance 계산
          if (n_groups == 2) {
            g1_cells <- y == levels(y)[1]
            g2_cells <- y == levels(y)[2]
            mean_g1 <- colMeans(X_xgb[g1_cells, top_genes, drop=FALSE])
            mean_g2 <- colMeans(X_xgb[g2_cells, top_genes, drop=FALSE])
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
            pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
            if (requireNamespace("pROC", quietly = TRUE)) {
              roc_obj <- pROC::roc(y, scores, quiet=TRUE)
              auc <- as.numeric(pROC::auc(roc_obj))
            } else { auc <- NA }
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
          } else {
            pred_class_indices <- apply(pred_probs, 1, which.max)
            pred <- levels(y)[pred_class_indices]
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, confusion = table(pred, y))
          }
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) xgb_model else NULL)
        },

        # --- 2. Regularization ---

        lasso = {
          if (!requireNamespace("glmnet", quietly = TRUE)) {
            stop("glmnet package required. Install with: install.packages('glmnet')")
          }
          
          # X는 스케일링된 데이터 (expr_mat_scaled)
          if (is.null(control_vars)) {
            X_model <- X
            penalty_vec <- rep(1, ncol(X_model))
          } else {
            X_model <- cbind(X, covariate_mat_model)
            penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
          }
          
          # alpha=1 (LASSO)
          if (n_groups == 2) {
            cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=1, 
                                        penalty.factor = penalty_vec, ...)
          } else {
            cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=1, 
                                        penalty.factor = penalty_vec, ...)
          }
          
          coefs <- coef(cv_fit, s = lambda_selection) 
          
          if (n_groups == 2) {
            weights_all <- as.numeric(coefs[2:(ncol(X)+1)]) 
            names(weights_all) <- rownames(coefs)[2:(ncol(X)+1)]
          } else {
            coef_list <- lapply(coefs, function(x) as.numeric(x[2:(ncol(X)+1)]))
            weights_all <- rowMeans(do.call(cbind, coef_list))
            names(weights_all) <- rownames(coefs[[1]])[2:(ncol(X)+1)]
          }
          
          if (length(weights_all) == 0) {
            stop("LASSO returned no gene coefficients.")
          }

          top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
          weights <- weights_all[top_genes]
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          pred_probs <- predict(cv_fit, newx=X_model, s = lambda_selection, type="response")
          if (n_groups == 2) {
            pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
            if (requireNamespace("pROC", quietly = TRUE)) {
              roc_obj <- pROC::roc(y, scores, quiet=TRUE) 
              auc <- as.numeric(pROC::auc(roc_obj))
            } else { auc <- NA }
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
          } else {
            pred_class_indices <- apply(pred_probs[,,1], 1, which.max)
            pred <- levels(y)[pred_class_indices]
            acc <- mean(pred == y)
            perf <- list(accuracy = acc, confusion = table(pred, y))
          }
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) cv_fit else NULL)
        },
        
        ridge = { # [Req 4] 신규
          if (!requireNamespace("glmnet", quietly = TRUE)) { stop("glmnet required.") }
          
          # X는 스케일링된 데이터 (expr_mat_scaled)
          if (is.null(control_vars)) {
            X_model <- X
            penalty_vec <- rep(1, ncol(X_model))
          } else {
            X_model <- cbind(X, covariate_mat_model)
            penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
          }
          
          # alpha=0 (Ridge)
          if (n_groups == 2) {
            cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=0, 
                                        penalty.factor = penalty_vec, ...)
          } else {
            cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=0, 
                                        penalty.factor = penalty_vec, ...)
          }
          
          # ... (v4.1 lasso와 동일한 계수 추출, top_genes 선별, performance 로직) ...
          list(...)
        },
        
        elastic_net = { # [Req 4] 신규
          if (!requireNamespace("glmnet", quietly = TRUE)) { stop("glmnet required.") }
          
          # X는 스케일링된 데이터 (expr_mat_scaled)
          if (is.null(control_vars)) {
            X_model <- X
            penalty_vec <- rep(1, ncol(X_model))
          } else {
            X_model <- cbind(X, covariate_mat_model)
            penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
          }
          
          # alpha=enet.alpha (기본 0.5)
          if (n_groups == 2) {
            cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=enet.alpha, 
                                        penalty.factor = penalty_vec, ...)
          } else {
            cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=enet.alpha, 
                                        penalty.factor = penalty_vec, ...)
          }
          
          # ... (v4.1 lasso와 동일한 계수 추출, top_genes 선별, performance 로직) ...
          list(...)
        },

        # --- 3. Loadings ---

        pca_loadings = {
          # [Req 2] n_pcs 로직 적용
          
          # expr_mat_method는 보정O, 스케일링X (또는 스케일링O, 보정X)
          # v5.2는 prcomp 전에 다시 스케일링을 보장
          pca_res <- prcomp(X, center=TRUE, scale.=TRUE) 
          
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
          top_n_pcs_indices <- order(pc_cors, decreasing=TRUE)[1:min(pca.n_pcs, n_pcs_to_test)]
          
          # 3. 중요도 계산: Top N개 PC의 로딩 절대값 합
          top_pc_loadings <- pca_res$rotation[, top_n_pcs_indices, drop=FALSE]
          weights_magnitude_all <- rowSums(abs(top_pc_loadings))
          
          # 4. 방향성 계산: Best PC 1개의 부호있는 로딩
          weights_directional_all <- pca_res$rotation[, best_pc_idx]
          
          # 5. 유전자 선별: Top N개 기준
          top_genes <- names(sort(weights_magnitude_all, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_all))])
          
          # 6. 최종 가중치: Best PC 1개의 방향성 적용
          weights <- weights_directional_all[top_genes]
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          perf <- list(Best_PC = best_pc_idx, 
                       Top_N_PCs = top_n_pcs_indices,
                       Top_N_Correlations = pc_cors[top_n_pcs_indices],
                       Best_PC_VarExplained = summary(pca_res)$importance[2, best_pc_idx])
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) pca_res else NULL)
        },

        nmf_loadings = { # [Req 3] 이름 변경
          if (!requireNamespace("NMF", quietly = TRUE)) {
            stop("NMF package required.")
          }
          
          # expr_mat_method는 'expr_mat_base'
          expr_mat_nmf <- expr_mat_method 
          
          if (!is.null(control_vars)) {
            warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s...", m))
            covariate_mat <- model.matrix(~ . - 1, data = meta.data_clean[, control_vars, drop=FALSE])
            expr_mat_nmf <- limma::removeBatchEffect(expr_mat_nmf, covariates = covariate_mat)
          }
          
          expr_mat_pos <- expr_mat_nmf - min(expr_mat_nmf) + 0.01 
          
          rank <- min(n_groups + 2, 10)
          dots <- list(...)
          
          # [Req 6 FIX] '...' (dots)에서 NMF 유효 인자만 필터링
          valid_nmf_args <- c("nrun", ".options", ".pbackend") 
          filtered_dots <- dots[names(dots) %in% valid_nmf_args]
          
          nmf_args <- c(list(x = expr_mat_pos, 
                             rank = rank, 
                             seed = fgs_seed, 
                             method = "brunet"), # 기본 알고리즘 명시
                        filtered_dots)
          
          nmf_res <- do.call(NMF::nmf, nmf_args)
          
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
          
          top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, length(weights_magnitude))])
          weights_magnitude <- weights_magnitude[top_genes]
          
          if (n_groups == 2) {
            g1_cells <- y == levels(y)[1]
            g2_cells <- y == levels(y)[2]
            mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
            mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
            effect_size <- mean_g2 - mean_g1
            weights <- weights_magnitude * sign(effect_size)
          } else {
            warning("nmf_loadings: n_groups > 2. Score represents magnitude (importance), not direction.")
            weights <- weights_magnitude
          }
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          perf <- list(component = best_component, correlation = component_cors[best_component])
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) nmf_res else NULL)
        },

        # --- 4. Statistical Modelling ---
        
        gam = {
          if (!requireNamespace("mgcv", quietly = TRUE)) { ... }
          
          deviance_explained <- numeric(nrow(expr_mat_method))
          names(deviance_explained) <- rownames(expr_mat_method)
          
          gam_data <- data.frame(y_var_numeric = as.numeric(target_binary) - 1)
          if (!is.null(control_vars)) {
            gam_data <- cbind(gam_data, covariate_mat_model)
            formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse=" + "))
          } else {
            formula_base <- ""
          }
          
          # [Req 1] gam.k 인자 사용
          full_formula_str_base <- paste("y_var_numeric ~ s(gene_expr, k=", gam.k, ", bs='cr')", formula_base)

          convergence_warnings <- 0
          genes_skipped <- 0
          
          for (i in 1:nrow(expr_mat_method)) {
            gam_data$gene_expr <- expr_mat_method[i, ]
            
            # [Req 1] 고유 값 개수 필터링
            n_unique_vals <- length(unique(gam_data$gene_expr))
            if (n_unique_vals < gam.min_unique) {
              deviance_explained[i] <- 0
              genes_skipped <- genes_skipped + 1
              next
            }

            fit_result <- tryCatch({
              if (n_groups == 2) {
                # [Point 2] 'gam' 대신 'bam' 사용
                mgcv::bam(as.formula(full_formula_str_base), data = gam_data, family="binomial", ...)
              } else {
                mgcv::bam(as.formula(full_formula_str_base), data = gam_data, ...)
              }
            }, warning = function(w) {
              if (grepl("did not converge", w$message)) {
                convergence_warnings <<- convergence_warnings + 1
              }
              invokeRestart("muffleWarning")
            }, error = function(e) {
              warning(sprintf("GAM failed for gene %s: %s", rownames(expr_mat_method)[i], e$message))
              return(NULL) 
            })
            
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
          
          weights_magnitude <- deviance_explained
          top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0, na.rm=TRUE))])
          
          if (length(top_genes) == 0) {
            warning("GAM method found no genes with deviance > 0.")
            return(list(genes=character(0), weights=numeric(0), scores=numeric(0), performance=list()))
          }
          
          weights_magnitude <- weights_magnitude[top_genes]

          if (n_groups == 2) {
            g1_cells <- y == levels(y)[1]
            g2_cells <- y == levels(y)[2]
            mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
            mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
            effect_size <- mean_g2 - mean_g1
            weights <- weights_magnitude * sign(effect_size)
          } else {
            warning("gam: n_groups > 2. Score represents magnitude (importance), not direction.")
            weights <- weights_magnitude
          }
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          perf <- list(deviance_explained = weights_magnitude)
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = NULL)
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
            design <- model.matrix(~0 + target_binary_limma, data = meta.data_clean)
            colnames(design)[1:n_groups] <- levels(target_binary_limma)
          } else {
            control_formula <- paste(control_vars, collapse = " + ")
            full_formula <- as.formula(paste("~0 + target_binary_limma +", control_formula))
            design <- model.matrix(full_formula, data = meta.data_clean)
            colnames(design)[1:n_groups] <- levels(target_binary_limma) 
          }

          fit <- limma::lmFit(expr_mat_method, design)
          
          if (n_groups == 2) {
            contrast_str <- paste(levels(target_binary_limma)[2], levels(target_binary_limma)[1], sep="-")
          } else {
            contrast_pairs <- combn(levels(target_binary_limma), 2, function(x) paste(x[2], x[1], sep="-"))
            contrast_str <- contrast_pairs
            warning("limma: n_groups > 2. Using all pairwise contrasts. Weights based on average t-statistic.")
          }
          
          contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
          fit2 <- limma::contrasts.fit(fit, contrast_mat)
          fit2 <- limma::eBayes(fit2)
          
          if(n_groups == 2) {
              top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
              weights_all <- top_table$t
          } else {
              top_table <- limma::topTable(fit2, number=Inf, sort.by="F")
              weights_all <- rowMeans(abs(fit2$t[, contrast_str, drop=FALSE]))
          }
          
          names(weights_all) <- rownames(top_table)
          
          top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
          
          if(n_groups == 2) {
              weights <- top_table[top_genes, "t"]
          } else {
              weights <- weights_all[top_genes]
          }

          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          perf <- list(top_table = top_table[top_genes, ])
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = if(return_model) fit2 else NULL)
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
              
              test <- try(wilcox.test(group1, group2), silent=TRUE) 
              
              if(inherits(test, "try-error")) {
                pvals[i] <- 1.0
                effect_sizes[i] <- 0
              } else {
                pvals[i] <- test$p.value
                effect_sizes[i] <- median(group2) - median(group1)
              }
            } else {
              test <- try(kruskal.test(expr_mat_method[i, ] ~ target_binary), silent=TRUE) 
              
              if(inherits(test, "try-error")) {
                 pvals[i] <- 1.0
                 effect_sizes[i] <- 0
              } else {
                 pvals[i] <- test$p.value
                 effect_sizes[i] <- var(tapply(expr_mat_method[i, ], target_binary, median))
              }
            }
          }
          
          padj <- p.adjust(pvals, method="BH")
          
          ranks <- rank(pvals) + rank(-abs(effect_sizes))
          top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
          
          weights <- effect_sizes[top_genes]
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          perf <- list(pvalues = pvals[top_genes], padj = padj[top_genes],
                       effect_sizes = effect_sizes[top_genes])
          
          list(genes = top_genes, weights = weights, scores = scores,
               performance = perf, model = NULL)
        }
        
      ) # --- switch 끝 ---

      # === 5. 결과 저장 ===
      result$method <- m
      result$target_var <- target_var
      result$n_groups <- n_groups
      result$n_cells <- ncol(expr_mat_base)
      
      class(result) <- c("gene_signature", "list")
      results_list[[m]] <- result
      
    }, error = function(e) {
      warning(sprintf("Method '%s' failed with error: %s", m, e$message))
      results_list[[m]] <- list(method = m, error = e$message)
    })
    
  } # --- for 루프 끝 ---

  return(results_list)
}


#' Train a Meta-Learner (L2 Stacking Model)
#'
#' Trains multiple L2 candidate models on L1 scores and selects the best
#' one based on cross-validation performance.
#'
#' @param l1_signatures A named list of trained L1 signatures (outputs from
#'                      find_gene_signature). e.g., list(lasso=sig1, rf=sig2)
#' @param holdout_data A Seurat object or matrix (genes x cells) NOT used
#'                      for training L1 models.
#' @param target_var The target variable column name in holdout_data@meta.data.
#' @param l2_methods A vector of model methods supported by caret 
#'                   (e.g., c("glm", "ranger", "xgbTree", "svmRadial")).
#' @param k_folds Number of folds for cross-validation to select the best L2 model.
#' @param metric The metric to optimize (e.g., "AUC", "Accuracy"). 
#'               (Note: For "AUC", target must be binary factor).
#' @param fgs_seed Seed for reproducibility.
#'
#' @return A list containing:
#'    - $best_model: The final L2 model (a caret 'train' object) 
#'                     trained on all l2_train_df.
#'    - $model_comparison: Results from caret::resamples comparing L2 candidates.
#'    - $l2_train_df: The dataframe used for L2 training (L1 scores + target).
#'    - $l1_signatures: The provided L1 signatures (for reference).
#' @export
train_meta_learner_v1 <- function(l1_signatures,
                               holdout_data,
                               target_var,
                               l2_methods = c("glm", "ranger", "xgbTree"),
                               k_folds = 5,
                               metric = "AUC",
                               fgs_seed = 42) {
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("caret package required. Install with: install.packages('caret')")
  }
  set.seed(fgs_seed)

  # === 1. L2 Feature 생성 ===
  message("--- Generating L2 features from holdout_data ---")
  
  l2_features_list <- lapply(l1_signatures, function(sig) {
    score_signature(expr_data = holdout_data,
                    signature = sig,
                    normalize = TRUE) # 스케일 통일을 위해 정규화
  })
  
  l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  
  # L1 모델 이름이 없으면 기본 이름 할당
  if (is.null(names(l1_signatures))) {
    names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures))
  }
  colnames(l2_train_df) <- names(l1_signatures)

  # === 2. L2 훈련 데이터 준비 ===
  
  # Seurat 객체에서 메타데이터 추출
  if (inherits(holdout_data, "Seurat")) {
    meta.data <- holdout_data@meta.data
    # holdout_data의 셀 순서와 l2_train_df의 행 순서가 일치하는지 확인
    if(!all(rownames(l2_train_df) == rownames(meta.data))) {
        meta.data <- meta.data[rownames(l2_train_df), ]
    }
    l2_target <- meta.data[[target_var]]
  } else {
    # 행렬인 경우, target_var는 별도의 벡터로 제공되어야 함 (이 예제에서는 단순화)
    stop("holdout_data is not a Seurat object. Please provide target_var as a separate vector matching cell order.")
  }

  # caret은 factor형 타겟 변수를 선호 (특히 분류 및 AUC 계산 시)
  if (!is.factor(l2_target)) {
    l2_target <- factor(l2_target)
  }
  
  # caret이 R 변수명으로 부적합한 이름을 처리하도록 함 (예: "limma-fold")
  colnames(l2_train_df) <- make.names(colnames(l2_train_df))
  
  # === 3. L2 모델 교차 검증 및 선택 ===
  message("--- Training and comparing L2 models via k-fold CV ---")

  # CV 설정
  train_control <- caret::trainControl(
    method = "cv",
    number = k_folds,
    summaryFunction = caret::twoClassSummary, # AUC, Sens, Spec 계산
    classProbs = TRUE, # AUC 계산을 위해 확률 예측 필요
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  # metric이 AUC가 아니면 기본 summary 사용
  if (metric != "AUC") {
    train_control$summaryFunction <- caret::defaultSummary
  }
  
  # 여러 L2 모델을 훈련시키기 위한 리스트
  model_list <- list()
  
  for (method in l2_methods) {
    tryCatch({
      message(sprintf("Training L2 candidate: %s", method))
      model_list[[method]] <- caret::train(
        x = l2_train_df,
        y = l2_target,
        method = method,
        trControl = train_control,
        metric = metric,
        tuneLength = 5 # 각 모델의 기본 하이퍼파라미터 튜닝
      )
    }, error = function(e) {
      warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
    })
  }

  if (length(model_list) == 0) {
    stop("No L2 models were successfully trained.")
  }

  # 훈련된 모델들의 성능 비교
  model_comparison <- caret::resamples(model_list)
  print(summary(model_comparison))

  # === 4. Best 모델 선택 및 재훈련 (caret이 자동으로 처리) ===
  # caret::train은 CV를 사용해 최적의 하이퍼파라미터를 찾고,
  # 그 파라미터로 "전체 L2 데이터"를 재훈련시킨 모델($finalModel)을 반환함.
  # 우리는 이 훈련된 모델들 중 'Best'를 선택하기만 하면 됨.
  
  best_model_name <- names(which.max(
    sapply(model_list, function(m) max(m$results[, metric]))
  ))
  
  best_model <- model_list[[best_model_name]]
  
  message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
                  best_model_name, metric, max(best_model$results[, metric])))
  
  # === 5. 결과 반환 ===
  return(list(
    best_model = best_model,
    best_model_name = best_model_name,
    model_comparison = model_comparison,
    l2_train_df = data.frame(l2_train_df, target = l2_target),
    l1_signatures = l1_signatures
  ))
}


# ---
# 2. train_meta_learner (수정본)
#    - (FIX 1): GetAssayData 1회 호출
#    - (FIX 2): NA target 명시적 제거
# ---
# (caret, ranger, xgboost 등 필요 패키지 로드 가정)
#' Train a Meta-Learner (L2 Stacking Model)
#'
#' Trains multiple L2 candidate models on L1 scores and selects the best
#' one based on cross-validation performance.
#'
#' @param l1_signatures A named list of trained L1 signatures (outputs from
#'                      find_gene_signature). e.g., list(lasso=sig1, rf=sig2)
#' @param holdout_data A Seurat object or matrix (genes x cells) NOT used
#'                      for training L1 models.
#' @param target_var The target variable column name in holdout_data@meta.data.
#' @param l2_methods A vector of model methods supported by caret 
#'                   (e.g., c("glm", "ranger", "xgbTree", "svmRadial")).
#' @param k_folds Number of folds for cross-validation to select the best L2 model.
#' @param metric The metric to optimize (e.g., "AUC", "Accuracy"). 
#'               (Note: For "AUC", target must be binary factor).
#' @param fgs_seed Seed for reproducibility.
#'
#' @return A list containing:
#'    - $best_model: The final L2 model (a caret 'train' object) 
#'                     trained on all l2_train_df.
#'    - $model_comparison: Results from caret::resamples comparing L2 candidates.
#'    - $l2_train_df: The dataframe used for L2 training (L1 scores + target).
#'    - $l1_signatures: The provided L1 signatures (for reference).
#' @export
train_meta_learner_v2 <- function(l1_signatures,
                               holdout_data,
                               target_var,
                               l2_methods = c("glm", "ranger", "xgbTree"),
                               k_folds = 5,
                               metric = "AUC",
                               fgs_seed = 42) {
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("caret package required.")
  }
  set.seed(fgs_seed)

  # === 1. 데이터 추출 (PERFORMANCE FIX 1) ===
  message("--- Preparing holdout data (Extracting once) ---")
  if (inherits(holdout_data, "Seurat")) {
    # GetAssayData를 루프 밖에서 단 한 번만 호출 (희소 행렬로)
    expr_mat <- Seurat::GetAssayData(holdout_data, layer="data")
    meta_data <- holdout_data@meta.data
  } else {
    # Seurat 객체가 아니면, holdout_data가 이미 (희소) 행렬이라고 가정
    expr_mat <- holdout_data
    # 이 경우 target_var는 반드시 메타데이터 프레임이거나 벡터여야 함 (로직 단순화)
    if (!is.vector(target_var) && !is.factor(target_var)) {
        stop("If holdout_data is a matrix, target_var must be the actual target vector.")
    }
    meta_data <- NULL # 메타데이터 객체가 없음
  }

  # === 2. L2 Feature 생성 ===
  message("--- Generating L2 features from holdout data ---")
  
  l2_features_list <- lapply(l1_signatures, function(sig) {
    # Seurat 객체 대신 추출된 'expr_mat'를 전달
    score_signature(expr_data = expr_mat, 
                    signature = sig,
                    normalize = TRUE) 
  })
  
  l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  
  if (is.null(names(l1_signatures))) {
    names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures))
  }
  # caret이 인식하도록 R 변수명으로 정리
  colnames(l2_train_df) <- make.names(names(l1_signatures))

  # === 3. L2 타겟 변수 준비 ===
  if (!is.null(meta_data)) {
    # l2_train_df의 행 순서(cell ID)와 meta_data의 행 순서가 일치하는지 확인/정렬
    if (!all(rownames(l2_train_df) == rownames(meta_data))) {
        meta_data <- meta_data[rownames(l2_train_df), , drop = FALSE]
    }
    l2_target <- meta_data[[target_var]]
  } else {
    l2_target <- target_var # 이미 벡터라고 가정
    if (length(l2_target) != nrow(l2_train_df)) {
      stop("Target vector length does not match expression data cell count.")
    }
  }

  if (!is.factor(l2_target)) {
    l2_target <- factor(l2_target)
  }
  
  # === 4. NA 타겟 제거 (BUG FIX 2) ===
  message(sprintf("Initial L2 target length: %d", length(l2_target)))
  na_indices <- is.na(l2_target)
  
  if (any(na_indices)) {
    n_na <- sum(na_indices)
    message(sprintf("Found and removing %d cells with NA target.", n_na))
    
    l2_target_filtered <- l2_target[!na_indices]
    l2_train_df_filtered <- l2_train_df[!na_indices, , drop = FALSE]
    
  } else {
    l2_target_filtered <- l2_target
    l2_train_df_filtered <- l2_train_df
  }
  
  message(sprintf("Final L2 training set size: %d cells", length(l2_target_filtered)))

  if (length(l2_target_filtered) == 0 || nrow(l2_train_df_filtered) == 0) {
      stop("No valid data remaining after removing NAs from target variable.")
  }

  # === 5. L2 모델 교차 검증 및 선택 ===
  message("--- Training and comparing L2 models via k-fold CV ---")

  # 타겟 유형(이진/다중)에 따라 CV 컨트롤 분기
  if (length(levels(l2_target_filtered)) == 2) {
    message("Binary target detected. Using twoClassSummary (for AUC).")
    train_control <- caret::trainControl(
      method = "cv", number = k_folds,
      summaryFunction = caret::twoClassSummary,
      classProbs = TRUE,
      savePredictions = "final",
      allowParallel = TRUE
    )
    # 이진 분류 시 metric이 AUC가 아니면 AUC로 강제
    if (metric != "AUC") {
      message("NOTE: Target is binary, forcing metric to 'AUC'.")
      metric <- "AUC"
    }
  } else {
    message("Multi-class target detected. Using defaultSummary.")
    train_control <- caret::trainControl(
      method = "cv", number = k_folds,
      summaryFunction = caret::defaultSummary,
      savePredictions = "final",
      allowParallel = TRUE
    )
    # 다중 분류 시 AUC 사용 불가
    if (metric == "AUC") {
      message("WARNING: Target is multi-class, 'AUC' metric is invalid. Defaulting to 'Accuracy'.")
      metric <- "Accuracy"
    }
  }
  
  model_list <- list()
  
  for (method in l2_methods) {
    tryCatch({
      message(sprintf("Training L2 candidate: %s", method))
      model_list[[method]] <- caret::train(
        x = l2_train_df_filtered, # <-- 필터링된 데이터 사용
        y = l2_target_filtered,   # <-- 필터링된 타겟 사용
        method = method,
        trControl = train_control,
        metric = metric,
        tuneLength = 5 # 기본 하이퍼파라미터 튜닝
      )
    }, error = function(e) {
      warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
    })
  }

  if (length(model_list) == 0) {
    stop("No L2 models were successfully trained.")
  }

  # 훈련된 모델들의 성능 비교
  model_comparison <- caret::resamples(model_list)
  print(summary(model_comparison))

  # Best 모델 선택
  best_metric_values <- sapply(model_list, function(m) max(m$results[, metric], na.rm = TRUE))
  best_model_name <- names(which.max(best_metric_values))
  best_model <- model_list[[best_model_name]]
  
  message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
                  best_model_name, metric, max(best_model$results[, metric], na.rm = TRUE)))
  
  # === 6. 결과 반환 ===
  return(list(
    best_model = best_model,
    best_model_name = best_model_name,
    model_comparison = model_comparison,
    l2_train_df_filtered = data.frame(l2_train_df_filtered, target = l2_target_filtered),
    l1_signatures = l1_signatures
  ))
}

#' Train a Meta-Learner (L2 Stacking Model)
#'
#' Trains multiple L2 candidate models on L1 scores and selects the best
#' one based on cross-validation performance.
#'
#' @param l1_signatures A named list of trained L1 signatures (outputs from
#'                      find_gene_signature). e.g., list(lasso=sig1, rf=sig2)
#' @param holdout_data A Seurat object or matrix (genes x cells) NOT used
#'                      for training L1 models.
#' @param target_var The target variable column name in holdout_data@meta.data.
#' @param l2_methods A vector of model methods supported by caret 
#'                   (e.g., c("glm", "ranger", "xgbTree", "svmRadial")).
#' @param k_folds Number of folds for cross-validation to select the best L2 model.
#' @param metric The metric to optimize (e.g., "AUC", "Accuracy"). 
#'               (Note: For "AUC", target must be binary factor).
#' @param fgs_seed Seed for reproducibility.
#'
#' @return A list containing:
#'    - $best_model: The final L2 model (a caret 'train' object) 
#'                     trained on all l2_train_df.
#'    - $model_comparison: Results from caret::resamples comparing L2 candidates.
#'    - $l2_train_df: The dataframe used for L2 training (L1 scores + target).
#'    - $l1_signatures: The provided L1 signatures (for reference).
#' @export
train_meta_learner_v3 <- function(l1_signatures,
                               holdout_data,
                               target_var,
                               l2_methods = c("glm", "ranger", "xgbTree"),
                               k_folds = 5,
                               metric = "AUC",
                               fgs_seed = 42) {
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("caret package required.")
  }
  set.seed(fgs_seed)

  # === 1. 데이터 추출 (PERFORMANCE FIX 1) ===
  message("--- Preparing holdout data (Extracting once) ---")
  if (inherits(holdout_data, "Seurat")) {
    expr_mat <- Seurat::GetAssayData(holdout_data, layer="data")
    meta_data <- holdout_data@meta.data
  } else {
    expr_mat <- holdout_data
    if (!is.vector(target_var) && !is.factor(target_var)) {
        stop("If holdout_data is a matrix, target_var must be the actual target vector.")
    }
    meta_data <- NULL 
  }

  # === 2. L2 Feature 생성 ===
  message("--- Generating L2 features from holdout data ---")
  
  l2_features_list <- lapply(l1_signatures, function(sig) {
    score_signature(expr_data = expr_mat, 
                    signature = sig,
                    normalize = TRUE) 
  })
  
  l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  
  if (is.null(names(l1_signatures))) {
    names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures))
  }
  colnames(l2_train_df) <- make.names(names(l1_signatures))

  # === 3. L2 타겟 변수 준비 ===
  if (!is.null(meta_data)) {
    if (!all(rownames(l2_train_df) == rownames(meta_data))) {
        meta_data <- meta_data[rownames(l2_train_df), , drop = FALSE]
    }
    l2_target <- meta_data[[target_var]]
  } else {
    l2_target <- target_var 
  }

  if (!is.factor(l2_target)) {
    l2_target <- factor(l2_target)
  }
  
  # === (BUG FIX 3): caret을 위한 Factor Level 이름 소독 ===
  original_levels <- levels(l2_target)
  sanitized_levels <- make.names(original_levels)
  
  if (!all(original_levels == sanitized_levels)) {
    message(sprintf("Sanitizing target levels for caret: '%s' -> '%s'", 
                    paste(original_levels, collapse="' / '"), 
                    paste(sanitized_levels, collapse="' / '")))
    levels(l2_target) <- sanitized_levels
  }
  # === End of fix ===

  # === 4. NA 타겟 제거 (BUG FIX 2) ===
  message(sprintf("Initial L2 target length: %d", length(l2_target)))
  na_indices <- is.na(l2_target)
  
  if (any(na_indices)) {
    n_na <- sum(na_indices)
    message(sprintf("Found and removing %d cells with NA target.", n_na))
    
    l2_target_filtered <- l2_target[!na_indices]
    l2_train_df_filtered <- l2_train_df[!na_indices, , drop = FALSE]
    
  } else {
    l2_target_filtered <- l2_target
    l2_train_df_filtered <- l2_train_df
  }
  
  message(sprintf("Final L2 training set size: %d cells", length(l2_target_filtered)))

  if (length(l2_target_filtered) == 0 || nrow(l2_train_df_filtered) == 0) {
      stop("No valid data remaining after removing NAs from target variable.")
  }

  # === 5. L2 모델 교차 검증 및 선택 ===
  message("--- Training and comparing L2 models via k-fold CV ---")

  if (length(levels(l2_target_filtered)) == 2) {
    message("Binary target detected. Using twoClassSummary (for AUC).")
    train_control <- caret::trainControl(
      method = "cv", number = k_folds,
      summaryFunction = caret::twoClassSummary,
      classProbs = TRUE,
      savePredictions = "final",
      allowParallel = TRUE
    )
    if (metric != "AUC") {
      message("NOTE: Target is binary, forcing metric to 'AUC'.")
      metric <- "AUC"
    }
  } else {
    message("Multi-class target detected. Using defaultSummary.")
    train_control <- caret::trainControl(
      method = "cv", number = k_folds,
      summaryFunction = caret::defaultSummary,
      savePredictions = "final",
      allowParallel = TRUE
    )
    if (metric == "AUC") {
      message("WARNING: Target is multi-class, 'AUC' metric is invalid. Defaulting to 'Accuracy'.")
      metric <- "Accuracy"
    }
  }
  
  model_list <- list()
  
  for (method in l2_methods) {
    tryCatch({
      message(sprintf("Training L2 candidate: %s", method))
      model_list[[method]] <- caret::train(
        x = l2_train_df_filtered, 
        y = l2_target_filtered,   
        method = method,
        trControl = train_control,
        metric = metric,
        tuneLength = 5 
      )
    }, error = function(e) {
      warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
    })
  }

  if (length(model_list) == 0) {
    stop("No L2 models were successfully trained.")
  }

  model_comparison <- caret::resamples(model_list)
  print(summary(model_comparison))

  best_metric_values <- sapply(model_list, function(m) max(m$results[, metric], na.rm = TRUE))
  best_model_name <- names(which.max(best_metric_values))
  best_model <- model_list[[best_model_name]]
  
  message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
                  best_model_name, metric, max(best_model$results[, metric], na.rm = TRUE)))
  
  # === 6. 결과 반환 ===
  return(list(
    best_model = best_model,
    best_model_name = best_model_name,
    model_comparison = model_comparison,
    l2_train_df_filtered = data.frame(l2_train_df_filtered, target = l2_target_filtered),
    l1_signatures = l1_signatures
  ))
}

# (score_signature 함수는 이전의 'Robust Fix' 버전을 사용한다고 가정)
# (caret, ranger, xgboost 등 필요 패키지 로드 가정)

train_meta_learner_v4 <- function(l1_signatures,
                               holdout_data,
                               target_var,
                               l2_methods = c("glm", "ranger", "xgbTree"),
                               k_folds = 5,
                               metric = "AUC",
                               fgs_seed = 42) {
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("caret package required.")
  }
  set.seed(fgs_seed)

  # === 1. 데이터 추출 ===
  message("--- Preparing holdout data (Extracting once) ---")
  if (inherits(holdout_data, "Seurat")) {
    expr_mat <- Seurat::GetAssayData(holdout_data, layer="data")
    meta_data <- holdout_data@meta.data
  } else {
    expr_mat <- holdout_data
    if (!is.vector(target_var) && !is.factor(target_var)) {
        stop("If holdout_data is a matrix, target_var must be the actual target vector.")
    }
    meta_data <- NULL 
  }

  # === 2. L2 Feature 생성 ===
  message("--- Generating L2 features from holdout data ---")
  l2_features_list <- lapply(l1_signatures, function(sig) {
    score_signature(expr_data = expr_mat, signature = sig, normalize = TRUE) 
  })
  l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  if (is.null(names(l1_signatures))) { names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures)) }
  colnames(l2_train_df) <- make.names(names(l1_signatures))

  # === 3. L2 타겟 변수 준비 및 소독 (BUG FIX 3) ===
  if (!is.null(meta_data)) {
    if (!all(rownames(l2_train_df) == rownames(meta_data))) {
        meta_data <- meta_data[rownames(l2_train_df), , drop = FALSE]
    }
    l2_target <- meta_data[[target_var]]
  } else { l2_target <- target_var }
  if (!is.factor(l2_target)) { l2_target <- factor(l2_target) }
  
  original_levels <- levels(l2_target)
  sanitized_levels <- make.names(original_levels)
  if (!all(original_levels == sanitized_levels)) {
    message(sprintf("Sanitizing target levels for caret: '%s' -> '%s'", 
                    paste(original_levels, collapse="' / '"), 
                    paste(sanitized_levels, collapse="' / '")))
    levels(l2_target) <- sanitized_levels
  }

  # === 4. NA 타겟 제거 (BUG FIX 2) ===
  message(sprintf("Initial L2 target length: %d", length(l2_target)))
  na_target_indices <- is.na(l2_target)
  if (any(na_target_indices)) {
    n_na <- sum(na_target_indices)
    message(sprintf("Found and removing %d cells with NA target.", n_na))
    l2_target <- l2_target[!na_target_indices]
    l2_train_df <- l2_train_df[!na_target_indices, , drop = FALSE]
  }
  
  # === 5. NA/NaN/Inf Feature 제거 (BUG FIX 4) ===
  message("Checking for NA/NaN/Inf in L2 features (scores)...")
  complete_feature_indices <- complete.cases(l2_train_df)
  if (!all(complete_feature_indices)) {
    n_na_features <- sum(!complete_feature_indices)
    message(sprintf("Found and removing %d cells with NA/NaN/Inf features.", n_na_features))
    l2_target_filtered <- l2_target[complete_feature_indices]
    l2_train_df_filtered <- l2_train_df[complete_feature_indices, , drop = FALSE]
  } else {
    l2_target_filtered <- l2_target
    l2_train_df_filtered <- l2_train_df
  }

  # === (BUG FIX 5): Zero-Variance Predictor 제거 ===
  message("Checking for zero-variance predictors (constant features)...")
  nzv_indices <- caret::nearZeroVar(l2_train_df_filtered, saveMetrics = FALSE)
  
  if (length(nzv_indices) > 0) {
    n_nzv <- length(nzv_indices)
    nzv_names <- colnames(l2_train_df_filtered)[nzv_indices]
    message(sprintf("Found and removing %d constant features: %s",
                    n_nzv, paste(nzv_names, collapse=", ")))
    
    l2_train_df_filtered <- l2_train_df_filtered[, -nzv_indices, drop = FALSE]
  }
  # === End of fix ===

  message(sprintf("Final L2 training set size: %d cells, %d features", 
                  length(l2_target_filtered), ncol(l2_train_df_filtered)))
                  
  if (length(l2_target_filtered) == 0 || nrow(l2_train_df_filtered) == 0) {
      stop("No valid data remaining after removing NAs from target and/or features.")
  }
  if (ncol(l2_train_df_filtered) == 0) {
      stop("No valid features remaining after removing zero-variance predictors.")
  }

  # === 6. L2 모델 교차 검증 및 선택 ===
  message("--- Training and comparing L2 models via k-fold CV ---")
  
  if (length(levels(l2_target_filtered)) == 2) {
    message("Binary target detected. Using twoClassSummary (for AUC).")
    train_control <- caret::trainControl(
      method = "cv", number = k_folds,
      summaryFunction = caret::twoClassSummary,
      classProbs = TRUE,
      savePredictions = "final",
      allowParallel = TRUE
    )
    if (metric != "AUC") { message("NOTE: Target is binary, forcing metric to 'AUC'."); metric <- "AUC" }
  } else {
    message("Multi-class target detected. Using defaultSummary.")
    train_control <- caret::trainControl(method = "cv", number = k_folds, summaryFunction = caret::defaultSummary, savePredictions = "final", allowParallel = TRUE)
    if (metric == "AUC") { message("WARNING: Target is multi-class, 'AUC' metric is invalid. Defaulting to 'Accuracy'."); metric <- "Accuracy" }
  }
  
  model_list <- list()
  for (method in l2_methods) {
    tryCatch({
      message(sprintf("Training L2 candidate: %s", method))
      model_list[[method]] <- caret::train(
        x = l2_train_df_filtered,
        y = l2_target_filtered,
        method = method,
        trControl = train_control,
        metric = metric,
        tuneLength = 5
        # 'na.action' 인수 제거됨
      )
    }, error = function(e) {
      warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
    })
  }

  if (length(model_list) < 2) {
    if (length(model_list) == 0) {
        stop("No L2 models were successfully trained.")
    }
    warning("Only one L2 model was successfully trained. Comparison skipped.")
    model_comparison <- NULL
    best_model <- model_list[[1]]
    best_model_name <- names(model_list)[1]
    
  } else {
    model_comparison <- caret::resamples(model_list)
    print(summary(model_comparison))
    best_metric_values <- sapply(model_list, function(m) max(m$results[, metric], na.rm = TRUE))
    best_model_name <- names(which.max(best_metric_values))
    best_model <- model_list[[best_model_name]]
  }
  
  # 최종 에러 방지 (best_model이 메트릭 생성에 실패한 경우)
  best_metric_val <- tryCatch({
    max(best_model$results[, metric], na.rm = TRUE)
  }, error = function(e) {
    warning(sprintf("Could not extract metric '%s' from best model '%s'.", metric, best_model_name))
    return(NA_real_)
  })

  message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
                  best_model_name, metric, best_metric_val))
  
  # === 7. 결과 반환 ===
  return(list(
    best_model = best_model,
    best_model_name = best_model_name,
    model_comparison = model_comparison,
    l2_train_df_filtered = data.frame(l2_train_df_filtered, target = l2_target_filtered),
    l1_signatures = l1_signatures
  ))
}

#' @export
train_meta_learner_v5 <- function(...) {
  .Deprecated("TML6", package = "myR")
  TML6(...)
}

#' @export
train_meta_learner_v6 <- function(...) {
  .Deprecated("TML6", package = "myR")
  TML6(...)
}

#' Train Meta-Learner (TML6): Stacked Ensemble Model for Signature Scores
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
#'   (default: `c("glm", "ranger", "xgbTree")`). See `caret::train()` for
#'   available methods.
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
#'
#' @return A list with components:
#'   \describe{
#'     \item{best_model}{Caret `train` object for the selected L2 model.}
#'     \item{best_model_name}{Character string identifier of the best model.}
#'     \item{best_metric_name}{Name of the metric used for selection.}
#'     \item{model_comparison}{A `resamples` object comparing all trained models
#'           (or `NULL` if only one model succeeded).}
#'     \item{l2_train}{Data frame with signature scores (columns) and target
#'           variable (`.target` column).}
#'     \item{l1_signatures}{Standardized signatures (named numeric weight vectors).}
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
#' @examples
#' \dontrun{
#' # Example: Combine multiple signature-finding methods
#' sigs <- list(
#'   lasso = find_gene_signature(data, target_var="g3", method="lasso"),
#'   rf = find_gene_signature(data, target_var="g3", method="tree_based"),
#'   limma = find_gene_signature(data, target_var="g3", method="limma")
#' )
#'
#' # Train meta-learner on holdout data
#' meta_model <- TML6(
#'   l1_signatures = sigs,
#'   holdout_data = seurat_obj,
#'   target_var = "g3",
#'   l2_methods = c("glm", "ranger"),
#'   metric = "AUC"
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
TML6 <- function(
  l1_signatures,
  holdout_data,
  target_var,
  l2_methods = c("glm","ranger","xgbTree"),
  k_folds   = 5,
  metric    = c("AUC","ROC","Accuracy","Kappa"),
  fgs_seed  = 42,
  layer     = "data",
  allow_parallel = FALSE,
  parallel_workers = NULL
){
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  # ---- NEW: 표준화 유틸 ----
  as_signature <- function(sig){
    # 문자형 유전자 벡터
    if (is.character(sig)) return(structure(rep(1, length(sig)), names = sig))

    # 이름달린 numeric (가중치)
    if (is.numeric(sig) && !is.null(names(sig))) return(sig)

    # list(up=..., down=...)
    if (is.list(sig) && all(c("up","down") %in% names(sig))) {
      up <- sig$up; down <- sig$down
      stopifnot(is.character(up) || is.null(up), is.character(down) || is.null(down))
      up   <- up   %||% character()
      down <- down %||% character()
      nm <- c(up, down)
      wt <- c(rep(1, length(up)), rep(-1, length(down)))
      return(structure(wt, names = nm))
    }

    # list(genes=..., weights=...)
    if (is.list(sig) && all(c("genes","weights") %in% names(sig))) {
      g <- sig$genes; w <- sig$weights
      stopifnot(is.character(g), is.numeric(w), length(g) == length(w))
      names(w) <- g
      return(w)
    }

    # data.frame/tibble: gene/feature + weight/w/score
    if (is.data.frame(sig)) {
      cn <- tolower(colnames(sig))
      gene_col   <- which(cn %in% c("gene","genes","feature","features","symbol","id"))[1]
      weight_col <- which(cn %in% c("weight","weights","w","score","scores","coef","coefs"))[1]
      if (!is.na(gene_col) && !is.na(weight_col)) {
        g <- as.character(sig[[gene_col]])
        w <- as.numeric(sig[[weight_col]])
        ok <- !is.na(g) & !is.na(w)
        w  <- w[ok]; g <- g[ok]
        names(w) <- g
        return(w)
      }
    }

    stop("Unsupported signature format.")
  }

  .score_signature <- function(expr_data, signature, normalize = TRUE) {
    # signature를 표준화해 항상 "이름달린 가중치 numeric"으로 맞춘다
    w <- as_signature(signature)

    genes <- intersect(rownames(expr_data), names(w))
    if (length(genes) == 0) stop("Signature has no overlap with expression data.")

    ww <- w[genes]
    # (표준) 가중평균
    s <- as.numeric(Matrix::t(expr_data[genes, , drop = FALSE]) %*% ww)
    den <- sum(abs(ww))
    if (!is.finite(den) || den == 0) den <- length(ww)
    s <- s / den
    if (normalize) s <- as.numeric(scale(s))
    s
  }

  .is_binary <- function(y) length(levels(y)) == 2

  metric_choices <- c("AUC", "ROC", "Accuracy", "Kappa")

  .metric_map <- function(user_metric, is_binary) {
    if (is_binary) {
      # caret twoClassSummary → columns: ROC, Sens, Spec
      if (user_metric %in% c("AUC","ROC")) return(list(train_metric="ROC",  summary="twoClassSummary"))
      if (user_metric %in% c("Accuracy","Kappa")) return(list(train_metric=user_metric, summary="twoClassSummary"))
    } else {
      # 멀티클래스 → defaultSummary (Accuracy, Kappa)
      if (user_metric == "AUC") {
        message("WARNING: 'AUC/ROC' not defined for multi-class defaultSummary. Falling back to 'Accuracy'.")
        return(list(train_metric="Accuracy", summary="defaultSummary"))
      }
      return(list(train_metric=user_metric, summary="defaultSummary"))
    }
    list(train_metric = if (is_binary) "ROC" else "Accuracy",
         summary      = if (is_binary) "twoClassSummary" else "defaultSummary")
  }

  .package_ok <- function(pkg){
    suppressWarnings(suppressMessages(requireNamespace(pkg, quietly = TRUE)))
  }

  if (!.package_ok("caret")) stop("caret package required.")
  if (!.package_ok("Matrix")) stop("Matrix package required.")

  register_sequential <- function() {
    if (.package_ok("foreach")) {
      try(foreach::registerDoSEQ(), silent = TRUE)
    }
  }

  allow_parallel <- isTRUE(allow_parallel)
  worker_count <- parallel_workers
  if (is.null(worker_count)) {
    fallback <- getOption("mylit.meta_learner.workers", 4L)
    cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    cores <- cores %||% fallback
    worker_count <- max(1L, min(8L, as.integer(cores)))
  }

  if (allow_parallel) {
    if (.package_ok("future") && .package_ok("doFuture")) {
      oplan <- future::plan()
      on.exit({
        try(future::plan(oplan), silent = TRUE)
      }, add = TRUE)

      message(sprintf("Enabling parallel training with %d workers (future::multisession).", worker_count))
      tryCatch(
        future::plan(future::multisession, workers = worker_count),
        error = function(e) {
          warning("Failed to configure future multisession plan: ", conditionMessage(e), call. = FALSE)
          allow_parallel <<- FALSE
        }
      )

      if (allow_parallel) {
        doFuture::registerDoFuture()
        on.exit(register_sequential(), add = TRUE)
      } else {
        register_sequential()
      }

    } else {
      message("NOTE: Parallel dependencies (future/doFuture) not available; falling back to sequential execution.")
      allow_parallel <- FALSE
      register_sequential()
    }
  } else {
    register_sequential()
  }

  set.seed(fgs_seed)
  RNGkind(sample.kind = "Rejection")

  # --- 데이터 추출 ---
  if (inherits(holdout_data, "Seurat")) {
    if (!.package_ok("Seurat")) stop("Seurat package required for Seurat input.")
    expr_mat <- Seurat::GetAssayData(holdout_data, layer = layer)
    meta_data <- holdout_data@meta.data
    if (!target_var %in% colnames(meta_data))
      stop(sprintf("Target column '%s' not found in Seurat meta.data.", target_var))
    l2_target <- meta_data[[target_var]]
  } else {
    if (is.data.frame(holdout_data)) holdout_data <- as.matrix(holdout_data)
    if (!is.matrix(holdout_data)) stop("holdout_data must be Seurat or a matrix/data.frame [features x cells].")
    expr_mat <- holdout_data
    l2_target <- target_var
  }
  if (is.null(rownames(expr_mat))) stop("Expression matrix must have rownames (feature IDs).")
  if (is.null(colnames(expr_mat))) colnames(expr_mat) <- paste0("cell_", seq_len(ncol(expr_mat)))

  # --- L2 특성 구성 (시그니처 점수 계산) ---
  if (is.null(names(l1_signatures))) names(l1_signatures) <- paste0("L1_", seq_along(l1_signatures))

  l2_features_list <- lapply(l1_signatures, function(sig) .score_signature(expr_mat, sig, normalize = TRUE))
  l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  colnames(l2_train_df) <- make.names(names(l1_signatures))

  # --- 타깃 준비/정리 ---
  if (!is.factor(l2_target)) l2_target <- factor(l2_target)
  orig_lvls <- levels(l2_target)
  safe_lvls <- make.names(orig_lvls)
  if (!all(orig_lvls == safe_lvls)) {
    message(sprintf("Sanitizing target levels: '%s' -> '%s'",
                    paste(orig_lvls, collapse="'/'"),
                    paste(safe_lvls, collapse="'/'")))
    levels(l2_target) <- safe_lvls
  }

  keep <- !is.na(l2_target)
  if (!all(keep)) {
    message(sprintf("Removing %d rows with NA target.", sum(!keep)))
    l2_target <- l2_target[keep]
    l2_train_df <- l2_train_df[keep, , drop = FALSE]
  }

  row_ok <- stats::complete.cases(l2_train_df) &
            apply(l2_train_df, 1, function(r) all(is.finite(r)))
  if (!all(row_ok)) {
    message(sprintf("Removing %d rows with NA/NaN/Inf features.", sum(!row_ok)))
    l2_target <- l2_target[row_ok]
    l2_train_df <- l2_train_df[row_ok, , drop = FALSE]
  }

  nzv <- caret::nearZeroVar(l2_train_df, saveMetrics = FALSE)
  if (length(nzv) > 0) {
    message(sprintf("Removing %d zero-variance features: %s",
                    length(nzv), paste(colnames(l2_train_df)[nzv], collapse=", ")))
    l2_train_df <- l2_train_df[, -nzv, drop = FALSE]
  }

  if (nrow(l2_train_df) == 0 || ncol(l2_train_df) == 0)
    stop("No usable data remains after cleaning (check signatures and target).")

  # --- metric 매핑 & trainControl ---
  is_bin <- .is_binary(l2_target)
  metric <- match.arg(metric, metric_choices)
  map <- .metric_map(metric, is_bin)
  caret_metric <- map$train_metric

  ctrl <- caret::trainControl(
    method = "cv", number = k_folds,
    classProbs = is_bin,
    summaryFunction = if (map$summary == "twoClassSummary") caret::twoClassSummary else caret::defaultSummary,
    savePredictions = "final",
    allowParallel = allow_parallel
  )

  if ("xgbTree" %in% l2_methods) {
    if (!.package_ok("xgboost")) {
      warning("xgboost not available; dropping 'xgbTree'.")
      l2_methods <- setdiff(l2_methods, "xgbTree")
    } else {
      ok <- TRUE
      tryCatch(utils::packageVersion("xgboost"), error = function(e) ok <<- FALSE)
      if (!ok) {
        warning("xgboost seems broken; dropping 'xgbTree'.")
        l2_methods <- setdiff(l2_methods, "xgbTree")
      }
    }
  }
  if (length(l2_methods) == 0) stop("No usable models remain.")

  model_list <- list()
  for (m in l2_methods) {
    message(sprintf("Training L2 candidate: %s (metric=%s)", m, caret_metric))
    fit <- try(
      caret::train(
        x = l2_train_df,
        y = l2_target,
        method = m,
        trControl = ctrl,
        metric = caret_metric,
        tuneLength = 5
      ), silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      warning(sprintf("Failed to train '%s': %s", m, as.character(fit)))
    } else {
      model_list[[m]] <- fit
    }
  }
  if (length(model_list) == 0) stop("No L2 models were successfully trained.")

  if (length(model_list) == 1) {
    best_name <- names(model_list)[1]
    best_fit  <- model_list[[1]]
    best_val  <- suppressWarnings(max(best_fit$results[[caret_metric]], na.rm = TRUE))
    message(sprintf("Only one model trained. Selected '%s' (CV %s=%.4f).",
                    best_name, caret_metric, best_val))
    resamp <- NULL
  } else {
    resamp <- caret::resamples(model_list)
    vals <- sapply(model_list, function(f) suppressWarnings(max(f$results[[caret_metric]], na.rm = TRUE)))
    best_name <- names(which.max(vals))
    best_fit  <- model_list[[best_name]]
    message(sprintf("Best model: %s (CV %s=%.4f).", best_name, caret_metric, max(vals, na.rm = TRUE)))
  }

  list(
    best_model       = best_fit,
    best_model_name  = best_name,
    best_metric_name = caret_metric,
    model_comparison = resamp,
    l2_train         = data.frame(l2_train_df, .target = l2_target),
    l1_signatures    = l1_signatures
  )
}



#' Derive gene-level importance from a trained meta learner
#'
#' Combines L2 signature importances with L1 signature weights to obtain
#' per-gene contribution scores. When the selected meta learner is a logistic
#' regression, signed contributions indicate whether higher gene expression
#' increases (`> 0`) or decreases (`< 0`) the odds of the positive class.
#' For non-linear models the magnitude is still meaningful, but signs should be
#' interpreted cautiously because the L2 model can capture non-linear effects.
#'
#' @param meta_result A result list returned by [TML6()].
#' @param normalize If `TRUE`, scale contributions so that the maximum absolute
#'   per-signature contribution is 1.
#' @return A list with elements:
#'   \describe{
#'     \item{signature_importance}{Named numeric vector of L1 signature importances.}
#'     \item{gene_importance}{Data frame with columns `signature`, `gene`,
#'           `contribution` (signed), and `abs_contribution`.}
#'     \item{gene_summary}{Aggregated contributions per gene across signatures.}
#'     \item{positive_class}{Name of the positive class (usually `X2`).}
#'     \item{model_type}{Name of the selected L2 model.}
#'   }
#' @export
compute_meta_gene_importance <- function(meta_result, normalize = TRUE) {
  if (!is.list(meta_result) ||
      !all(c("best_model", "best_model_name", "l1_signatures", "l2_train") %in% names(meta_result))) {
    stop("meta_result must be the object returned by TML6().")
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  model_type <- meta_result$best_model_name
  model <- meta_result$best_model
  sig_names <- setdiff(colnames(meta_result$l2_train), ".target")

  if (length(sig_names) == 0) {
    stop("No signature features found in l2_train.")
  }

  # helper to standardise signature definition into named numeric vector
  standardise_signature <- function(sig) {
    if (is.null(sig)) return(numeric(0))
    if (is.numeric(sig) && !is.null(names(sig))) return(sig)
    if (is.character(sig)) return(structure(rep(1, length(sig)), names = sig))
    if (is.list(sig) && all(c("up", "down") %in% names(sig))) {
      up <- sig$up %||% character()
      down <- sig$down %||% character()
      nm <- c(up, down)
      wt <- c(rep(1, length(up)), rep(-1, length(down)))
      names(wt) <- nm
      return(wt)
    }
    if (is.list(sig) && all(c("genes", "weights") %in% names(sig))) {
      g <- sig$genes
      w <- sig$weights
      if (length(g) != length(w)) stop("Signature genes/weights length mismatch.")
      names(w) <- g
      return(w)
    }
    if (is.data.frame(sig)) {
      cn <- tolower(colnames(sig))
      gene_col <- which(cn %in% c("gene","genes","feature","features","symbol","id"))[1]
      weight_col <- which(cn %in% c("weight","weights","w","score","scores","coef","coefs"))[1]
      if (!is.na(gene_col) && !is.na(weight_col)) {
        g <- as.character(sig[[gene_col]])
        w <- as.numeric(sig[[weight_col]])
        ok <- !is.na(g) & !is.na(w)
        w <- w[ok]; g <- g[ok]
        names(w) <- g
        return(w)
      }
    }
    stop("Unsupported signature format when standardising weights.")
  }

  signature_weights <- lapply(meta_result$l1_signatures[sig_names], standardise_signature)

  # obtain signature-level importance
  signature_importance <- NULL
  if (identical(model_type, "glm") && inherits(model$finalModel, "glm")) {
    coefs <- stats::coef(model$finalModel)
    coefs <- coefs[names(coefs) != "(Intercept)"]
    signature_importance <- coefs[sig_names]
  } else {
    vi <- try(caret::varImp(model, scale = FALSE), silent = TRUE)
    if (!inherits(vi, "try-error")) {
      importance_df <- vi$importance
      if ("Overall" %in% colnames(importance_df)) {
        signature_importance <- importance_df[sig_names, "Overall"]
      } else if (ncol(importance_df) >= 1) {
        signature_importance <- importance_df[sig_names, 1]
      }
    }
  }

  if (is.null(signature_importance)) {
    stop("Could not derive signature importance for model type: ", model_type)
  }

  if (normalize) {
    signature_importance <- signature_importance / (max(abs(signature_importance), na.rm = TRUE) %||% 1)
  }

  gene_tables <- lapply(sig_names, function(sig) {
    weights <- signature_weights[[sig]]
    if (length(weights) == 0) return(NULL)
    contrib <- signature_importance[[sig]] * weights
    if (normalize && any(is.finite(contrib))) {
      denom <- max(abs(contrib), na.rm = TRUE)
      if (is.finite(denom) && denom > 0) contrib <- contrib / denom
    }
    data.frame(
      signature = sig,
      gene = names(contrib),
      contribution = as.numeric(contrib),
      abs_contribution = abs(as.numeric(contrib)),
      stringsAsFactors = FALSE
    )
  })

  gene_importance <- do.call(rbind, gene_tables)
  if (is.null(gene_importance) || nrow(gene_importance) == 0) {
    stop("No gene contributions could be computed.")
  }

  gene_summary <- stats::aggregate(contribution ~ gene, data = gene_importance, sum)
  gene_summary$abs_contribution <- abs(gene_summary$contribution)

  positive_class <- tail(levels(meta_result$l2_train$.target), 1)

  list(
    signature_importance = signature_importance,
    gene_importance = gene_importance,
    gene_summary = gene_summary[order(-gene_summary$abs_contribution), ],
    positive_class = positive_class,
    model_type = model_type
  )
}



# LMM analysis suite -----
# PART 1: 데이터 준비 및 설정

#' 분석 설정 객체 생성
#'
#' 분석 파이프라인 전반에서 사용될 메타데이터의 컬럼명(변수명)을
#' 표준화된 리스트 객체로 생성합니다.
#'
#' @param patient 환자 ID 컬럼명. (기본값: "patient_id")
#' @param drug 약물 정보 컬럼명. (기본값: "drug")
#' @param timepoint 시간점(e.g., "pre", "post") 컬럼명. (기본값: "timepoint")
#' @param ck CK status (e.g., "CK+", "CK-") 컬럼명. (기본값: "ck_status")
#' @param response 반응 여부(e.g., "R", "NR") 컬럼명. (기본값: "response")
#' @param aoi AOI (Area of Interest) ID 컬럼명. (기본값: "aoi_id")
#'
#' @return 컬럼명 설정이 담긴 명명된 리스트(named list).
#' @export
#' @examples
#' config <- create_analysis_config(patient = "Subject", drug = "Treatment")
create_analysis_config <- function(
  patient = "patient_id",
  drug = "drug", 
  timepoint = "timepoint",
  ck = "ck_status",
  response = "response",
  aoi = "aoi_id"
) {
  list(
    patient = patient,
    drug = drug,
    timepoint = timepoint,
    ck = ck,
    response = response,
    aoi = aoi
  )
}

#' GeoMx 데이터 준비 및 Seurat 객체 생성
#'
#' Excel/CSV 파일 또는 R 객체(matrix, data.frame)로부터 count 및 metadata를 읽어
#' 분석을 위한 Seurat 객체를 생성합니다.
#'
#' @param count_file Count matrix 파일 경로 (Excel 또는 CSV).
#' @param metadata_file Metadata 파일 경로 (Excel 또는 CSV).
#' @param count_matrix (선택 사항) 파일 대신 사용할 count matrix 객체.
#' @param metadata (선택 사항) 파일 대신 사용할 metadata data.frame 객체.
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param q3_normalize Q3 정규화가 이미 수행되었는지 여부.
#'        `TRUE` (기본값)이면 별도 정규화를 생략하고, `FALSE`이면 LogNormalize를 수행합니다.
#'
#' @return 전처리된 Seurat 객체. `condition_full`과 `condition_simple` 메타데이터가 추가됩니다.
#' @importFrom Seurat CreateSeuratObject NormalizeData
#' @importFrom openxlsx read.xlsx
#' @export
prepare_geomx_data <- function(count_file = NULL, 
                             metadata_file = NULL,
                             count_matrix = NULL,
                             metadata = NULL,
                             config = create_analysis_config(),
                             q3_normalize = TRUE) {
  
  # 파일에서 읽기 또는 직접 입력
  if (!is.null(count_file)) {
    if (grepl("\\.xlsx$", count_file)) {
      count_matrix <- openxlsx::read.xlsx(count_file, sheet = 1, rowNames = TRUE)
    } else {
      count_matrix <- read.csv(count_file, row.names = 1)
    }
  }
  
  if (!is.null(metadata_file)) {
    if (grepl("\\.xlsx$", metadata_file)) {
      metadata <- openxlsx::read.xlsx(metadata_file, sheet = 2)
    } else {
      metadata <- read.csv(metadata_file)
    }
  }
  
  # AOI 이름 매칭
  rownames(metadata) <- metadata[[config$aoi]]
  
  # Seurat 객체 생성
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(count_matrix),
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
  
  # Q3 정규화는 이미 되어있다고 가정
  if (!q3_normalize) {
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  }
  
  # 조합 변수 생성 (빠른 스크리닝용)
  seurat_obj$condition_full <- paste0(
    seurat_obj@meta.data[[config$timepoint]], "_",
    seurat_obj@meta.data[[config$drug]], "_", 
    seurat_obj@meta.data[[config$response]]
  )
  
  seurat_obj$condition_simple <- paste0(
    seurat_obj@meta.data[[config$timepoint]], "_",
    seurat_obj@meta.data[[config$drug]]
  )
  
  return(seurat_obj)
}

# PART 2: 빠른 스크리닝 (Wilcoxon)

#' 빠른 유전자 스크리닝 (Wilcoxon Test)
#'
#' `FindAllMarkers` (Wilcoxon rank sum test)를 사용하여 그룹 간 차등 발현 유전자를
#' 신속하게 스크리닝합니다.
#'
#' @param seurat_obj Seurat 객체.
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param comparison_type 비교할 그룹 기준. "pre_post", "drug", "response", "combined" 중 선택.
#' @param ck_subset (선택 사항) "CK+" 또는 "CK-"를 지정하여 특정 CK status로 필터링.
#' @param min_pct `FindAllMarkers`의 `min.pct` 파라미터. (기본값: 0.1)
#' @param logfc_threshold `FindAllMarkers`의 `logfc.threshold` 파라미터. (기본값: 0.25)
#' @param top_n 반환할 상위 유전자 개수. (기본값: 1000)
#' @param use_adj `TRUE`이면 `p_val_adj` 기준, `FALSE` (기본값)이면 `p_val` 기준으로 정렬.
#'
#' @return 리스트:
#' \item{markers}{`FindAllMarkers` 결과 (data.frame)}
#' \item{top_genes}{상위 유전자 벡터}
#' \item{seurat_obj}{필터링이 적용된 Seurat 객체 (if `ck_subset` is used)}
#' @importFrom Seurat Idents FindAllMarkers
#' @importFrom dplyr arrange head pull unique
#' @importFrom rlang sym
#' @export
quick_screen_genes <- function(seurat_obj,
                               config = create_analysis_config(),
                               comparison_type = "pre_post",
                               ck_subset = NULL,
                               min_pct = 0.1,
                               logfc_threshold = 0.25,
                               top_n = 1000,
                               use_adj = FALSE) {
  
  # CK 서브셋 필터링
  if (!is.null(ck_subset)) {
    seurat_obj <- subset(seurat_obj, 
                         subset = !!sym(config$ck) == ck_subset)
  }
  
  # 비교 그룹 설정
  if (comparison_type == "pre_post") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$timepoint]]
  } else if (comparison_type == "drug") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$drug]]
  } else if (comparison_type == "response") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$response]]
  } else if (comparison_type == "combined") {
    Idents(seurat_obj) <- seurat_obj$condition_full
  }
  
  # FindAllMarkers 실행
  markers <- FindAllMarkers(
    seurat_obj,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = FALSE,
    return.thresh = 1  # 모든 유전자 반환
  )
  
  # p-value 선택 및 상위 유전자 선택
  if (use_adj) {
    top_genes <- markers %>%
      arrange(p_val_adj) %>%
      head(top_n) %>%
      pull(gene) %>%
      unique()
  } else {
    top_genes <- markers %>%
      arrange(p_val) %>%
      head(top_n) %>%
      pull(gene) %>%
      unique()
  }
  
  return(list(
    markers = markers,
    top_genes = top_genes,
    seurat_obj = seurat_obj
  ))
}

# PART 3: Linear Mixed Model 분석

#' 단일 유전자에 대한 LMM (Linear Mixed Model) 적합
#'
#' (내부 헬퍼 함수) `run_lmm_multiple_genes` 내에서 단일 유전자에 대해
#' `lmer`를 사용하여 LMM을 적합합니다.
#'
#' @param gene_expr 단일 유전자의 발현값 벡터 (numeric).
#' @param metadata Seurat 객체의 `meta.data` 슬롯 (data.frame).
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param formula_components (선택 사항) LMM 공식을 동적으로 생성하기 위한 리스트.
#'        (e.g., `list(fixed = c("drug", "timepoint"), random = "(1|patient)")`)
#' @param use_config_names `formula_components` 사용 시, "patient"와 같은 일반 용어를
#'        `config$patient` (`"patient_id"`)로 자동 매핑할지 여부. (기본값: `TRUE`)
#'
#' @return 모델 적합 결과 리스트.
#' \item{model}{`lmerMod` 객체}
#' \item{effects}{계수 요약 (data.frame)}
#' \item{anova}{ANOVA 테이블}
#' \item{converged}{모델 수렴 여부 (logical)}
#' \item{formula}{모델에 사용된 공식 (string)}
#' \item{all_drugs}{데이터에 있는 모든 약물}
#' \item{ref_drug}{모델의 참조(reference) 약물}
#' @importFrom lmerTest anova lmer
#' @importFrom stats as.formula relevel
fit_lmm_single_gene <- function(gene_expr,
                                metadata,
                                config = create_analysis_config(),
                                formula_components = NULL,
                                use_config_names = TRUE) {
  
  # 데이터프레임 생성
  df <- cbind(
    data.frame(Expression = gene_expr),
    metadata
  )
  
  # Formula 생성
  if (is.null(formula_components)) {
    # 기본값: config 변수명 사용
    fixed_effects <- c(config$drug, config$timepoint, config$response)
    interactions <- c(
      paste(config$drug, config$timepoint, sep = ":"),
      paste(config$drug, config$response, sep = ":"),
      paste(config$timepoint, config$response, sep = ":"),
      paste(config$drug, config$timepoint, config$response, sep = ":")
    )
    random <- paste0("(1|", config$patient, ")")
    
    formula_str <- paste0(
      "Expression ~ ",
      paste(c(fixed_effects, interactions), collapse = " + "),
      " + ", random
    )
  } else {
    # formula_components 제공 시
    if (use_config_names) {
      # config 매핑 사용
      mapping <- list(
        "patient" = config$patient,
        "drug" = config$drug,
        "timepoint" = config$timepoint,
        "response" = config$response,
        "ck" = config$ck
      )
      
      # fixed effects 매핑
      fixed_mapped <- sapply(formula_components$fixed, function(x) {
        if (x %in% names(mapping)) mapping[[x]] else x
      })
      
      # interactions 매핑
      if (!is.null(formula_components$interactions)) {
        interactions_mapped <- sapply(formula_components$interactions, function(x) {
          parts <- strsplit(x, ":")[[1]]
          mapped_parts <- sapply(parts, function(p) {
            if (p %in% names(mapping)) mapping[[p]] else p
          })
          paste(mapped_parts, collapse = ":")
        })
      } else {
        interactions_mapped <- NULL
      }
      
      # random effects 매핑
      random_mapped <- formula_components$random
      for (key in names(mapping)) {
        random_mapped <- gsub(key, mapping[[key]], random_mapped)
      }
      
      formula_str <- paste0(
        "Expression ~ ",
        paste(c(fixed_mapped, interactions_mapped), collapse = " + "),
        " + ", random_mapped
      )
    } else {
      # 직접 사용 (매핑 없이)
      fixed_effects <- formula_components$fixed
      interactions <- formula_components$interactions
      formula_str <- paste0(
        "Expression ~ ",
        paste(c(fixed_effects, interactions), collapse = " + "),
        " + ", formula_components$random
      )
    }
  }
  
  # 모델 적합
  tryCatch({
    # Factor 레벨 재정렬 (reference level 설정)
    if (config$drug %in% names(df)) {
      df[[config$drug]] <- relevel(as.factor(df[[config$drug]]), ref = levels(as.factor(df[[config$drug]]))[1])
    }
    if (config$response %in% names(df)) {
      df[[config$response]] <- relevel(as.factor(df[[config$response]]), ref = levels(as.factor(df[[config$response]]))[1])
    }
    
    model <- lmer(as.formula(formula_str), data = df, REML = FALSE)
    
    # 결과 추출
    coef_summary <- summary(model)$coefficients
    anova_result <- anova(model)
    
    # 모든 약물 대비 추출 (reference level 포함)
    all_drugs <- unique(df[[config$drug]])
    ref_drug <- levels(as.factor(df[[config$drug]]))[1]
    
    # 주요 효과 계산
    effects <- data.frame(
      term = rownames(coef_summary),
      estimate = coef_summary[, "Estimate"],
      std_error = coef_summary[, "Std. Error"],
      t_value = coef_summary[, "t value"],
      p_value = coef_summary[, "Pr(>|t|)"],
      ref_drug = ref_drug
    )
    
    return(list(
      model = model,
      effects = effects,
      anova = anova_result,
      converged = TRUE,
      formula = formula_str,
      all_drugs = all_drugs,
      ref_drug = ref_drug
    ))
    
  }, error = function(e) {
    return(list(
      converged = FALSE,
      error = e$message,
      formula = formula_str
    ))
  })
}

#' 다중 유전자에 대한 LMM 병렬 수행
#'
#' 지정된 유전자 목록에 대해 `fit_lmm_single_gene` 함수를 병렬로 실행하여
#' LMM 분석을 수행합니다.
#'
#' @param seurat_obj Seurat 객체.
#' @param genes (선택 사항) 분석할 유전자 이름 벡터.
#'        `NULL` (기본값)인 경우, 경고와 함께 첫 100개 유전자를 사용합니다.
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param formula_components (선택 사항) `fit_lmm_single_gene`로 전달될 LMM 공식 리스트.
#' @param use_config_names `formula_components` 사용 시, `config` 매핑 사용 여부. (기본값: `TRUE`)
#' @param n_cores 병렬 처리에 사용할 CPU 코어 수. (기본값: 16)
#' @param verbose 진행 상황 메시지를 출력할지 여부. (기본값: `TRUE`)
#'
#' @return 리스트:
#' \item{raw_results}{각 유전자에 대한 `fit_lmm_single_gene`의 원시 결과 리스트}
#' \item{summary}{`summarize_lmm_results`로 요약된 전체 결과 (data.frame)}
#' \item{converged_genes}{수렴에 성공한 유전자 수}
#' \item{total_genes}{시도한 총 유전자 수}
#' @importFrom Seurat GetAssayData
#' @importFrom parallel makeCluster clusterEvalQ clusterExport parLapply stopCluster
#' @export
run_lmm_multiple_genes <- function(seurat_obj,
                                 genes = NULL,
                                 config = create_analysis_config(),
                                 formula_components = NULL,
                                 use_config_names = TRUE,
                                 n_cores = 16,
                                 verbose = TRUE) {
  
  if (is.null(genes)) {
    genes <- rownames(seurat_obj)[1:100]  # 기본값: 상위 100개
    warning("No genes specified. Using first 100 genes.")
  }
  
  # 발현 매트릭스와 메타데이터 추출
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  metadata <- seurat_obj@meta.data
  
  # 진행상황 표시
  if (verbose) {
    message(sprintf("Running LMM for %d genes using %d cores...", 
                    length(genes), n_cores))
  }
  
  # 병렬 처리
  if (n_cores > 1) {
    cl <- makeCluster(n_cores, type="PSOCK")
    clusterEvalQ(cl, {
      library(lme4)
      library(lmerTest)
    })
    clusterExport(cl, c("fit_lmm_single_gene", "config", "metadata", 
                        "expr_matrix", "formula_components", "use_config_names"),
                  envir = environment())
    
    results <- parLapply(cl, genes, function(gene) {
      if (gene %in% rownames(expr_matrix)) {
        result <- fit_lmm_single_gene(
          gene_expr = as.numeric(expr_matrix[gene, ]),
          metadata = metadata,
          config = config,
          formula_components = formula_components,
          use_config_names = use_config_names
        )
        result$gene <- gene
        return(result)
      } else {
        return(list(gene = gene, converged = FALSE, 
                    error = "Gene not found"))
      }
    })
    
    stopCluster(cl)
  } else {
    # 순차 처리
    results <- lapply(genes, function(gene) {
      if (verbose && which(genes == gene) %% 100 == 0) {
        message(sprintf("  Processing gene %d/%d", 
                        which(genes == gene), length(genes)))
      }
      
      if (gene %in% rownames(expr_matrix)) {
        result <- fit_lmm_single_gene(
          gene_expr = as.numeric(expr_matrix[gene, ]),
          metadata = metadata,
          config = config,
          formula_components = formula_components
        )
        result$gene <- gene
        return(result)
      } else {
        return(list(gene = gene, converged = FALSE, 
                    error = "Gene not found"))
      }
    })
  }
  
  names(results) <- genes
  
  # 결과 요약
  summary_df <- summarize_lmm_results(results, config)
  
  return(list(
    raw_results = results,
    summary = summary_df,
    converged_genes = sum(sapply(results, function(x) x$converged)),
    total_genes = length(genes)
  ))
}

#' LMM 결과 요약 및 p-value 보정
#'
#' (내부 헬퍼 함수) `run_lmm_multiple_genes`에서 생성된 원시 결과 리스트를
#' 하나의 요약 data.frame으로 통합하고, 용어(term)별로 p-value를 보정(BH)합니다.
#'
#' @param lmm_results `run_lmm_multiple_genes`의 `$raw_results` 객체 (리스트).
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#'
#' @return 모든 유전자와 모든 효과(term)를 포함하는 요약 data.frame.
#'         `p_adj` 및 `significant` 컬럼이 추가됩니다.
#' @importFrom dplyr bind_rows group_by mutate ungroup
#' @importFrom stats p.adjust
summarize_lmm_results <- function(lmm_results, config) {
  # 수렴한 모델만 처리
  converged <- lmm_results[sapply(lmm_results, function(x) x$converged)]
  
  if (length(converged) == 0) {
    warning("No models converged!")
    return(NULL)
  }
  
  # 모든 효과 수집
  all_effects <- do.call(rbind, lapply(names(converged), function(gene) {
    effects <- converged[[gene]]$effects
    effects$gene <- gene
    return(effects)
  }))
  
  # p-value 보정
  all_effects <- all_effects %>%
    group_by(term) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05
    ) %>%
    ungroup()
  
  # Reference drug 정보 추가
  if ("ref_drug" %in% names(converged[[1]]$effects)) {
    ref_drugs <- unique(sapply(converged, function(x) x$ref_drug))
    attr(all_effects, "ref_drug") <- ref_drugs[1]
  }
  
  return(all_effects)
}


# validation tools (not using) -----

#' Validate Seurat Object
#'
#' Checks if the input is a valid Seurat object with expected components.
#'
#' @param obj Object to validate
#' @param assay Optional assay name to check for
#' @param reduction Optional reduction name to check for
#' @param min_cells Minimum number of cells required
#' @param min_features Minimum number of features required
#'
#' @return TRUE if valid, stops with error message otherwise
#' @keywords internal
#' @export
validate_seurat <- function(obj, 
                            assay = NULL, 
                            reduction = NULL,
                            min_cells = 0,
                            min_features = 0) {
  
  # Check if it's a Seurat object
  if (!inherits(obj, "Seurat")) {
    stop("Input must be a Seurat object. Current class: ", 
         paste(class(obj), collapse = ", "))
  }
  
  # Check cell count
  n_cells <- ncol(obj)
  if (n_cells < min_cells) {
    stop("Seurat object has ", n_cells, " cells, but ", min_cells, " required")
  }
  
  # Check feature count
  n_features <- nrow(obj)
  if (n_features < min_features) {
    stop("Seurat object has ", n_features, " features, but ", min_features, " required")
  }
  
  # Check assay if specified
  if (!is.null(assay)) {
    available_assays <- Seurat::Assays(obj)
    if (!assay %in% available_assays) {
      stop("Assay '", assay, "' not found. Available assays: ", 
           paste(available_assays, collapse = ", "))
    }
  }
  
  # Check reduction if specified
  if (!is.null(reduction)) {
    available_reductions <- names(obj@reductions)
    if (!reduction %in% available_reductions) {
      # Try case-insensitive match
      matched_reduction <- available_reductions[tolower(available_reductions) == tolower(reduction)]
      if (length(matched_reduction) == 1) {
        message("Note: Using case-insensitive match: '", matched_reduction, "' for '", reduction, "'")
        return(TRUE)
      }
      stop("Reduction '", reduction, "' not found. Available reductions: ", 
           paste(available_reductions, collapse = ", "))
    }
  }
  
  return(TRUE)
}

#' Validate Metadata Column
#'
#' Checks if a metadata column exists and optionally validates its type.
#'
#' @param obj Seurat object or data frame
#' @param column_name Column name to validate
#' @param required_type Optional required type ("numeric", "factor", "character")
#' @param allow_na Whether NA values are allowed
#'
#' @return TRUE if valid, stops with error message otherwise
#' @keywords internal
#' @export
validate_metadata_column <- function(obj, 
                                     column_name, 
                                     required_type = NULL,
                                     allow_na = TRUE) {
  
  # Get metadata
  if (inherits(obj, "Seurat")) {
    metadata <- obj@meta.data
  } else if (is.data.frame(obj)) {
    metadata <- obj
  } else {
    stop("Input must be a Seurat object or data frame")
  }
  
  # Check if column exists
  if (!column_name %in% colnames(metadata)) {
    stop("Column '", column_name, "' not found in metadata. Available columns: ",
         paste(head(colnames(metadata), 10), collapse = ", "),
         if (ncol(metadata) > 10) "..." else "")
  }
  
  column_data <- metadata[[column_name]]
  
  # Check for NA values if not allowed
  if (!allow_na && any(is.na(column_data))) {
    n_na <- sum(is.na(column_data))
    stop("Column '", column_name, "' contains ", n_na, " NA values, but NA not allowed")
  }
  
  # Check type if specified
  if (!is.null(required_type)) {
    is_correct_type <- switch(
      required_type,
      "numeric" = is.numeric(column_data),
      "factor" = is.factor(column_data),
      "character" = is.character(column_data),
      stop("Unknown required_type: ", required_type)
    )
    
    if (!is_correct_type) {
      stop("Column '", column_name, "' must be of type '", required_type, 
           "', but is: ", paste(class(column_data), collapse = ", "))
    }
  }
  
  return(TRUE)
}

#' Validate Gene List
#'
#' Checks if genes exist in the object and optionally filters to valid genes.
#'
#' @param obj Seurat object or character vector of available genes
#' @param genes Character vector of genes to validate
#' @param min_present Minimum number of genes that must be present (default: all)
#' @param assay Assay to check genes in (for Seurat objects)
#' @param warn_missing Whether to warn about missing genes
#'
#' @return Character vector of valid genes present in the object
#' @keywords internal
#' @export
validate_genes <- function(obj, 
                          genes, 
                          min_present = NULL,
                          assay = NULL,
                          warn_missing = TRUE) {
  
  if (!is.character(genes) || length(genes) == 0) {
    stop("genes must be a non-empty character vector")
  }
  
  # Get available genes
  if (inherits(obj, "Seurat")) {
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(obj)
    }
    available_genes <- rownames(obj[[assay]])
  } else if (is.character(obj)) {
    available_genes <- obj
  } else {
    stop("obj must be a Seurat object or character vector of gene names")
  }
  
  # Find present and missing genes
  genes_present <- intersect(genes, available_genes)
  genes_missing <- setdiff(genes, available_genes)
  
  # Set default minimum
  if (is.null(min_present)) {
    min_present <- length(genes)
  }
  
  # Check if enough genes are present
  if (length(genes_present) < min_present) {
    stop("Only ", length(genes_present), " of ", length(genes), 
         " genes found, but ", min_present, " required. ",
         "First missing genes: ", 
         paste(head(genes_missing, 5), collapse = ", "),
         if (length(genes_missing) > 5) "..." else "")
  }
  
  # Warn about missing genes if requested
  if (warn_missing && length(genes_missing) > 0) {
    warning("Genes not found (", length(genes_missing), "/", length(genes), "): ",
            paste(head(genes_missing, 10), collapse = ", "),
            if (length(genes_missing) > 10) "..." else "")
  }
  
  return(genes_present)
}

#' Validate Numeric Range
#'
#' Checks if a numeric parameter is within acceptable range.
#'
#' @param value Numeric value to validate
#' @param param_name Name of parameter (for error messages)
#' @param min Minimum allowed value (inclusive)
#' @param max Maximum allowed value (inclusive)
#' @param allow_na Whether NA is allowed
#'
#' @return TRUE if valid, stops with error message otherwise
#' @keywords internal
#' @export
validate_numeric_range <- function(value, 
                                   param_name, 
                                   min = -Inf, 
                                   max = Inf,
                                   allow_na = FALSE) {
  
  # Check for NA
  if (is.na(value)) {
    if (allow_na) {
      return(TRUE)
    } else {
      stop(param_name, " cannot be NA")
    }
  }
  
  # Check if numeric
  if (!is.numeric(value) || length(value) != 1) {
    stop(param_name, " must be a single numeric value")
  }
  
  # Check range
  if (value < min || value > max) {
    stop(param_name, " must be between ", min, " and ", max, 
         ", but is: ", value)
  }
  
  return(TRUE)
}

#' Validate Choice
#'
#' Validates that a value is one of allowed choices (like match.arg but more informative).
#'
#' @param value Value to validate
#' @param param_name Name of parameter (for error messages)
#' @param choices Vector of allowed values
#' @param multiple Whether multiple choices are allowed
#'
#' @return The validated value (or values if multiple=TRUE)
#' @keywords internal
#' @export
validate_choice <- function(value, param_name, choices, multiple = FALSE) {
  
  if (is.null(value)) {
    stop(param_name, " cannot be NULL")
  }
  
  if (multiple) {
    if (!all(value %in% choices)) {
      invalid <- setdiff(value, choices)
      stop(param_name, " contains invalid choices: ", 
           paste(invalid, collapse = ", "),
           ". Allowed choices: ", 
           paste(choices, collapse = ", "))
    }
  } else {
    if (length(value) != 1) {
      stop(param_name, " must be a single value, not ", length(value))
    }
    if (!value %in% choices) {
      stop(param_name, " must be one of: ", 
           paste(choices, collapse = ", "),
           ". Got: ", value)
    }
  }
  
  return(value)
}

#' Validate File Path
#'
#' Checks if a file path exists and optionally validates extension.
#'
#' @param path File path to validate
#' @param must_exist Whether file must already exist
#' @param extensions Optional vector of allowed extensions (e.g., c("csv", "txt"))
#' @param type Type of path ("file" or "directory")
#'
#' @return Normalized path if valid, stops with error otherwise
#' @keywords internal
#' @export
validate_path <- function(path, 
                         must_exist = TRUE, 
                         extensions = NULL,
                         type = c("file", "directory")) {
  
  type <- match.arg(type)
  
  if (!is.character(path) || length(path) != 1) {
    stop("path must be a single character string")
  }
  
  # Check existence if required
  if (must_exist) {
    if (type == "file") {
      if (!file.exists(path)) {
        stop("File not found: ", path)
      }
      if (dir.exists(path)) {
        stop("Expected a file but got a directory: ", path)
      }
    } else if (type == "directory") {
      if (!dir.exists(path)) {
        stop("Directory not found: ", path)
      }
    }
  }
  
  # Check extension if specified
  if (!is.null(extensions) && type == "file") {
    file_ext <- tools::file_ext(path)
    if (!file_ext %in% extensions) {
      stop("File extension must be one of: ", 
           paste(extensions, collapse = ", "),
           ". Got: ", file_ext)
    }
  }
  
  # Return normalized path
  return(normalizePath(path, mustWork = must_exist))
}

#' Create Informative Error Message
#'
#' Helper to create consistent, informative error messages.
#'
#' @param context Context where error occurred (e.g., function name)
#' @param message Main error message
#' @param suggestion Optional suggestion for fixing the error
#'
#' @return Formatted error message
#' @keywords internal
#' @export
create_error_message <- function(context, message, suggestion = NULL) {
  msg <- paste0("[", context, "] ", message)
  if (!is.null(suggestion)) {
    msg <- paste0(msg, "\nSuggestion: ", suggestion)
  }
  return(msg)
}

#' Check Package Dependencies
#'
#' Checks if required packages are installed and optionally loads them.
#'
#' @param packages Character vector of package names
#' @param load Whether to load the packages (default: FALSE)
#'
#' @return TRUE if all packages available, stops with error otherwise
#' @keywords internal
#' @export
check_packages <- function(packages, load = FALSE) {
  
  missing <- character(0)
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    } else if (load) {
      library(pkg, character.only = TRUE)
    }
  }
  
  if (length(missing) > 0) {
    stop("Required packages not installed: ", 
         paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c('", 
         paste(missing, collapse = "', '"), "'))")
  }
  
  return(TRUE)
}


# plots -----

scatter_smooth_colored3 <- function(object,
                                   feature,
                                   group.by   = "sample_no",
                                   x_var      = "nih_change",
                                   transpose  = FALSE,
                                   color_by   = NULL,
                                   split.by   = NULL,
                                   palette    = NULL,
                                   transparency    = TRUE,
                                   transparency_desc = FALSE,
                                   fitted_line = c("linear", "loess", "lasso", NULL)) {
    fitted_line <- match.arg(fitted_line)
    stopifnot(is.character(feature), length(feature) == 1)
    
    # 1. Build per‑cell tibble ------------------------------------------------
    
    if (inherits(object, "Seurat")) {
        expr_vec <- Seurat::FetchData(object, vars = feature)[, 1]
        meta_df  <- tibble::as_tibble(object@meta.data)
        cell_df  <- dplyr::mutate(meta_df, !!feature := expr_vec)
    } else {
        cell_df <- tibble::as_tibble(object)
        if (!feature %in% names(cell_df))
            stop("In data.frame mode the column '", feature, "' must exist.")
    }
    
    for (col in c(group.by, x_var, color_by, split.by)) {
        if (!is.null(col) && !col %in% names(cell_df))
            stop("Column '", col, "' not found in data.")
    }
    
    if (!is.null(split.by)) {
        if (is.numeric(cell_df[[split.by]])) {
            stop("`split.by` column '", split.by, "' must be categorical (character or factor), not numeric.")
        }
        cell_df[[split.by]] <- as.factor(cell_df[[split.by]])
    }
    
    
    # 2. Aggregate by group ---------------------------------------------------
    
    is_color_numeric <- if (!is.null(color_by)) {
        is.numeric(cell_df[[color_by]])
    } else {
        FALSE
    }
    
    agg_df <- cell_df %>%
        dplyr::group_by(.data[[group.by]]) %>%
        dplyr::summarise(
            avg_expr = mean(.data[[feature]], na.rm = TRUE),
            x_val    = mean(.data[[x_var]], na.rm = TRUE),
            colour   = if (!is.null(color_by)) {
                if (is_color_numeric) {
                    mean(.data[[color_by]], na.rm = TRUE)
                } else {
                    dplyr::first(.data[[color_by]], na_rm = TRUE)
                }
            } else {
                NA
            },
            split_col = if (!is.null(split.by)) {
                dplyr::first(.data[[split.by]], na_rm = TRUE)
            } else {
                NA
            },
            .groups  = "drop"
        )
    
    
    # 3. Aesthetics -----------------------------------------------------------
    
    x_col <- if (transpose) "avg_expr" else "x_val"
    y_col <- if (transpose) "x_val"   else "avg_expr"
    
    p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
    
    if (!is.null(split.by)) {
        if (!is.null(color_by)) {
            warning("`color_by` argument is ignored when `split.by` is provided.")
        }
        p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data[["split_col"]]), size = 3)
        
        pal <- palette %||% RColorBrewer::brewer.pal(max(3, length(unique(agg_df$split_col))), "Set1")
        p <- p + ggplot2::scale_colour_manual(values = pal, name = split.by)
        
    } else if (!is.null(color_by)) {
        if (is.numeric(cell_df[[color_by]])) {
            p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour,
                                                      alpha  = colour), size = 3)
            alpha_range <- if (transparency_desc) c(1, 0.2) else c(0.2, 1)
            if (transparency) {
                p <- p + ggplot2::scale_alpha(range = alpha_range, guide = "none")
            } else {
                p <- p + ggplot2::guides(alpha = "none")
            }
            pal <- if (is.null(palette)) viridisLite::viridis(256) else palette
            p <- p + ggplot2::scale_colour_gradientn(colours = pal, name = color_by)
        } else {
            p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour), size = 3)
            pal <- palette %||% RColorBrewer::brewer.pal(max(3, length(unique(agg_df$colour))), "Set1")
            p <- p + ggplot2::scale_colour_manual(values = pal, name = color_by)
        }
    } else {
        p <- p + ggplot2::geom_point(size = 3)
    }
    
    
    # 4. Smoothing line -------------------------------------------------------
    
    if (!is.null(fitted_line)) {
        
        if (!is.null(split.by)) {
            
            # *** 수정된 부분 시작 ***
            # `fitted_line` 인수를 `geom_smooth`가 이해하는 `method` 문자열로 변환
            
            method_val <- NULL # 초기화
            
            if (fitted_line == "linear") {
                method_val <- "lm"
            } else if (fitted_line == "loess") {
                method_val <- "loess"
            } else if (fitted_line == "lasso") {
                warning("ggplot2::geom_smooth does not support 'lasso'. Falling back to 'lm' method.")
                method_val <- "lm"
            }
            # *** 수정된 부분 끝 ***
            
            warning("Custom statistical annotations (p-value, equation) are disabled when `split.by` is used.")
            
            p <- p + ggplot2::geom_smooth(ggplot2::aes(colour = .data[["split_col"]]),
                                          method = method_val, # "linear" 대신 "lm"이 전달됨
                                          se = TRUE,
                                          show.legend = FALSE) 
            
        } else {
            # (원본 로직: split.by가 NULL일 때만 실행)
            if (fitted_line == "linear") {
                p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black")
                fit <- stats::lm(agg_df[[y_col]] ~ agg_df[[x_col]])
                coef <- round(stats::coef(fit), 3)
                pval <- signif(summary(fit)$coefficients[2, 4], 3)
                annot <- paste0("y = ", coef[1], " + ", coef[2], " * x\np = ", pval)
                p <- p + ggplot2::annotate("text", x = min(agg_df[[x_col]], na.rm = TRUE),
                                           y = max(agg_df[[y_col]], na.rm = TRUE),
                                           label = annot, hjust = 0, vjust = 1, size = 4)
            } else if (fitted_line == "loess") {
                p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "black")
            } else if (fitted_line == "lasso") {
                if (!requireNamespace("glmnet", quietly = TRUE)) {
                    warning("glmnet not installed; falling back to linear fit.")
                    p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE)
                } else {
                    xmat <- as.matrix(agg_df[[x_col]])
                    fit  <- glmnet::cv.glmnet(xmat, agg_df[[y_col]], alpha = 1)
                    preds <- as.numeric(glmnet::predict.glmnet(fit$glmnet.fit, newx = xmat,
                                                               s = fit$lambda.min))
                    pred_df <- agg_df %>% dplyr::mutate(pred = preds)
                    p <- p + ggplot2::geom_line(data = pred_df[order(pred_df[[x_col]]), ],
                                                ggplot2::aes(x = .data[[x_col]], y = pred),
                                                colour = "red", linewidth = 1)
                }
            }
        }
    }
    
    
    # 5. Labels & theme -------------------------------------------------------
    
    p <- p + ggplot2::theme_bw() +
        ggplot2::labs(x = if (transpose) paste("Average", feature, "expression") else x_var,
                      y = if (transpose) x_var else paste("Average", feature, "expression"))
    
    p
    return(p)
}

plot_volcano <- function(lmm_summary,
                         x_col="estimate",
                         y_col="p_value",
                         filter_col="term",
                        filter_pattern = NULL,
                        title = "Volcano Plot",
                        effect_threshold = 0.5,
                        p_threshold = 0.05,
                        resize_x=TRUE,
                        resize_y=TRUE) {
  #ver2: generalized
  plot_data <- lmm_summary
  
  if (!is.null(filter_pattern)) {
    plot_data <- plot_data %>%
      filter(grepl(filter_pattern, term))
  }
  if (!is.null(x_col)){plot_data=plot_data%>%mutate(effect_size= !!sym(x_col))}
  if (!is.null(y_col)){plot_data=plot_data%>%mutate(p_value= !!sym(y_col))}

  # plot label definition
  x_label="Effect Size"
  y_label="Significance"
  
  # plot resizing
  if(resize_x){xlim_val=quantile(abs(plot_data$effect_size),0.95,na.rm=TRUE)}
  if(resize_y){ylim_val=quantile(abs(plot_data$effect_size),0.95,na.rm=TRUE)}
  
  # gene significance label arrange
  plot_data <- plot_data %>%
    mutate(
      log_p = -log10(p_value),
      category = case_when(
        abs(effect_size) > effect_threshold & p_value < p_threshold ~ "Significant",
        abs(effect_size) > effect_threshold ~ "Large effect",
        p_value  < p_threshold ~ "Small effect",
        TRUE ~ "Not significant"
      )
    )
  ggplot(plot_data, aes(x = effect_size, y = log_p, color = category)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold), 
               linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Significant" = "red",
                                 "Large effect" = "orange", 
                                 "Small effect" = "blue",
                                 "Not significant" = "gray")) +
    labs(title = title,
         x = x_label,
         y = y_label) +
    theme_bw()+
    coord_cartesian(xlim = c(-xlim_val, xlim_val))
}



