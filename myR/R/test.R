

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

#' Seurat 객체에서 유전자별 CLMM 분석 수행 (v2: 기능 추가)
#'
#' @param sobj Seurat 객체
#' @param ordinal_var 메타데이터의 순서형 반응 변수 (예: "heal")
#' @param patient_col 메타데이터의 환자 ID(random effect) (예: "emrid")
#' @param assay 사용할 Seurat assay (기본값: "RNA")
#' @param layer 사용할 Seurat layer (v4/v3의 'slot') (기본값: "data")
#' @param top_n 분석할 상위 유전자 개수 (기본값: NULL = 모든 유전자)
#' @param sort_by 'top_n'을 사용할 경우 정렬 기준 ("variance" 또는 "mean")
#' @param decreasing 'top_n' 정렬 시 내림차순 (기본값: TRUE)
#' @param ... clmm() 함수에 전달할 추가 인자 (예: control = clmm.control(maxIter = 200))
#'
#' @return 'gene' 및 CLMM 결과 (Estimate, p-value 등) data.frame
#'
run_clmm_analysis <- function(sobj, ordinal_var, patient_col, assay = "RNA", layer = "data", 
                              top_n = NULL, sort_by = "variance", decreasing = TRUE, ...) {
  
  # 0. 필수 패키지 확인
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("'ordinal' 패키지가 필요합니다: install.packages('ordinal')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("'dplyr' 패키지가 필요합니다: install.packages('dplyr')")
  }
  
  # 1. 메타데이터 준비
  meta <- sobj@meta.data
  unique_levels <- sort(unique(meta[[ordinal_var]]))
  meta[[ordinal_var]] <- factor(meta[[ordinal_var]], levels = unique_levels, ordered = TRUE)
  cat("Note: '", ordinal_var, "'를 ordered factor로 변환했습니다: ", paste(unique_levels, collapse = "<"), "\n")
  meta[[patient_col]] <- as.factor(meta[[patient_col]])
  
  # 2. 유전자 발현 데이터 및 유전자 목록 준비
  expr_matrix <- GetAssayData(sobj, assay = assay, slot = layer)
  
  if (!is.null(top_n)) {
    cat(paste("'top_n'이 설정되었습니다.", top_n, "개 유전자를", sort_by, "기준으로 정렬합니다...\n"))
    
    if (sort_by == "variance") {
      metrics <- apply(expr_matrix, 1, var)
    } else if (sort_by == "mean") {
      metrics <- rowMeans(expr_matrix)
    } else {
      stop("'sort_by'는 'variance' 또는 'mean'이어야 합니다.")
    }
    
    ordered_indices <- order(metrics, decreasing = decreasing)
    genes_to_test <- rownames(expr_matrix)[ordered_indices][1:min(top_n, length(metrics))]
    
  } else {
    genes_to_test <- rownames(expr_matrix)
  }
  
  total_genes <- length(genes_to_test)
  cat("총", total_genes, "개의 유전자에 대해 CLMM 모델을 실행합니다...\n")
  
  # 3. clmm 모델 실행 (for loop로 변경: 진행 상황 및 시간 측정)
  clmm_results_list <- list()
  overall_start_time <- Sys.time()
  
  for (i in seq_along(genes_to_test)) {
    gene <- genes_to_test[i]
    
    model_data <- data.frame(
      ord_var = meta[[ordinal_var]],
      pat_var = meta[[patient_col]],
      gene_expr = as.numeric(expr_matrix[gene, ])
    )
    
    # 모델 피팅 (오류 발생 시 NULL 반환)
    model_fit <- tryCatch({
      # '...' 인자를 통해 clmm.control(maxIter=200) 등 전달 가능
      ordinal::clmm(ord_var ~ gene_expr + (1|pat_var), data = model_data, ...)
    }, error = function(e) {
      NULL
    })
    
    # 결과 추출
    if (!is.null(model_fit)) {
      coef_summary <- summary(model_fit)$coefficients
      if ("gene_expr" %in% rownames(coef_summary)) {
        clmm_results_list[[gene]] <- coef_summary["gene_expr", ]
      } else {
        clmm_results_list[[gene]] <- rep(NA, 4)
      }
    } else {
      clmm_results_list[[gene]] <- rep(NA, 4)
    }
    
    # 4. 진행 상황 및 예상 시간 리포트
    if (i == 100) {
      first_100_time <- Sys.time()
      time_per_100 <- as.numeric(difftime(first_100_time, overall_start_time, units = "secs"))
      total_estimated_secs <- (time_per_100 / 100) * total_genes
      total_estimated_mins <- total_estimated_secs / 60
      
      cat(paste0("[진행 상황] 첫 100개 유전자 완료 (", round(time_per_100, 1), "초 소요)\n"))
      cat(paste0("   > 총 예상 소요 시간: 약 ", round(total_estimated_mins, 1), " 분\n"))
      
    } else if (i > 100 && i %% 100 == 0) {
      current_time <- Sys.time()
      elapsed_mins <- as.numeric(difftime(current_time, overall_start_time, units = "mins"))
      cat(paste0("[진행 상황] ", i, " / ", total_genes, " (", round(i/total_genes*100), "%) 완료. (", round(elapsed_mins, 1), "분 경과)\n"))
    }
  } # end for loop
  
  # 5. 결과 data.frame으로 변환
  clmm_results_df <- do.call(rbind, clmm_results_list)
  rownames(clmm_results_df) <- names(clmm_results_list) # genes_to_test 순서
  colnames(clmm_results_df) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  clmm_results_df <- as.data.frame(clmm_results_df)
  clmm_results_df$gene <- rownames(clmm_results_df)
  
  clmm_results_df <- clmm_results_df[!is.na(clmm_results_df$`Pr(>|z|)`), ]
  clmm_results_df <- clmm_results_df %>% 
    dplyr::arrange(`Pr(>|z|)`) %>%
    dplyr::select(gene, dplyr::everything())
  
  cat("CLMM 분석 완료.\n")
  return(clmm_results_df)
}

#' Seurat 객체에서 유전자별 스피어만 상관관계 분석 (Pseudobulk) 수행 (v2: 버그 수정 및 기능 추가)
#'
#' @param sobj Seurat 객체
#' @param ordinal_var 메타데이터의 순서형 변수 (예: "heal")
#' @param patient_col 메타데이터의 환자 ID (예: "emrid")
#' @param assay 사용할 Seurat assay (기본값: "RNA")
#' @param layer 사용할 Seurat layer (v4/v3의 'slot') (기본값: "data")
#' @param top_n 분석할 상위 유전자 개수 (기본값: NULL = 모든 유전자)
#' @param sort_by 'top_n'을 사용할 경우 정렬 기준 ("variance" 또는 "mean")
#' @param decreasing 'top_n' 정렬 시 내림차순 (기본값: TRUE)
#'
#' @return 'gene', 'rho', 'p.value'를 포함하는 data.frame
#'
run_spearman_pseudobulk <- function(sobj, ordinal_var, patient_col, assay = "RNA", layer = "data",
                                    top_n = NULL, sort_by = "variance", decreasing = TRUE) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("'dplyr' 패키지가 필요합니다: install.packages('dplyr')")
  }

  # 1. (오류 수정) 숫자로 시작하는 ID에 'g'를 덧붙인 임시 컬럼 생성
  #    AggregateExpression이 내부적으로 하는 작업을 미리 수행하여 ID 불일치 방지
  temp_group_by_col <- "temp_spearman_group_by"
  sobj@meta.data[[temp_group_by_col]] <- ifelse(
    grepl("^[0-9]", sobj@meta.data[[patient_col]]),
    paste0("g", sobj@meta.data[[patient_col]]),
    as.character(sobj@meta.data[[patient_col]])
  )

  # 2. Pseudobulk 발현 매트릭스 생성
  cat("환자별 Pseudobulk 프로필을 생성합니다 (임시 ID 사용)...\n")
  avg_expr <- AggregateExpression(
    sobj,
    group.by = temp_group_by_col, # 수정된 임시 컬럼 사용
    assays = assay,
    slot = layer, 
    return.seurat = FALSE
  )
  
  pseudobulk_matrix <- avg_expr[[assay]]
  
  # (참고) AggregateExpression의 'Inf' 경고 메시지:
  # 'layer = "data"' (보통 log-normalized)를 사용할 때 발생할 수 있습니다.
  # 이는 'AggregateExpression'이 내부적으로 'exp(x)-1'을 시도할 때 발생할 수 있으나,
  # 'slot = "data"'의 경우 평균(mean)을 계산하므로 결과(평균 발현 값)는 유효합니다.

  # 3. 환자별 메타데이터(ordinal_var) 준비
  meta <- sobj@meta.data
  
  # 임시 ID를 기준으로 유니크한 메타데이터 추출
  patient_meta <- meta %>%
    dplyr::select(dplyr::all_of(c(temp_group_by_col, ordinal_var))) %>%
    dplyr::distinct(!!dplyr::sym(temp_group_by_col), .keep_all = TRUE)
  
  # pseudobulk_matrix의 컬럼 이름(환자 ID) 순서대로 patient_meta 정렬
  pseudobulk_patients <- colnames(pseudobulk_matrix)
  match_indices <- match(pseudobulk_patients, patient_meta[[temp_group_by_col]])
  
  # 오류 검사 (NA가 없는지)
  if (any(is.na(match_indices))) {
      stop(paste("Patient ID 매칭 실패. 임시 ID 생성 후에도 NA가 발생:",
                 paste(pseudobulk_patients[is.na(match_indices)], collapse = ",")))
  }
  
  ordered_patient_meta <- patient_meta[match_indices, ]
  
  # (오류 수정) 원본 코드의 if문 검사를 더 강력하게 수정
  if (!all(pseudobulk_patients == ordered_patient_meta[[temp_group_by_col]])) {
    stop("환자 ID 정렬 오류: Pseudobulk 순서와 메타데이터 순서가 일치하지 않습니다.")
  }
  
  heal_values <- ordered_patient_meta[[ordinal_var]]
  
  # 4. 유전자 필터링 (Pseudobulk Matrix 기준)
  if (!is.null(top_n)) {
    cat(paste("'top_n'이 설정되었습니다.", top_n, "개 유전자를", sort_by, "기준으로 정렬합니다 (Pseudobulk 기준)...\n"))
    
    if (sort_by == "variance") {
      metrics <- apply(pseudobulk_matrix, 1, var)
    } else if (sort_by == "mean") {
      metrics <- rowMeans(pseudobulk_matrix)
    } else {
      stop("'sort_by'는 'variance' 또는 'mean'이어야 합니다.")
    }
    
    ordered_indices <- order(metrics, decreasing = decreasing)
    genes_to_test <- rownames(pseudobulk_matrix)[ordered_indices][1:min(top_n, length(metrics))]
    
    # 분석할 매트릭스 서브셋
    pseudobulk_matrix_subset <- pseudobulk_matrix[genes_to_test, ]
    
  } else {
    genes_to_test <- rownames(pseudobulk_matrix)
    pseudobulk_matrix_subset <- pseudobulk_matrix
  }

  # 5. 스피어만 상관관계 계산 (for loop 사용)
  total_genes <- length(genes_to_test)
  cat("총", total_genes, "개의 유전자에 대해 스피어만 상관관계를 계산합니다...\n")
  
  spearman_results_list <- list()
  overall_start_time <- Sys.time()
  
  for (i in seq_along(genes_to_test)) {
    gene <- genes_to_test[i]
    gene_avg_expr <- pseudobulk_matrix_subset[gene, ]
    
    # 발현에 분산이 있는지 확인
    if (var(gene_avg_expr, na.rm=TRUE) == 0) {
      spearman_results_list[[gene]] <- c(rho = NA, p.value = NA)
      next # 다음 유전자로
    }
    
    cor_test_result <- tryCatch({
      stats::cor.test(gene_avg_expr, heal_values, method = "spearman")
    }, error = function(e) { NULL })
    
    if (!is.null(cor_test_result)) {
      spearman_results_list[[gene]] <- c(
        rho = cor_test_result$estimate,
        p.value = cor_test_result$p.value
      )
    } else {
      spearman_results_list[[gene]] <- c(rho = NA, p.value = NA)
    }
    
    # 진행 상황 리포트
    if (i == 100) {
      first_100_time <- Sys.time()
      time_per_100 <- as.numeric(difftime(first_100_time, overall_start_time, units = "secs"))
      total_estimated_secs <- (time_per_100 / 100) * total_genes
      total_estimated_mins <- total_estimated_secs / 60
      
      cat(paste0("[진행 상황] 첫 100개 유전자 완료 (", round(time_per_100, 1), "초 소요)\n"))
      cat(paste0("   > 총 예상 소요 시간: 약 ", round(total_estimated_mins, 1), " 분\n"))
      
    } else if (i > 100 && i %% 100 == 0) {
      current_time <- Sys.time()
      elapsed_mins <- as.numeric(difftime(current_time, overall_start_time, units = "mins"))
      cat(paste0("[진행 상황] ", i, " / ", total_genes, " (", round(i/total_genes*100), "%) 완료. (", round(elapsed_mins, 1), "분 경과)\n"))
    }
  } # end for loop

  # 6. 결과 data.frame으로 변환
  spearman_results_df <- do.call(rbind, spearman_results_list)
  spearman_results_df <- as.data.frame(spearman_results_df)
  spearman_results_df$gene <- rownames(spearman_results_df)
  
  spearman_results_df <- spearman_results_df[!is.na(spearman_results_df$p.value), ]
  spearman_results_df <- spearman_results_df %>%
    dplyr::arrange(p.value) %>%
    dplyr::select(gene, rho, p.value)
  
  cat("스피어만 상관관계 분석 완료.\n")
  return(spearman_results_df)
}

#' @title Limma-Dream-SVA (LDS) 파이프라인
#' @description SVA로 숨겨진 공변량을 찾고, limma-dream으로 다중 임의 효과
#'              (예: 환자, 배치)를 포함한 LMM을 피팅합니다. 기본적으로 GeoMx data와 같이 샘플 수가 많지 않은 경우에 사용됩니다.
#'
#' @param sobj Seurat 객체, DGEList 객체, 또는 Raw Count Matrix.
#' @param formula (formula) 모델 포뮬러. LMM 수식(lme4)을 따릅니다.
#'                (예: ~ response + (1|emrid) + (1|set))
#' @param meta.data (data.frame) 메타데이터. `sobj`가 Matrix일 경우 필수.
#' @param layer (character) Seurat 객체 사용 시, Raw Count가 저장된 layer.
#' @param n_sv (numeric) 사용할 SV의 개수. `NULL` (기본값)일 경우,
#'             `sv_var_cutoff`에 따라 자동으로 결정됩니다.
#' @param sv_var_cutoff (numeric) `n_sv=NULL`일 때, SV가 설명해야 할
#'                      *잔차 분산(residual variance)*의 누적 비율 (0 ~ 1).
#' @param n_cores (numeric) 병렬 처리에 사용할 CPU 코어 수.
#'
#' @return (list) 'fit' (MArrayLM 객체), 'voom' (EList 객체),
#'         'sva' (SVA 객체), 'final_formula' (사용된 최종 포뮬러)
#'
LDS <- function(sobj,
                formula,
                meta.data = NULL,
                layer = "counts",
                n_sv = NULL,
                sv_var_cutoff = 0.5,
                n_cores = max(1, parallel::detectCores() - 2)) {

  # --- 1. 입력값 검증 및 데이터 추출 ---
  message("1/7: 입력 데이터 처리 중...")
  if (inherits(sobj, "Seurat")) {
    counts_matrix <- GetAssayData(sobj, layer = layer)
    if (is.null(meta.data)) {
      meta.data <- sobj@meta.data
    }
  } else if (inherits(sobj, "DGEList")) {
    counts_matrix <- sobj$counts
    if (is.null(meta.data)) {
      meta.data <- sobj$samples
    }
  } else if (is.matrix(sobj) || inherits(sobj, "dgCMatrix")) {
    counts_matrix <- sobj
    if (is.null(meta.data)) {
      stop("`sobj`가 Matrix일 경우, `meta.data`를 반드시 제공해야 합니다.")
    }
  } else {
    stop("`sobj`는 Seurat, DGEList, 또는 Matrix여야 합니다.")
  }
  
  if (ncol(counts_matrix) != nrow(meta.data)) {
    stop("Count matrix의 열 수(샘플)와 meta.data의 행 수(샘플)가 일치하지 않습니다.")
  }
  
  # 병렬 백엔드 설정
  BPPARAM_SETUP <- MulticoreParam(n_cores)

  # --- 2. SVA를 위한 포뮬러 파싱 ---
  message("2/7: 포뮬러 파싱 중...")
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("LMM 포뮬러 파싱을 위해 'lme4' 패키지가 필요합니다.")
  }
  
  # `formula`에서 고정 효과(Fixed Effects)만 추출
  # SVA의 'mod' 인자 및 `filterByExpr`에 사용됨
  fixed_effects_formula <- lme4::nobars(formula)
  
  if (is.null(fixed_effects_formula)) {
     # 예: ~ (1|emrid) + (1|set) 처럼 고정 효과가 아예 없는 경우
     fixed_effects_formula <- ~ 1 
  }

  # --- 3. DGEList 생성, 필터링, 정규화 ---
  message("3/7: DGEList 생성 및 필터링 중...")
  
  # `filterByExpr`를 위한 디자인 행렬 (고정 효과 기반)
  design_for_filter <- model.matrix(fixed_effects_formula, data = meta.data)
  
  dge <- DGEList(counts_matrix, samples = meta.data)
  
  keep_genes <- filterByExpr(dge, design = design_for_filter)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  
  if (sum(keep_genes) == 0) {
    stop("모든 유전자가 필터링되었습니다. `filterByExpr` 조건을 확인하십시오.")
  }
  
  dge <- calcNormFactors(dge)
  
  message(sprintf("... 유전자 필터링 완료: %d / %d 개 통과", 
                  sum(keep_genes), nrow(counts_matrix)))

  # --- 4. SVA 실행 (임시 voom 기반) ---
  message("4/7: SVA 실행 (숨겨진 변동성 탐색)...")
  
  mod_sva <- model.matrix(fixed_effects_formula, data = meta.data)
  mod0_sva <- model.matrix(~ 1, data = meta.data)
  
  # SVA는 voom 변환된 데이터에 실행
  v_sva <- voom(dge, mod_sva, plot = FALSE)
  
  # 1) n.sv=NULL로 SVA를 실행하여 *최대* SV 개수 및 SV *모두* 찾기
  # (n_sv=NULL일 때만 SVD를 수행하여 sva_obj$svd를 반환함. 
  #  -> 이 부분이 sva 버전마다 다를 수 있어, 수동 잔차 계산이 더 안정적임)
  sva_obj <- sva(v_sva$E, mod = mod_sva, mod0 = mod0_sva, n.sv = NULL)
  
  n_sv_max <- sva_obj$n.sv
  
  if (n_sv_max == 0) {
    message("... SVA가 유의미한 대리 변수(SV)를 찾지 못했습니다.")
    svs_final <- NULL
    n_sv_final <- 0
  } else {
    message(sprintf("... SVA가 최대 %d개의 SV를 탐지했습니다.", n_sv_max))
    
    # 2) 사용할 SV 개수 결정
    if (!is.null(n_sv)) {
      # 사용자가 SV 개수 명시
      n_sv_final <- min(n_sv, n_sv_max)
      
    } else {
      # (사용자 요청) 잔차 분산 기반으로 SV 개수 자동 결정
      message(sprintf("... 잔차 분산의 %.0f%%를 설명하는 SV 개수 자동 탐색 중...", 
                      sv_var_cutoff * 100))
                      
      # SVA가 사용한 것과 동일한 잔차(Residuals)를 수동 계산
      # lm.fit이 (샘플 x 유전자) 형태의 t(v_sva$E)를 입력받음
      res_matrix <- t(resid(lm.fit(mod_sva, t(v_sva$E))))
      
      # 잔차 행렬의 SVD (Singular Value Decomposition)
      svd_res <- svd(res_matrix)
      
      # 각 SV가 설명하는 분산 비율
      percent_var_explained <- (svd_res$d^2) / sum(svd_res$d^2)
      cumulative_var <- cumsum(percent_var_explained)
      
      # Cutoff를 만족하는 최소 SV 개수 찾기
      n_sv_auto <- which(cumulative_var >= sv_var_cutoff)[1]
      
      if (is.na(n_sv_auto)) { # 모든 SV를 합쳐도 cutoff 미만일 경우
        n_sv_final <- n_sv_max
      } else {
        n_sv_final <- n_sv_auto
      }
      
      message(sprintf("... SV %d개가 잔차 분산의 %.1f%%를 설명합니다.", 
                      n_sv_final, cumulative_var[n_sv_final] * 100))
    }
    
    # 3) 최종 SV 추출 및 메타데이터에 추가
    if (n_sv_final > 0) {
      svs_final <- sva_obj$sv[, 1:n_sv_final, drop = FALSE]
      colnames(svs_final) <- paste0("SV", 1:n_sv_final)
      meta.data <- cbind(meta.data, svs_final)
    }
  }

  # --- 5. 최종 포뮬러 생성 ---
  message("5/7: 최종 모델 포뮬러 생성 중...")
  
  original_formula_str <- paste(deparse(formula), collapse = "")
  
  if (n_sv_final > 0) {
    sv_str <- paste(colnames(svs_final), collapse = " + ")
    # 포뮬러에 SV 추가
    final_formula_str <- paste(original_formula_str, sv_str, sep = " + ")
  } else {
    final_formula_str <- original_formula_str
  }
  
  final_formula <- as.formula(final_formula_str)
  message(sprintf("... 최종 포뮬러: %s", final_formula_str))

  # --- 6. limma-dream 파이프라인 실행 ---
  message(sprintf("6/7: limma-dream 실행 (Core: %d개)...", n_cores))
  
  # 1) 유효 가중치 계산 (voomWithDreamWeights)
  v_dream <- voomWithDreamWeights(
    dge, 
    final_formula, 
    meta.data, 
    BPPARAM = BPPARAM_SETUP
  )
  
  # 2) LMM 피팅 (dream)
  fit_dream <- dream(
    v_dream, 
    final_formula, 
    meta.data, 
    BPPARAM = BPPARAM_SETUP
  )
  
  # 3) Empirical Bayes 조정
  fit_ebayes <- eBayes(fit_dream)
  
  message("7/7: 분석 완료.")

  # --- 7. 결과 반환 ---
  return(list(
    fit = fit_ebayes,       # 최종 MArrayLM 객체 (topTable 사용 가능)
    voom = v_dream,         # 가중치가 포함된 EList 객체
    sva_obj = sva_obj,      # 원본 SVA 객체
    svs_used = svs_final,   # 모델에 실제 사용된 SV 매트릭스
    final_formula = final_formula, # 최종 사용된 포뮬러
    dge = dge               # 필터링/정규화된 DGEList
  ))
}



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


#' @title Muscat (Pseudo-bulking) 파이프라인 함수 (V5 - prepData 사용)
#'
#' @description prepData를 사용하여 SCE 객체를 포맷하고, pbDS로 GLM(edgeR/DESeq2) 실행
#' @note 이 함수는 'pbDS'의 한계로 인해 혼합 효과(Random Effects/block)를 지원하지 않습니다.
#'
#' @param sobj (Seurat) Seurat 객체
#' @param cluster_id (character) 세포 유형 컬럼명
#' @param sample_id (character) 환자/샘플 ID 컬럼명
#' @param group_id (character) 비교할 그룹 컬럼명
#' @param formula_str (character) 포뮬러 (예: "~ 0 + group + set")
#'                 'group' 키워드는 'group_id'로 자동 변환됩니다.
#' @param contrast (character) Contrast (예: "groupA-groupB")
#' @param method (character) "edgeR", "DESeq2", "limma-trend", "limma-voom" (block 없음)
#'
#' @export
runMUSCAT <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,                 # ex) "exp_batch"
  contrast   = NULL,                 # ex) "IS - SAH"
  method     = "edgeR",
  pb_min_cells = 3,
  filter_genes = c("none","genes","both","edgeR"),
  keep_clusters = NULL,
  cluster_label_map = NULL
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")
  filter_genes <- match.arg(filter_genes)

  # deps
  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # 1) Seurat -> SCE, prepSCE
  sce <- Seurat::as.SingleCellExperiment(sobj)
  sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

  # factor 보장
  sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
  sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
  sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
    sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
  }

  # 2) Pseudobulk
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  # (선택) 특정 클러스터만
  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  # 2-1) pb 메타 보강 (sample_id / group_id / batch)
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))

  # sample_id 없으면 assay의 colnames로 복구
  if (!"sample_id" %in% names(pb_meta)) {
    first_assay <- names(SummarizedExperiment::assays(pb))[1]
    sid_guess <- colnames(SummarizedExperiment::assays(pb)[[first_assay]])
    if (is.null(sid_guess)) stop("pb에 sample_id가 없습니다.")
    pb_meta$sample_id <- sid_guess
    rownames(pb_meta) <- pb_meta$sample_id
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)
  }

  # sce에서 (sample_id -> group_id / batch) map
  sce_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  map_cols <- c("sample_id","group_id")
  if (!is.null(batch_id) && batch_id %in% names(sce_meta)) map_cols <- c(map_cols, batch_id)
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])

  # pb에 group_id / batch 보강
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  need_fix <- (!"group_id" %in% names(pb_meta)) ||
              (length(unique(pb_meta$group_id)) < 2) ||
              (all(unique(pb_meta$group_id) %in% c("type","group","group_id", NA, "")))
  if (need_fix || (!is.null(batch_id) && !batch_id %in% names(pb_meta))) {
    pb_meta2 <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
    if ("group_id.x" %in% names(pb_meta2) && "group_id.y" %in% names(pb_meta2)) {
      pb_meta2$group_id <- ifelse(is.na(pb_meta2$group_id.y), pb_meta2$group_id.x, pb_meta2$group_id.y)
      pb_meta2$group_id.x <- NULL; pb_meta2$group_id.y <- NULL
    }
    rownames(pb_meta2) <- rownames(pb_meta)
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta2)
  }

  # factor화
  pb$sample_id <- droplevels(factor(SummarizedExperiment::colData(pb)$sample_id))
  pb$group_id  <- droplevels(factor(SummarizedExperiment::colData(pb)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb))) {
    pb[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[batch_id]]))
  }

  # 3) contrast 그룹만 자동 subset
  extract_groups <- function(contrast_str, levels_available){
    z <- gsub("\\s+", "", contrast_str)
    toks <- unique(gsub("^group(_id)?", "", unlist(strsplit(z, "[^A-Za-z0-9_]+"))))
    toks <- toks[nchar(toks) > 0]
    keep <- intersect(toks, levels_available)
    if (length(keep) < 1) {
      g2 <- levels_available[vapply(levels_available, function(g) grepl(g, z), logical(1))]
      keep <- unique(g2)
    }
    keep
  }
  grp_lvls <- levels(SummarizedExperiment::colData(pb)$group_id)
  tg <- extract_groups(contrast, grp_lvls)
  if (length(tg) < 2) stop(sprintf("contrast에서 추출한 그룹이 부족합니다. contrast='%s', 사용가능레벨=%s",
                                   contrast, paste(grp_lvls, collapse=", ")))

  keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
  pb_sub <- pb[, keep_idx]
  pb_sub$group_id <- droplevels(factor(SummarizedExperiment::colData(pb_sub)$group_id))

  # **sce도 동일 기준으로 subset (resDS용 필수)**
  sce_sub <- sce[, sce$sample_id %in% SummarizedExperiment::colData(pb_sub)$sample_id &
                    sce$group_id  %in% tg]
  sce_sub$cluster_id <- droplevels(factor(sce_sub$cluster_id))
  sce_sub$sample_id  <- droplevels(factor(sce_sub$sample_id))
  sce_sub$group_id   <- droplevels(factor(sce_sub$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce_sub))) {
    sce_sub[[batch_id]] <- droplevels(factor(sce_sub[[batch_id]]))
  }

  # 4) design/contrast (batch는 'batch'로 복사해서 사용)
  pb_sub$group <- pb_sub$group_id
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb_sub))) {
    pb_sub$batch <- droplevels(factor(SummarizedExperiment::colData(pb_sub)[[batch_id]]))
    design <- stats::model.matrix(~ 0 + group + batch,
                                  data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
  } else {
    design <- stats::model.matrix(~ 0 + group,
                                  data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
  }

  fix_contrast <- function(contrast_str, design_cols){
    z <- gsub("\\s+", "", contrast_str)
    toks <- unlist(strsplit(z, "([+\\-])", perl=TRUE))
    ops  <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl=TRUE)))
    rebuild <- function(tok){
      tok <- gsub("^group(_id)?", "group", tok)
      if (!grepl("^group", tok)) tok <- paste0("group", tok)
      tok
    }
    toks2 <- vapply(toks, rebuild, character(1))
    out <- toks2[1]; if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i+1])
    out
  }
  contrast_fixed <- fix_contrast(contrast, colnames(design))
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)

  # 5) pbDS
  res <- muscat::pbDS(
    pb_sub,
    design    = design,
    method    = method,
    contrast  = contrast_matrix,
    min_cells = pb_min_cells,
    filter    = filter_genes,
    verbose   = TRUE
  )

  # 6) 결과 평탄화: **sce_sub가 먼저, res가 다음**
  combined <- muscat::resDS(sce_sub, res)

  # cluster_id 정리 + 라벨
  if ("cluster" %in% names(combined) && !"cluster_id" %in% names(combined)) {
    combined$cluster_id <- combined$cluster
  }
  if (!"cluster_id" %in% names(combined)) stop("resDS 결과에 'cluster_id'가 없습니다.")
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  # 7) 요약 헬퍼
  .pick_cols_for <- function(tab, contrast_fixed) {
    cstr <- gsub("\\s+", "", contrast_fixed)
    patt <- paste0("(^|\\.)", gsub("([+\\-])", "\\\\\\1", cstr), "$")

    # 1) 접미사 있는 형태 먼저 시도
    logfc <- grep("^logFC(\\.|_)?", names(tab), value=TRUE); logfc <- logfc[grep(patt, logfc)]
    padj  <- grep("^(p_adj(\\.loc|\\.glb)?|FDR)(\\.|_)?", names(tab), value=TRUE); padj <- padj[grep(patt, padj)]
    pcol  <- grep("^(p_val|PValue)(\\.|_)?", names(tab), value=TRUE); pcol <- pcol[grep(patt, pcol)]

    # 2) 접미사 매칭 실패 시 기본 컬럼으로 폴백
    if (!length(logfc) && "logFC" %in% names(tab)) logfc <- "logFC"
    if (!length(pcol)  && "p_val" %in% names(tab)) pcol  <- "p_val"
    if (!length(padj)) {
      if ("p_adj.loc" %in% names(tab))      padj <- "p_adj.loc"
      else if ("p_adj.glb" %in% names(tab)) padj <- "p_adj.glb"
      else if ("FDR" %in% names(tab))       padj <- "FDR"
    }

    if (!length(logfc) || !length(padj)) {
      stop("결과 테이블에서 logFC/adj.p 컬럼을 찾지 못했습니다. resDS 출력 컬럼명을 확인하세요.")
    }
    list(logfc = logfc[1], padj = padj[1], p = if (length(pcol)) pcol[1] else NULL)
  }

  # --- [REPLACE] 클러스터별 Top-N
  top_by_cluster <- function(n=25){
    tab_use <- combined
    # contrast 컬럼이 있으면 현재 contrast만 사용
    if ("contrast" %in% names(tab_use)) {
      target <- gsub("\\s+", "", contrast_fixed)
      tab_use <- tab_use[gsub("\\s+", "", tab_use$contrast) == target, , drop=FALSE]
    }
    cols <- .pick_cols_for(tab_use, contrast_fixed)

    tab_use |>
      dplyr::mutate(
        logFC_view = .data[[cols$logfc]],
        padj_view  = .data[[cols$padj]],
        p_view     = if (!is.null(cols$p)) .data[[cols$p]] else NA_real_
      ) |>
      dplyr::arrange(padj_view) |>
      dplyr::group_by(cluster_label) |>
      dplyr::slice_head(n = n) |>
      dplyr::ungroup()
  }

  # --- [REPLACE] 전체 Top-N
  top_overall <- function(n=100){
    tab_use <- combined
    if ("contrast" %in% names(tab_use)) {
      target <- gsub("\\s+", "", contrast_fixed)
      tab_use <- tab_use[gsub("\\s+", "", tab_use$contrast) == target, , drop=FALSE]
    }
    cols <- .pick_cols_for(tab_use, contrast_fixed)

    tab_use |>
      dplyr::mutate(
        logFC_view = .data[[cols$logfc]],
        padj_view  = .data[[cols$padj]],
        p_view     = if (!is.null(cols$p)) .data[[cols$p]] else NA_real_
      ) |>
      dplyr::arrange(padj_view) |>
      dplyr::slice_head(n = n)
  }

  list(
    sce        = sce,
    sce_sub    = sce_sub,
    pb         = pb,
    pb_sub     = pb_sub,
    design     = design,
    contrast_fixed = contrast_fixed,
    res_raw    = res,
    combined   = combined,
    top_by_cluster = top_by_cluster,
    top_overall    = top_overall
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