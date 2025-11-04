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




.fetch_plot_data <- function(data, features, sample_col = NULL, group.by = NULL, split.by = NULL, 
                             assay = NULL, layer = "data") {
  
  `%||%` <- rlang::`%||%`
  
  # 1. Seurat 객체 처리
  if (inherits(data, "Seurat")) {
    sobj <- data
    meta_cols <- names(sobj@meta.data)
    
    meta_features <- intersect(features, meta_cols)
    assay_features <- setdiff(features, meta_cols)
    
    # 필요한 모든 컬럼 정의
    grouping_cols <- unique(c(sample_col, group.by, split.by))
    grouping_cols <- grouping_cols[!is.null(grouping_cols)]
    
    vars_to_fetch <- unique(c(grouping_cols, meta_features, assay_features))
    vars_to_fetch <- vars_to_fetch[vars_to_fetch %in% c(meta_cols, rownames(sobj[[assay %||% DefaultAssay(sobj)]]))] # 존재하는 것만
    
    # 데이터 추출
    fetched_df <- data.frame(matrix(ncol = 0, nrow = ncol(sobj)))
    rownames(fetched_df) <- colnames(sobj)
    
    if (length(vars_to_fetch) > 0) {
        if (length(assay_features) > 0) {
          assay_to_use <- assay %||% DefaultAssay(sobj)
          fetched_df <- FetchData(sobj, vars = vars_to_fetch, assay = assay_to_use, layer = layer)
        } else {
          fetched_df <- FetchData(sobj, vars = vars_to_fetch)
        }
    } else {
       warning("가져올 유효한 피처나 메타데이터가 없습니다.")
       return(data.frame())
    }

  # 2. 데이터 프레임 처리
  } else if (is.data.frame(data)) {
    fetched_df <- data
    req_cols <- unique(c(features, sample_col, group.by, split.by))
    req_cols <- req_cols[!is.null(req_cols)]
    
    missing_cols <- setdiff(req_cols, names(fetched_df))
    if (length(missing_cols) > 0) {
      stop("데이터 프레임에 다음 컬럼이 없습니다: ", paste(missing_cols, collapse=", "))
    }
  } else {
    stop("`data` 인수는 Seurat 객체 또는 data.frame이어야 합니다.")
  }
  
  return(fetched_df)
}

#' 세포 단위 피처 발현 박스플롯 (plot_gene_boxplot 스타일)
#'
#' @param data Seurat 객체 또는 data.frame
#' @param features 플로팅할 피처 (유전자 또는 메타데이터)
#' @param group.by 패널(facet)로 분리할 그룹 변수 (예: "cell_type")
#' @param split.by X축에서 나눌 그룹 변수 (예: "timepoint")
#' @param idents 'group.by'에서 포함할 특정 값들
#' @param assay Seurat 객체용 assay
#' @param layer Seurat 객체용 layer
#' @param ncol 플롯 컬럼 수
#' @param pt.size 점 크기 (0이면 숨김)
#' @param violin 바이올린 플롯 오버레이
#' @param add_stats 통계 검정 추가 (ggpubr::stat_compare_means)
#' @param ... stat_compare_means에 전달할 추가 인수 (예: paired = TRUE, method = "t.test")
#'
#' @return ggplot 객체
#' @export
mybox_cell <- function(data, 
                                      features, 
                                      group.by, 
                                      split.by, 
                                      idents = NULL, 
                                      assay = NULL, 
                                      layer = "data",
                                      ncol = 3, 
                                      pt.size = 0.2, 
                                      violin = FALSE,
                                      add_stats = TRUE,
                                      ...) {

  plot_df_base <- .fetch_plot_data(data, features, NULL, group.by, split.by, assay, layer)
  
  all_plots <- list()
  
  for (feature_name in features) {
    current_df <- plot_df_base
    
    # 1. group.by가 NULL일 경우 "Overall"로 처리
    group.by.internal <- group.by
    if (is.null(group.by)) {
      group.by.internal <- ".internal_placeholder_group"
      current_df[[group.by.internal]] <- "Overall"
    }
    
    # 2. idents 필터링 (mybox 로직과 동일)
    if (!is.null(group.by) && !is.null(idents)) {
      original_col_type <- class(current_df[[group.by.internal]])
      current_df[[group.by.internal]] <- as.character(current_df[[group.by.internal]])
      current_df <- current_df[current_df[[group.by.internal]] %in% idents, ]
      
      if (nrow(current_df) == 0) {
         warning(feature_name, "에서 idents 필터링 후 데이터가 없습니다. 건너뜁니다.")
         next
      }
      
      if ("factor" %in% original_col_type) {
        original_levels <- levels(if(inherits(data, "Seurat")) data@meta.data[[group.by]] else data[[group.by]])
        valid_subset_levels <- intersect(original_levels, idents)
        current_df[[group.by.internal]] <- factor(current_df[[group.by.internal]], levels = valid_subset_levels)
      } else {
        current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
      }
    }
    
    current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
    current_df[[split.by]] <- as.factor(current_df[[split.by]])
    
    # 3. ggplot 생성 (Req 3 레이아웃)
    p <- ggplot(current_df, aes(x = .data[[split.by]], y = .data[[feature_name]], fill = .data[[split.by]]))
    
    if (violin) {
      p <- p + geom_violin(alpha = 0.5, scale = "width", trim = FALSE, na.rm = TRUE)
    }
    
    # outlier.shape=NA로 설정 (Req 2)
    p <- p + geom_boxplot(na.rm = TRUE, outlier.shape = NA, alpha = 0.7)
    
    # pt.size > 0 일 때만 jitter 추가 (Req 2: 색상 "black" 고정)
    if (pt.size > 0) {
      p <- p + geom_jitter(color = "black", size = pt.size, alpha = 0.3, 
                           width = 0.2, height = 0, na.rm = TRUE)
    }
    
    # 통계 추가 (ggpubr)
    if (add_stats) {
      # `...` 인수로 paired=TRUE 등을 받을 수 있음
      p <- p + stat_compare_means(label = "p.format", ...) 
    }
    
    # 패널(facet) 적용 및 테마 설정 (Req 3)
    p <- p + facet_wrap(as.formula(paste("~", group.by.internal))) +
      labs(title = feature_name,
           x = split.by,
           y = feature_name) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none") # 범례 제거
    
    all_plots[[feature_name]] <- p
  }
  
  if (length(all_plots) == 0) {
    message("생성된 플롯이 없습니다.")
    return(invisible(NULL))
  }
  
  return(patchwork::wrap_plots(all_plots, ncol = ncol))
}
#' 샘플 단위(Pseudobulk) 피처 평균 박스플롯 (mybox 로직 + plot_gene_boxplot 스타일)
#'
#' @param data Seurat 객체 또는 data.frame
#' @param features 플로팅할 피처 (유전자 또는 메타데이터)
#' @param sample_col 샘플/환자 식별 컬럼 (필수)
#' @param group.by 패널(facet)로 분리할 그룹 변수 (예: "cell_type")
#' @param split.by X축에서 나눌 그룹 변수 (예: "timepoint")
#' @param idents 'group.by'에서 포함할 특정 값들
#' @param assay Seurat 객체용 assay
#' @param layer Seurat 객체용 layer
#' @param ncol 플롯 컬럼 수
#' @param pt.size 샘플 평균 점 크기 (0이면 숨김)
#' @param violin 바이올린 플롯 오버레이
#' @param add_stats 통계 검정 추가 (요약된 데이터 대상)
#' @param ... stat_compare_means에 전달할 추가 인수 (예: paired = TRUE)
#'
#' @return ggplot 객체
#' @export
mybox_pb <- function(data, 
                                    features, 
                                    sample_col, 
                                    group.by, 
                                    split.by, 
                                    idents = NULL, 
                                    assay = NULL, 
                                    layer = "data",
                                    ncol = 3, 
                                    pt.size = 1.0, 
                                    violin = FALSE,
                                    add_stats = FALSE,
                                    ...) {

  if (is.null(sample_col)) {
    stop("`sample_col` 인수는 Pseudobulk 플롯에 필수입니다.")
  }

  plot_df_base <- .fetch_plot_data(data, features, sample_col, group.by, split.by, assay, layer)
  
  all_plots <- list()
  
  for (feature_name in features) {
    current_df <- plot_df_base
    
    # 1. group.by가 NULL일 경우 "Overall"로 처리
    group.by.internal <- group.by
    if (is.null(group.by)) {
      group.by.internal <- ".internal_placeholder_group"
      current_df[[group.by.internal]] <- "Overall"
    }
    
    # 2. idents 필터링
    if (!is.null(group.by) && !is.null(idents)) {
      original_col_type <- class(current_df[[group.by.internal]])
      current_df[[group.by.internal]] <- as.character(current_df[[group.by.internal]])
      current_df <- current_df[current_df[[group.by.internal]] %in% idents, ]
      
      if (nrow(current_df) == 0) {
         warning(feature_name, "에서 idents 필터링 후 데이터가 없습니다. 건너뜁니다.")
         next
      }
      
      if ("factor" %in% original_col_type) {
        original_levels <- levels(if(inherits(data, "Seurat")) data@meta.data[[group.by]] else data[[group.by]])
        valid_subset_levels <- intersect(original_levels, idents)
        current_df[[group.by.internal]] <- factor(current_df[[group.by.internal]], levels = valid_subset_levels)
      } else {
        current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
      }
    }
    
    current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
    current_df[[split.by]] <- as.factor(current_df[[split.by]])
    current_df[[sample_col]] <- as.factor(current_df[[sample_col]])
    
    # 3. *** 데이터 요약 (Aggregation) ***
    grouping_vars <- unique(c(sample_col, group.by.internal, split.by))
    
    aggregated_df <- current_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::summarise(mean_value = mean(.data[[feature_name]], na.rm = TRUE), .groups = 'drop')
    
    if (nrow(aggregated_df) == 0) {
      warning(feature_name, "에서 요약 후 데이터가 없습니다. 건너뜁니다.")
      next
    }
    
    # 4. ggplot 생성 (Req 3 레이아웃)
    p <- ggplot(aggregated_df, aes(x = .data[[split.by]], y = mean_value, fill = .data[[split.by]]))
    
    if (violin) {
      p <- p + geom_violin(alpha = 0.5, scale = "width", trim = FALSE, na.rm = TRUE)
    }
    
    # outlier.shape=NA로 설정 (Req 2)
    p <- p + geom_boxplot(na.rm = TRUE, outlier.shape = NA, alpha = 0.7)
    
    # pt.size > 0 일 때만 jitter 추가 (Req 2: 색상 "black" 고정)
    if (pt.size > 0) {
      # PB 플롯은 점이 샘플을 의미하므로 더 진하게 (alpha=0.7)
      p <- p + geom_jitter(color = "black", size = pt.size, alpha = 0.7, 
                           width = 0.2, height = 0, na.rm = TRUE)
    }

    # 통계 추가 (요약된 데이터를 대상으로 함)
    if (add_stats) {
      # paired=TRUE를 `...`로 전달받아 샘플 페어 테스트 가능
      p <- p + stat_compare_means(label = "p.format", ...) 
    }

    # 패널(facet) 적용 및 테마 설정 (Req 3)
    p <- p + facet_wrap(as.formula(paste("~", group.by.internal))) +
      labs(title = feature_name,
           x = split.by,
           y = paste("Average", feature_name)) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none") # 범례 제거
            
    all_plots[[feature_name]] <- p
  }
  
  if (length(all_plots) == 0) {
    message("생성된 플롯이 없습니다.")
    return(invisible(NULL))
  }
  
  return(patchwork::wrap_plots(all_plots, ncol = ncol))
}