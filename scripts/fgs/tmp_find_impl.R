find_gene_signature_v5_impl <- function(data, 
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
                                 return_model = FALSE,
                                 fgs_seed = 42,
                                 # --- 신규/수정된 인자 ---
                                 lambda_selection = "lambda.1se",
                                 enet.alpha = 0.5,        # (Elastic Net용)
                                 pca.n_pcs = 1,           # (PCA용)
                                 gam.min_unique = 15,     # (GAM용)
                                 gam.k = NULL,
                                 gam.k_dynamic_factor = 5,
                                 method_impl = c("v5.3","v5.4"),
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
  methods_requiring_scale <- c("lasso", "ridge", "elastic_net", 
                               "gam", "pca_loadings", "xgboost")
  methods_requiring_correction <- c("wilcoxon", "pca_loadings")
  
  preprocessed_data <- fgs_preprocess_data_v5(
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
      } else { auc <- NA }
      acc <- mean(pred == y)
      perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
    } else {
      pred_class_indices <- apply(pred_probs[,,1], 1, which.max)
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
        
        random_forest_ranger = {
          if (!requireNamespace("ranger", quietly = TRUE)) {
            stop("ranger package required.")
          }
          
          X_ranger <- t(expr_mat_base) 
          if (!is.null(control_vars)) { X_ranger <- cbind(X_ranger, covariate_mat_model) }
          
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
          if (!is.null(control_vars)) {
            weights_magnitude_all <- weights_magnitude_all[!names(weights_magnitude_all) %in% colnames(covariate_mat_model)]
          }
          weights_magnitude_all <- weights_magnitude_all[is.finite(weights_magnitude_all)]
          
          top_genes <- names(sort(weights_magnitude_all, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_all))])
          weights_magnitude <- weights_magnitude_all[top_genes]

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
          
          scores <- as.numeric(X_ranger[, top_genes, drop=FALSE] %*% weights)
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

        nmf_loadings = {
          if (!requireNamespace("NMF", quietly = TRUE)) {
            stop("NMF package required.")
          }
          
          expr_mat_nmf <- expr_mat_method 
          
          if (!is.null(control_vars)) {
            warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s...", m))
            covariate_mat <- model.matrix(~ . - 1, data = meta.data_clean[, control_vars, drop=FALSE])
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
          if (!requireNamespace("mgcv", quietly = TRUE)) {
            stop("mgcv package required. Install with: install.packages('mgcv')")
          }
          
          deviance_explained <- numeric(nrow(expr_mat_method))
          names(deviance_explained) <- rownames(expr_mat_method)
          
          gam_data <- data.frame(y_var_numeric = as.numeric(target_binary) - 1)
          if (!is.null(control_vars)) {
            gam_data <- cbind(gam_data, covariate_mat_model)
            formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse=" + "))
          } else {
            formula_base <- ""
          }

          convergence_warnings <- 0
          genes_skipped <- 0
          genes_with_dynamic_k <- 0
          
          for (i in 1:nrow(expr_mat_method)) {
            gam_data$gene_expr <- expr_mat_method[i, ]
            
            n_unique_vals <- length(unique(gam_data$gene_expr))
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

            fit_result <- tryCatch({
              if (n_groups == 2) {
                mgcv::bam(as.formula(formula_str), data = gam_data, family="binomial", ...)
              } else {
                mgcv::bam(as.formula(formula_str), data = gam_data, ...)
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
          if (use_dynamic_k && genes_with_dynamic_k > 0) {
            message(sprintf("GAM: Applied dynamic k for %d genes.", genes_with_dynamic_k))
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
          
          scores <- as.numeric(X[, top_genes, drop=FALSE] %*% weights)
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
