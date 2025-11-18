
# # analysis -----
# NEBULA=function(){}


# #' @title MAST 파이프라인 함수
# #'
# #' @description [수정] Seurat -> SCE로 변환하여 MAST를 실행
# #'
# #' @param sobj (Seurat) Seurat 객체
# #' @param formula (formula or character) lme4 문법의 포뮬러
# #' @param min_cells_expr (numeric) 유전자 필터링 기준 (최소 발현 세포 수)
# #' @param n_cores (numeric) 병렬 처리 코어 수
# #' @param lrt_variable (character) LRT 검정을 수행할 변수명 (예: "type")
# #'
# #' @export
# runMAST <- function(sobj,
#                     formula,
#                     min_cells_expr = 10,
#                     n_cores = 4,
#                     lrt_variable = NULL) {
  
#   if (!requireNamespace("MAST", quietly = TRUE)) stop("MAST 패키지가 필요합니다.")
#   if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("SCE 패키지가 필요합니다.")
#   if (is.null(lrt_variable)) stop("'lrt_variable' 인자 (예: 'g3')를 지정해야 합니다.")
  
#   if (is.character(formula)) {
#     formula_obj <- as.formula(formula)
#   } else if (inherits(formula, "formula")) {
#     formula_obj <- formula
#   } else {
#     stop("'formula'는 문자열 또는 formula 객체여야 합니다.")
#   }
  
#   # --- 1. SCA 객체 생성 및 정규화 ---
#   message("1/5: Seurat -> SingleCellExperiment(SCE) 객체 변환 중...")
  
#   # [수정] MAST Problem 1 해결: FromSeurat 대신 as.SingleCellExperiment 사용
#   sca <- as.SingleCellExperiment(sobj)
  
#   # --- 2. 유전자 필터링 ---
#   message(sprintf("2/5: 유전자 필터링 (min %d cells)...", min_cells_expr))
#   # freq()는 발현 비율 (0~1)을 반환
#   keep_genes <- (MAST::freq(sca) * ncol(sca)) >= min_cells_expr
#   sca_filtered <- sca[keep_genes, ]
#   message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(sca)))

#   # --- 3. 정규화 (MAST는 log2cpm 사용) ---
#   message("3/5: Log2(CPM+1) 정규화 중...")
#   SummarizedExperiment::assay(sca_filtered, "logcpm") <- MAST::cpm(sca_filtered, log = TRUE)
  
#   # --- 4. zlm (Hurdle LMM) 실행 ---
#   message(sprintf("4/5: MAST::zlm 실행 (Cores: %d). 시간이 오래 걸릴 수 있습니다...", n_cores))
  
#   zfit <- MAST::zlm(formula_obj, 
#                     sca = sca_filtered, 
#                     method = "glmer", 
#                     parallel = TRUE,
#                     nCores = n_cores)
  
#   # --- 5. LRT 결과 요약 ---
#   message(sprintf("5/5: LRT 검정 수행 (변수: %s)...", lrt_variable))
#   summary_res <- summary(zfit, doLRT = lrt_variable)
#   summary_dt <- summary_res$datatable
  
#   results_df <- merge(
#       summary_dt[component == 'H', .(primerid, `Pr(>Chisq)`)], # Hurdle (logistic)
#       summary_dt[component == 'logcpm', .(primerid, coef, ci.hi, ci.lo)], # Continuous
#       by = 'primerid'
#   )
#   colnames(results_df)[2] <- "p_value_hurdle"
#   results_df <- results_df[order(p_value_hurdle), ]
  
#   message("MAST 분석 완료.")
#   return(results_df)
# }

# # [V4] find_gene_signature_v4
# #
# # 변경점:
# # 1. control_vars = NULL 인자 추가 (교란 변수 보정)
# # 2. test_n = NULL 인자 추가 (속도 향상을 위한 사전 필터링)
# # 3. limma, lasso, tree_based, gam: Native 방식 보정
# # 4. wilcoxon, nmf, pca_loadings: limma::removeBatchEffect 방식 보정

# #' @export
# find_gene_signature_v4 <- function(data, 
#                                  meta.data = NULL,
#                                  target_var,
#                                  target_group = NULL,
#                                  control_vars = NULL,   # <<< [V4 UPGRADE]
#                                  method = c("tree_based", "lasso", "limma", 
#                                             "nmf", "wilcoxon", "gam", "pca_loadings"),
#                                  n_features = 50,
#                                  test_n = NULL,         # <<< [V4 UPGRADE]
#                                  preprocess = TRUE,
#                                  min_cells = 10,
#                                  min_pct = 0.01,
#                                  return_model = FALSE,
#                                  fgs_seed = 42,
#                                  lambda_selection = "lambda.1se",
#                                  ...) {
  
#   # ===
#   # 0. Batch Processing (v3와 동일)
#   # ===
  
#   all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
#   if (is.null(method)) {
#     method <- all_methods
#   }
  
#   if (length(method) > 1) {
#     current_call <- match.call()
    
#     results_list <- lapply(method, function(m) {
#       single_call <- current_call
#       single_call$method <- m
#       tryCatch({
#         eval(single_call)
#       }, error = function(e) {
#         warning(sprintf("Method '%s' failed with error: %s", m, e$message))
#         return(list(method = m, error = e$message))
#       })
#     })
    
#     names(results_list) <- method
#     return(results_list)
#   }
  
#   if (!method %in% all_methods) {
#     stop(sprintf("Invalid method '%s'. Choose from: %s", 
#                  method, paste(all_methods, collapse=", ")))
#   }
  
#   set.seed(fgs_seed)
  
#   # ===
#   # 1. Input validation and data extraction (v3와 동일)
#   # ===
  
#   is_seurat <- inherits(data, "Seurat")
  
#   if (is_seurat) {
#     if (!requireNamespace("Seurat", quietly = TRUE)) {
#       stop("Seurat package required but not installed")
#     }
#     if (is.null(meta.data)) {
#       meta.data <- data@meta.data
#     }
#     if ("data" %in% names(data@assays[[1]])) {
#       expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
#     } else {
#       expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
#     }
#   } else {
#     if (is.null(meta.data)) {
#       stop("meta.data must be provided when data is not a Seurat object")
#     }
#     expr_mat <- as.matrix(data)
#     if (nrow(expr_mat) == nrow(meta.data)) {
#       expr_mat <- t(expr_mat)
#     }
#   }
  
#   if (!target_var %in% colnames(meta.data)) {
#     stop(sprintf("target_var '%s' not found", target_var))
#   }
  
#   # [V4 UPGRADE] control_vars 유효성 검사
#   if (!is.null(control_vars)) {
#     missing_vars <- control_vars[!control_vars %in% colnames(meta.data)]
#     if (length(missing_vars) > 0) {
#       stop(sprintf("control_vars not found in metadata: %s", 
#                    paste(missing_vars, collapse=", ")))
#     }
#   }
  
#   common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
#   if (length(common_cells) == 0) {
#     stop("No common cells found between expression matrix and metadata")
#   }
#   expr_mat <- expr_mat[, common_cells]
#   meta.data <- meta.data[common_cells, ]
  
#   # ===
#   # 2. Process target variable (v3와 동일)
#   # ===
  
#   target_values <- meta.data[[target_var]]
#   target_type <- class(target_values)[1]
  
#   if (is.numeric(target_values)) {
#     # ... (numeric 처리 로직, v3와 동일) ...
#      if (is.null(target_group)) { target_group <- 0.25 }
#      if (is.list(target_group)) {
#        low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
#        high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
#      } else if (length(target_group) == 1 && target_group < 1) {
#        low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
#        high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
#      } else {
#        low_cutoff <- high_cutoff <- target_group
#      }
#      group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                             ifelse(target_values >= high_cutoff, "High", NA))
#      keep_cells <- !is.na(group_labels)
#   } else {
#     if (!is.null(target_group)) {
#       keep_cells <- target_values %in% target_group
#     } else {
#       keep_cells <- TRUE # 모든 팩터 레벨 사용
#     }
#     group_labels <- target_values
#   }
  
#   expr_mat <- expr_mat[, keep_cells]
#   meta.data <- meta.data[keep_cells, ]
#   target_binary <- factor(group_labels[keep_cells]) # factor로 변환
  
#   if (length(unique(target_binary)) < 2) {
#     stop("Target variable must have at least 2 groups after processing")
#   }
  
#   n_groups <- length(unique(target_binary))
  
#   # ===
#   # 3. Filter and preprocess genes (v3와 동일)
#   # ===
  
#   n_cells_expr <- rowSums(expr_mat > 0)
#   pct_cells_expr <- n_cells_expr / ncol(expr_mat)
#   keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
#   expr_mat <- expr_mat[keep_genes, ]
  
#   if (nrow(expr_mat) == 0) {
#     stop("No genes pass filtering criteria")
#   }
  
#   # ===
#   # 3.5 [V4 UPGRADE] Speed-up: Pre-filter genes using test_n
#   # ===
#   if (!is.null(test_n) && nrow(expr_mat) > test_n) {
#     message(sprintf("Pre-filtering to top %d genes based on limma p-value...", test_n))
#     if (!requireNamespace("limma", quietly = TRUE)) {
#       stop("limma package required for 'test_n' pre-filtering")
#     }
    
#     # limma는 교란변수를 '미리' 보정하고 p-value를 뽑을 수 있음
#     if (is.null(control_vars)) {
#       design_test <- model.matrix(~ target_binary, data = meta.data)
#     } else {
#       formula_test <- as.formula(paste("~ target_binary +", paste(control_vars, collapse="+")))
#       design_test <- model.matrix(formula_test, data = meta.data)
#     }
    
#     fit_test <- limma::lmFit(expr_mat, design_test)
#     fit_test <- limma::eBayes(fit_test)
    
#     # p-value 추출 (target_binary 그룹 비교)
#     top_table_test <- limma::topTable(fit_test, 
#                                       coef = grep("target_binary", colnames(design_test)), 
#                                       number = Inf, sort.by = "P")
    
#     top_gene_names <- rownames(top_table_test)[1:min(test_n, nrow(top_table_test))]
#     expr_mat <- expr_mat[top_gene_names, ]
#     message(sprintf("... reduced to %d genes.", nrow(expr_mat)))
#   }
  
#   if (preprocess) {
#     if (max(expr_mat) > 100) {
#       expr_mat <- log1p(expr_mat)
#     }
#     if (method %in% c("lasso", "gam", "pca_loadings")) {
#       gene_means <- rowMeans(expr_mat)
#       gene_sds <- apply(expr_mat, 1, sd)
#       gene_sds[gene_sds == 0] <- 1
#       expr_mat <- (expr_mat - gene_means) / gene_sds
#     }
#   }
  
#   # ===
#   # 3.8 [V4 UPGRADE] Confounder pre-correction (for non-native methods)
#   # ===
#   if (!is.null(control_vars) && method %in% c("wilcoxon", "nmf", "pca_loadings")) {
#     if (!requireNamespace("limma", quietly = TRUE)) {
#       stop("limma package required for removeBatchEffect")
#     }
    
#     message(sprintf("Method %s: Applying limma::removeBatchEffect for: %s",
#                     method, paste(control_vars, collapse=", ")))
    
#     # removeBatchEffect는 design matrix가 아닌 'batch' 또는 'covariates'를 받음
#     covariates_df <- meta.data[, control_vars, drop = FALSE]
    
#     # Factor 변수들을 model.matrix로 변환
#     covariate_mat <- model.matrix(~ . - 1, data = covariates_df) 
    
#     expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
#   }
  
  
#   # ===
#   # 4. Method-specific signature discovery
#   # ===
  
#   # X (samples x genes)는 여러 메서드에서 공통으로 사용
#   X <- t(expr_mat)
#   y <- target_binary
  
#   # [V4 UPGRADE] Native-supporting methods에서 사용할 교란변수 매트릭스
#   if (!is.null(control_vars) && method %in% c("tree_based", "lasso", "gam")) {
#      covariate_mat_model <- model.matrix(~ . - 1, 
#                                 data = meta.data[, control_vars, drop=FALSE])
#   }

  
#   result <- switch(method,
                 
#                  tree_based = {
#                    if (!requireNamespace("randomForest", quietly = TRUE)) { ... }
                   
#                    X_rf <- X
                   
#                    # [V4 UPGRADE] 교란 변수를 X 매트릭스에 추가
#                    if (!is.null(control_vars)) {
#                      X_rf <- cbind(X_rf, covariate_mat_model)
#                    }
                   
#                    if (ncol(X_rf) > 2000) {
#                      # ... (v3의 top_var_genes 로직은 test_n으로 대체됨) ...
#                      # ... 단, control_vars가 추가되었으므로 컬럼 이름으로 필터링
#                      gene_vars <- apply(X, 2, var) # 원본 X에서만 var 계산
#                      top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
                     
#                      # 교란 변수는 항상 포함
#                      X_rf <- X_rf[, c(top_var_genes, colnames(covariate_mat_model))]
#                    }
                   
#                    rf_model <- randomForest::randomForest(x = X_rf, y = y, ntree = 500, importance = TRUE, ...)
                   
#                    importance_scores <- randomForest::importance(rf_model)
                   
#                    # ... (v3의 MeanDecreaseGini 로직) ...
#                    if (n_groups == 2) {
#                      weights_magnitude_all <- importance_scores[, "MeanDecreaseGini"]
#                    } else {
#                      weights_magnitude_all <- rowMeans(importance_scores[, grep("MeanDecreaseGini", colnames(importance_scores))])
#                    }
                   
#                    # [V4 UPGRADE] 교란 변수를 제외하고 '유전자'만 선택
#                    gene_names_in_model <- colnames(X_rf)[!colnames(X_rf) %in% colnames(covariate_mat_model)]
#                    weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]
                   
#                    top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
#                    weights_magnitude <- weights_magnitude_genes[top_genes]
                   
#                    # ... (v3의 방향성 보정 및 나머지 로직 동일) ...
#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]; g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      effect_size <- mean_g2 - mean_g1 
#                      weights <- weights_magnitude * sign(effect_size)
#                    } else {
#                      weights <- weights_magnitude
#                    }
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    # ... (v3의 performance 로직) ...
#                    list(...) 
#                  },
                 
#                  lasso = {
#                    if (!requireNamespace("glmnet", quietly = TRUE)) { ... }
                   
#                    # [V4 UPGRADE] X 매트릭스 구성 및 penalty.factor 설정
#                    if (is.null(control_vars)) {
#                      X_model <- X
#                      penalty_vec <- rep(1, ncol(X_model)) # 모든 유전자에 페널티 1
#                    } else {
#                      X_model <- cbind(X, covariate_mat_model)
#                      # 유전자는 페널티 1, 교란 변수는 페널티 0 (Lasso에서 제외)
#                      penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
#                    }
                   
#                    # [V4 UPGRADE] cv.glmnet 호출 시 penalty.factor 전달
#                    if (n_groups == 2) {
#                      cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=1, 
#                                                  penalty.factor = penalty_vec, ...)
#                    } else {
#                      cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=1, 
#                                                  penalty.factor = penalty_vec, ...)
#                    }
                   
#                    coefs <- coef(cv_fit, s = lambda_selection) 
                   
#                    # [V4 UPGRADE] weights 추출 시 유전자 부분만 선택
#                    if (n_groups == 2) {
#                      # 1번은 intercept, (ncol(X)+1)까지가 유전자
#                      weights_all <- as.numeric(coefs[2:(ncol(X)+1)]) 
#                      names(weights_all) <- rownames(coefs)[2:(ncol(X)+1)]
#                    } else {
#                      coef_list <- lapply(coefs, function(x) as.numeric(x[2:(ncol(X)+1)]))
#                      weights_all <- rowMeans(do.call(cbind, coef_list))
#                      names(weights_all) <- rownames(coefs[[1]])[2:(ncol(X)+1)]
#                    }
                   
#                    # ... (v3의 non-zero 유전자 선택 및 나머지 로직 동일) ...
#                    nonzero_genes <- names(weights_all)[weights_all != 0]
#                    if (length(nonzero_genes) == 0) {
#                      top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
#                    } else {
#                      nonzero_weights <- weights_all[nonzero_genes]
#                      top_genes <- names(sort(abs(nonzero_weights), decreasing=TRUE)[1:min(n_features, length(nonzero_weights))])
#                    }
#                    weights <- weights_all[top_genes]
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    # ... (v3의 performance 로직) ...
#                    list(...)
#                  },
                 
#                  limma = {
#                    if (!requireNamespace("limma", quietly = TRUE)) { ... }
                   
#                    # [V4 UPGRADE] control_vars를 포함하는 formula 생성
#                    if (is.null(control_vars)) {
#                      design <- model.matrix(~0 + target_binary, data = meta.data)
#                      colnames(design)[1:n_groups] <- levels(target_binary)
#                    } else {
#                      control_formula <- paste(control_vars, collapse = " + ")
#                      full_formula <- as.formula(paste("~0 + target_binary +", control_formula))
#                      design <- model.matrix(full_formula, data = meta.data)
#                      colnames(design)[1:n_groups] <- levels(target_binary) 
#                    }

#                    fit <- limma::lmFit(expr_mat, design)
                   
#                    # ... (v3의 contrast 생성 및 나머지 로직 동일) ...
#                    if (n_groups == 2) {
#                      contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
#                    } else {
#                       # ... (v3의 n_groups > 2 로직) ...
#                    }
#                    contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
#                    fit2 <- limma::contrasts.fit(fit, contrast_mat)
#                    fit2 <- limma::eBayes(fit2)
#                    top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
                   
#                    weights_all <- top_table$t
#                    names(weights_all) <- rownames(top_table)
#                    top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
#                    weights <- weights_all[top_genes]
#                    scores <- as.numeric(X[, top_genes] %*% weights)
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = list(top_table = top_table[top_genes, ]), 
#                         model = if(return_model) fit2 else NULL)
#                  },
                 
#                  wilcoxon = {
#                    # [V4] 이 블록은 수정할 필요가 없습니다.
#                    # 만약 control_vars가 있었다면, 함수 상단의
#                    # '3.8 Pre-correction' 블록에서 expr_mat이 이미 보정되었습니다.
                   
#                    # ... (v3의 wilcoxon 로직 동일) ...
#                    pvals <- numeric(nrow(expr_mat))
#                    effect_sizes <- numeric(nrow(expr_mat))
#                    for (i in 1:nrow(expr_mat)) { ... }
#                    # ... (v3 로직 계속) ...
#                    list(...)
#                  },
                 
#                  nmf = {
#                    # [V4] 이 블록도 수정할 필요가 없습니다.
#                    # control_vars가 있었다면 expr_mat이 이미 보정되었습니다.
                   
#                    # ... (v3의 nmf 로직 동일) ...
#                    list(...)
#                  },
                 
#                  gam = {
#                    if (!requireNamespace("mgcv", quietly = TRUE)) { ... }
                   
#                    # ... (v3의 n_test_genes 로직 동일) ...
                   
#                    deviance_explained <- numeric(nrow(expr_mat))
#                    names(deviance_explained) <- rownames(expr_mat)
                   
#                    # [V4 UPGRADE] 교란 변수를 포함하는 GAM formula 생성
#                    if (is.null(control_vars)) {
#                      formula_base <- ""
#                      gam_data <- data.frame(y_var = target_binary)
#                    } else {
#                      formula_base <- paste(" +", paste(control_vars, collapse=" + "))
#                      gam_data <- data.frame(y_var = target_binary, 
#                                             meta.data[, control_vars, drop=FALSE])
#                    }
                   
#                    # ... (v3의 gam 루프) ...
#                    for (i in test_genes_idx) {
#                      gam_data$gene_expr <- expr_mat[i, ]
                     
#                      if (n_groups == 2) {
#                        gam_data$y_var_numeric <- as.numeric(gam_data$y_var) - 1
#                        full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
#                        gam_fit <- mgcv::gam(as.formula(full_formula_str), 
#                                             data = gam_data, family="binomial", ...)
#                      } else {
#                        gam_data$y_var_numeric <- as.numeric(gam_data$y_var)
#                        full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
#                        gam_fit <- mgcv::gam(as.formula(full_formula_str), 
#                                             data = gam_data, ...)
#                      }
#                      deviance_explained[i] <- summary(gam_fit)$dev.expl
#                    }
                   
#                    # ... (v3의 나머지 로직 동일) ...
#                    weights_magnitude <- deviance_explained
#                    top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0))])
#                    # ... (방향성 보정) ...
#                    list(...)
#                  },
                 
#                  pca_loadings = {
#                    # [V4] 이 블록도 수정할 필요가 없습니다.
#                    # control_vars가 있었다면 expr_mat이 이미 보정되었습니다.
                   
#                    # ... (v3의 pca_loadings 로직 동일) ...
#                    list(...)
#                  }
#   )
  
#   # ===
#   # 5. Return results (v3와 동일)
#   # ===
  
#   result$method <- method
#   result$target_var <- target_var
#   result$n_groups <- n_groups
#   result$n_cells <- ncol(expr_mat)
#   result$formula <- paste(deparse(match.call()), collapse = " ")
  
#   class(result) <- c("gene_signature", "list")
#   return(result)
# }

# # [V(4.1)] find_gene_signature_v4.1
# #
# # 변경점 (2025-11-10):
# # 1. [FIX-NA] 2.Process target: meta.data 서브셋 후 'droplevels()'를 호출하여
# #    control_vars/target_binary의 미사용 레벨을 제거 (model.matrix NA 에러 방지)
# # 2. [FIX-Lasso] 4.Method (lasso): 'nonzero_genes' 필터링 로직 제거.
# #    항상 'abs(weights_all)' 기준 상위 n_features 반환 (n_features 개수 보장)
# # 3. [FIX-NMF] 3.8 Pre-correction: 'nmf'를 'removeBatchEffect' 대상에서 제외.
# # 4. [FIX-NMF] 4.Method (nmf): 'nmf' 블록 *내부*에서 control_vars를
# #    'removeBatchEffect'로 보정하고, 즉시 'min-shift'로 양수화.
# # 5. [FIX-GAM] 4.Method (gam): 'gam' 블록 내부의 중복된 유전자 필터링
# #    ('test_genes_idx') 로직 제거. 'test_n'으로 필터링된 모든 유전자 사용.
# # 6. [FIX-Wilcox] 4.Method (wilcoxon): 'wilcox.test'/'kruskal.test'에서 '...' 인자 제거.
# #' @export
# find_gene_signature_v4.1 <- function(data, 
#                                  meta.data = NULL,
#                                  target_var,
#                                  target_group = NULL,
#                                  control_vars = NULL,   
#                                  method = c("tree_based", "lasso", "limma", 
#                                             "nmf", "wilcoxon", "gam", "pca_loadings"),
#                                  n_features = 50,
#                                  test_n = NULL,         
#                                  preprocess = TRUE,
#                                  min_cells = 10,
#                                  min_pct = 0.01,
#                                  return_model = FALSE,
#                                  fgs_seed = 42,
#                                  lambda_selection = "lambda.1se",
#                                  ...) {
  
#   # ===
#   # 0. Batch Processing (v4와 동일)
#   # ===
  
#   all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
#   if (is.null(method)) {
#     method <- all_methods
#   }
  
#   if (length(method) > 1) {
#     current_call <- match.call()
    
#     results_list <- lapply(method, function(m) {
#       single_call <- current_call
#       single_call$method <- m
#       tryCatch({
#         eval(single_call)
#       }, error = function(e) {
#         warning(sprintf("Method '%s' failed with error: %s", m, e$message))
#         return(list(method = m, error = e$message))
#       })
#     })
    
#     names(results_list) <- method
#     return(results_list)
#   }
  
#   if (!method %in% all_methods) {
#     stop(sprintf("Invalid method '%s'. Choose from: %s", 
#                  method, paste(all_methods, collapse=", ")))
#   }
  
#   set.seed(fgs_seed)
  
#   # ===
#   # 1. Input validation and data extraction (v4와 동일)
#   # ===
  
#   is_seurat <- inherits(data, "Seurat")
  
#   if (is_seurat) {
#     if (!requireNamespace("Seurat", quietly = TRUE)) {
#       stop("Seurat package required but not installed")
#     }
#     if (is.null(meta.data)) {
#       meta.data <- data@meta.data
#     }
#     # Seurat v5/v4 호환
#     if ("data" %in% Seurat::Layers(data)) {
#         expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
#     } else if ("counts" %in% Seurat::Layers(data)) {
#         expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
#     } else if ("data" %in% names(data@assays[[Seurat::DefaultAssay(data)]])) {
#         expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
#     } else {
#         expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
#     }
#   } else {
#     if (is.null(meta.data)) {
#       stop("meta.data must be provided when data is not a Seurat object")
#     }
#     expr_mat <- as.matrix(data)
#     if (nrow(expr_mat) == nrow(meta.data)) {
#       expr_mat <- t(expr_mat)
#     }
#   }
  
#   if (!target_var %in% colnames(meta.data)) {
#     stop(sprintf("target_var '%s' not found in metadata columns: %s", 
#                  target_var, paste(colnames(meta.data), collapse=", ")))
#   }
  
#   if (!is.null(control_vars)) {
#     missing_vars <- control_vars[!control_vars %in% colnames(meta.data)]
#     if (length(missing_vars) > 0) {
#       stop(sprintf("control_vars not found in metadata: %s", 
#                    paste(missing_vars, collapse=", ")))
#     }
#   }
  
#   common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
#   if (length(common_cells) == 0) {
#     stop("No common cells found between expression matrix and metadata")
#   }
#   expr_mat <- expr_mat[, common_cells]
#   meta.data <- meta.data[common_cells, ]
  
#   # ===
#   # 2. Process target variable (v4.1 수정)
#   # ===
  
#   target_values <- meta.data[[target_var]]
#   target_type <- class(target_values)[1]
  
#   if (is.numeric(target_values)) {
#     if (is.null(target_group)) {
#       target_group <- 0.25
#     }
#     if (is.list(target_group)) {
#       low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
#       high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
#       group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                              ifelse(target_values >= high_cutoff, "High", NA))
#     } else if (length(target_group) == 1 && target_group < 1) {
#       low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
#       high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
#       group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                              ifelse(target_values >= high_cutoff, "High", NA))
#     } else {
#       group_labels <- ifelse(target_values < target_group, "Low", "High")
#     }
#     keep_cells <- !is.na(group_labels)
#   } else {
#     if (!is.null(target_group)) {
#       keep_cells <- target_values %in% target_group
#     } else {
#       keep_cells <- rep(TRUE, length(target_values))
#     }
#     group_labels <- target_values
#     keep_cells <- keep_cells & !is.na(group_labels)
#   }
  
#   expr_mat <- expr_mat[, keep_cells]
#   meta.data <- meta.data[keep_cells, ]
  
#   # [FIX-NA] 사용하지 않는 팩터 레벨 제거 (model.matrix 에러 방지)
#   target_binary <- factor(group_labels[keep_cells])
  
#   if (!is.null(control_vars)) {
#     for (cv in control_vars) {
#       if (is.factor(meta.data[[cv]])) {
#         meta.data[[cv]] <- droplevels(meta.data[[cv]])
#       }
#     }
#   }

#   if (length(unique(target_binary)) < 2) {
#     stop("Target variable must have at least 2 groups after processing")
#   }
#   n_groups <- length(unique(target_binary))
  
#   # ===
#   # 3. Filter and preprocess genes (v4와 동일)
#   # ===
  
#   n_cells_expr <- rowSums(expr_mat > 0)
#   pct_cells_expr <- n_cells_expr / ncol(expr_mat)
#   keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
#   expr_mat <- expr_mat[keep_genes, ]
  
#   if (nrow(expr_mat) == 0) {
#     stop("No genes pass filtering criteria")
#   }
  
#   # ===
#   # 3.5 [V4] Speed-up: Pre-filter genes using test_n
#   # ===
#   if (!is.null(test_n) && nrow(expr_mat) > test_n) {
#     message(sprintf("Pre-filtering to top %d genes based on limma p-value...", test_n))
#     if (!requireNamespace("limma", quietly = TRUE)) {
#       stop("limma package required for 'test_n' pre-filtering")
#     }
    
#     if (is.null(control_vars)) {
#       design_test <- model.matrix(~ target_binary, data = meta.data)
#     } else {
#       formula_test <- as.formula(paste("~ target_binary +", paste(control_vars, collapse="+")))
#       design_test <- model.matrix(formula_test, data = meta.data)
#     }
    
#     fit_test <- limma::lmFit(expr_mat, design_test)
#     fit_test <- limma::eBayes(fit_test)
    
#     coef_indices <- grep("target_binary", colnames(design_test))
#     if (length(coef_indices) == 0) {
#         warning("Could not find target_binary coefficients for test_n pre-filtering.")
#         # Fallback: use all coefficients (less ideal)
#         coef_indices <- 1:ncol(design_test)
#     }
    
#     top_table_test <- limma::topTable(fit_test, 
#                                       coef = coef_indices, 
#                                       number = Inf, sort.by = "P")
    
#     top_gene_names <- rownames(top_table_test)[1:min(test_n, nrow(top_table_test))]
#     expr_mat <- expr_mat[top_gene_names, ]
#     message(sprintf("... reduced to %d genes.", nrow(expr_mat)))
#   }
  
#   if (preprocess) {
#     if (max(expr_mat) > 100) {
#       expr_mat <- log1p(expr_mat)
#     }
#     if (method %in% c("lasso", "gam", "pca_loadings")) {
#       gene_means <- rowMeans(expr_mat)
#       gene_sds <- apply(expr_mat, 1, sd)
#       gene_sds[gene_sds == 0] <- 1
#       expr_mat <- (expr_mat - gene_means) / gene_sds
#     }
#   }
  
#   # ===
#   # 3.8 [V4.1] Confounder pre-correction (non-native methods)
#   # ===
  
#   # [FIX-NMF] 'nmf'는 removeBatchEffect 대상에서 제외 (음수 값 방지)
#   if (!is.null(control_vars) && method %in% c("wilcoxon", "pca_loadings")) { 
#     if (!requireNamespace("limma", quietly = TRUE)) {
#       stop("limma package required for removeBatchEffect")
#     }
    
#     message(sprintf("Method %s: Applying limma::removeBatchEffect for: %s",
#                     method, paste(control_vars, collapse=", ")))
    
#     covariates_df <- meta.data[, control_vars, drop = FALSE]
    
#     # Check for factors and create model matrix
#     if(any(sapply(covariates_df, is.factor)) || any(sapply(covariates_df, is.character))) {
#         covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
#     } else {
#         covariate_mat <- as.matrix(covariates_df)
#     }
    
#     expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
#   }
  
  
#   # ===
#   # 4. Method-specific signature discovery (v4.1 수정)
#   # ===
  
#   X <- t(expr_mat)
#   y <- target_binary
  
#   covariate_mat_model <- NULL # Initialize
#   if (!is.null(control_vars) && method %in% c("tree_based", "lasso", "gam")) {
#      covariates_df_model <- meta.data[, control_vars, drop = FALSE]
#      if(any(sapply(covariates_df_model, is.factor)) || any(sapply(covariates_df_model, is.character))) {
#         covariate_mat_model <- model.matrix(~ . - 1, data = covariates_df_model)
#      } else {
#         covariate_mat_model <- as.matrix(covariates_df_model)
#      }
#   }

  
#   result <- switch(method,
                 
#                  tree_based = {
#                    if (!requireNamespace("randomForest", quietly = TRUE)) {
#                      stop("randomForest package required. Install with: install.packages('randomForest')")
#                    }
                   
#                    X_rf <- X
                   
#                    if (!is.null(control_vars)) {
#                      X_rf <- cbind(X_rf, covariate_mat_model)
#                    }
                   
#                    # Pre-filter genes if matrix is still too large
#                    if (ncol(X_rf) > 2000) {
#                      gene_vars <- apply(X, 2, var) # var from original X
#                      top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
                     
#                      control_var_names <- if(is.null(control_vars)) character(0) else colnames(covariate_mat_model)
#                      X_rf <- X_rf[, c(top_var_genes, control_var_names)]
#                    }
                   
#                    rf_model <- randomForest::randomForest(x = X_rf, y = y, ntree = 500, importance = TRUE, ...)
                   
#                    importance_scores <- randomForest::importance(rf_model)
#                    if (n_groups == 2) {
#                      weights_magnitude_all <- importance_scores[, "MeanDecreaseGini"]
#                    } else {
#                      weights_magnitude_all <- rowMeans(importance_scores[, grep("MeanDecreaseGini", colnames(importance_scores))])
#                    }
                   
#                    gene_names_in_model <- colnames(X_rf)[!colnames(X_rf) %in% colnames(covariate_mat_model)]
#                    weights_magnitude_genes <- weights_magnitude_all[gene_names_in_model]
                   
#                    top_genes <- names(sort(weights_magnitude_genes, decreasing=TRUE)[1:min(n_features, length(weights_magnitude_genes))])
#                    weights_magnitude <- weights_magnitude_genes[top_genes]
                   
#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]
#                      g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      effect_size <- mean_g2 - mean_g1 
#                      weights <- weights_magnitude * sign(effect_size)
#                    } else {
#                      warning("tree_based: n_groups > 2. Score represents magnitude (importance), not direction.")
#                      weights <- weights_magnitude
#                    }
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    pred <- rf_model$predicted
#                    if (n_groups == 2) {
#                      if (requireNamespace("pROC", quietly = TRUE)) {
#                        roc_obj <- pROC::roc(y, scores, quiet=TRUE)
#                        auc <- as.numeric(pROC::auc(roc_obj))
#                      } else { auc <- NA }
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
#                    } else {
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, confusion = table(pred, y))
#                    }
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) rf_model else NULL)
#                  },
                 
#                  lasso = {
#                    if (!requireNamespace("glmnet", quietly = TRUE)) {
#                      stop("glmnet package required. Install with: install.packages('glmnet')")
#                    }
                   
#                    if (is.null(control_vars)) {
#                      X_model <- X
#                      penalty_vec <- rep(1, ncol(X_model))
#                    } else {
#                      X_model <- cbind(X, covariate_mat_model)
#                      penalty_vec <- c(rep(1, ncol(X)), rep(0, ncol(covariate_mat_model)))
#                    }
                   
#                    if (n_groups == 2) {
#                      cv_fit <- glmnet::cv.glmnet(X_model, y, family="binomial", alpha=1, 
#                                                  penalty.factor = penalty_vec, ...)
#                    } else {
#                      cv_fit <- glmnet::cv.glmnet(X_model, y, family="multinomial", alpha=1, 
#                                                  penalty.factor = penalty_vec, ...)
#                    }
                   
#                    coefs <- coef(cv_fit, s = lambda_selection) 
                   
#                    if (n_groups == 2) {
#                      weights_all <- as.numeric(coefs[2:(ncol(X)+1)]) 
#                      names(weights_all) <- rownames(coefs)[2:(ncol(X)+1)]
#                    } else {
#                      coef_list <- lapply(coefs, function(x) as.numeric(x[2:(ncol(X)+1)]))
#                      weights_all <- rowMeans(do.call(cbind, coef_list))
#                      names(weights_all) <- rownames(coefs[[1]])[2:(ncol(X)+1)]
#                    }
                   
#                    if (length(weights_all) == 0) {
#                      stop("LASSO returned no gene coefficients.")
#                    }

#                    # [FIX-Lasso] 'nonzero' 필터링 제거. 항상 abs(weight)로 정렬.
#                    top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
                   
#                    weights <- weights_all[top_genes]
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    pred_probs <- predict(cv_fit, newx=X_model, s = lambda_selection, type="response")
#                    if (n_groups == 2) {
#                      pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
#                      if (requireNamespace("pROC", quietly = TRUE)) {
#                        roc_obj <- pROC::roc(y, scores, quiet=TRUE) 
#                        auc <- as.numeric(pROC::auc(roc_obj))
#                      } else { auc <- NA }
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
#                    } else {
#                      pred_class_indices <- apply(pred_probs[,,1], 1, which.max)
#                      pred <- levels(y)[pred_class_indices]
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, confusion = table(pred, y))
#                    }
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) cv_fit else NULL)
#                  },
                 
#                  limma = {
#                    if (!requireNamespace("limma", quietly = TRUE)) {
#                      stop("limma package required. Install with: BiocManager::install('limma')")
#                    }
                   
#                    if (is.null(control_vars)) {
#                      design <- model.matrix(~0 + target_binary, data = meta.data)
#                      colnames(design)[1:n_groups] <- levels(target_binary)
#                    } else {
#                      control_formula <- paste(control_vars, collapse = " + ")
#                      full_formula <- as.formula(paste("~0 + target_binary +", control_formula))
#                      design <- model.matrix(full_formula, data = meta.data)
#                      colnames(design)[1:n_groups] <- levels(target_binary) 
#                    }

#                    fit <- limma::lmFit(expr_mat, design)
                   
#                    if (n_groups == 2) {
#                      contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
#                    } else {
#                      # Create all pairwise contrasts if n_groups > 2
#                      contrast_pairs <- combn(levels(target_binary), 2, function(x) paste(x[2], x[1], sep="-"))
#                      contrast_str <- contrast_pairs
#                      warning("limma: n_groups > 2. Using all pairwise contrasts. Weights based on average t-statistic.")
#                    }
                   
#                    contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
#                    fit2 <- limma::contrasts.fit(fit, contrast_mat)
#                    fit2 <- limma::eBayes(fit2)
                   
#                    if(n_groups == 2) {
#                        top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
#                        weights_all <- top_table$t
#                    } else {
#                        # For multiple contrasts, sort by F-statistic
#                        top_table <- limma::topTable(fit2, number=Inf, sort.by="F")
#                        # Use average absolute t-stat as weight magnitude
#                        weights_all <- rowMeans(abs(fit2$t[, contrast_str, drop=FALSE]))
#                    }
                   
#                    names(weights_all) <- rownames(top_table)
                   
#                    top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
                   
#                    # Re-fetch signed weights for top genes
#                    if(n_groups == 2) {
#                        weights <- top_table[top_genes, "t"]
#                    } else {
#                        # For n_groups > 2, weights are just magnitude (no clear direction)
#                        weights <- weights_all[top_genes]
#                    }

#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(top_table = top_table[top_genes, ])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) fit2 else NULL)
#                  },
                 
#                  wilcoxon = {
#                    # [V4.1] control_vars가 있었다면 3.8에서 expr_mat이 이미 보정됨.
                   
#                    pvals <- numeric(nrow(expr_mat))
#                    effect_sizes <- numeric(nrow(expr_mat))
#                    names(pvals) <- rownames(expr_mat)
#                    names(effect_sizes) <- rownames(expr_mat)

#                    for (i in 1:nrow(expr_mat)) {
#                      if (n_groups == 2) {
#                        group1 <- expr_mat[i, target_binary == levels(target_binary)[1]]
#                        group2 <- expr_mat[i, target_binary == levels(target_binary)[2]]
                       
#                        # [FIX-Wilcox] '...' 인자 제거
#                        test <- try(wilcox.test(group1, group2), silent=TRUE) 
                       
#                        if(inherits(test, "try-error")) {
#                          pvals[i] <- 1.0
#                          effect_sizes[i] <- 0
#                        } else {
#                          pvals[i] <- test$p.value
#                          effect_sizes[i] <- median(group2) - median(group1)
#                        }
#                      } else {
#                        # [FIX-Wilcox] '...' 인자 제거
#                        test <- try(kruskal.test(expr_mat[i, ] ~ target_binary), silent=TRUE) 
                       
#                        if(inherits(test, "try-error")) {
#                           pvals[i] <- 1.0
#                           effect_sizes[i] <- 0
#                        } else {
#                           pvals[i] <- test$p.value
#                           # Use variance of medians as effect size (magnitude only)
#                           effect_sizes[i] <- var(tapply(expr_mat[i, ], target_binary, median))
#                        }
#                      }
#                    }
                   
#                    padj <- p.adjust(pvals, method="BH")
                   
#                    ranks <- rank(pvals) + rank(-abs(effect_sizes))
#                    top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
                   
#                    weights <- effect_sizes[top_genes]
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(pvalues = pvals[top_genes], padj = padj[top_genes],
#                                 effect_sizes = effect_sizes[top_genes])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = NULL)
#                  },
                 
#                  nmf = {
#                    if (!requireNamespace("NMF", quietly = TRUE)) {
#                      stop("NMF package required. Install with: install.packages('NMF')")
#                    }
                   
#                    # [FIX-NMF] 'nmf' 블록 내부에서 보정 실행
#                    if (!is.null(control_vars)) {
#                      if (!requireNamespace("limma", quietly = TRUE)) {
#                        stop("limma package required for NMF confounder correction")
#                      }
#                      warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s. Re-shifting to non-negative.",
#                                      paste(control_vars, collapse=", ")))
                     
#                      covariates_df <- meta.data[, control_vars, drop = FALSE]
#                      if(any(sapply(covariates_df, is.factor)) || any(sapply(covariates_df, is.character))) {
#                         covariate_mat <- model.matrix(~ . - 1, data = covariates_df)
#                      } else {
#                         covariate_mat <- as.matrix(covariates_df)
#                      }
                     
#                      expr_mat <- limma::removeBatchEffect(expr_mat, covariates = covariate_mat)
#                    }
                   
#                    # [FIX-NMF] 보정 후 (혹은 보정 없이) 무조건 양수로 이동
#                    expr_mat_pos <- expr_mat - min(expr_mat) + 0.01 
                   
#                    rank <- min(n_groups + 2, 10)
#                    dots <- list(...)
#                    dots$seed <- NULL 
                   
#                    nmf_res <- do.call(NMF::nmf, c(list(x = expr_mat_pos, rank = rank, seed=fgs_seed), dots))
#                    W <- NMF::basis(nmf_res)
#                    H <- NMF::coef(nmf_res)
                   
#                    component_cors <- numeric(rank)
#                    for (k in 1:rank) {
#                      if (n_groups == 2) {
#                        component_cors[k] <- abs(cor(H[k, ], as.numeric(target_binary)))
#                      } else {
#                        component_cors[k] <- summary(aov(H[k, ] ~ target_binary))[[1]][1, "F value"]
#                      }
#                    }
                   
#                    best_component <- which.max(component_cors)
#                    weights_magnitude <- W[, best_component]
#                    names(weights_magnitude) <- rownames(expr_mat_pos)
                   
#                    top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, length(weights_magnitude))])
#                    weights_magnitude <- weights_magnitude[top_genes]
                   
#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]
#                      g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      effect_size <- mean_g2 - mean_g1
#                      weights <- weights_magnitude * sign(effect_size)
#                    } else {
#                      warning("nmf: n_groups > 2. Score represents magnitude (importance), not direction.")
#                      weights <- weights_magnitude
#                    }
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(component = best_component, correlation = component_cors[best_component])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) nmf_res else NULL)
#                  },
                 
#                  gam = {
#                    if (!requireNamespace("mgcv", quietly = TRUE)) {
#                      stop("mgcv package required. Install with: install.packages('mgcv')")
#                    }
                   
#                    deviance_explained <- numeric(nrow(expr_mat))
#                    names(deviance_explained) <- rownames(expr_mat)
                   
#                    if (is.null(control_vars)) {
#                      formula_base <- ""
#                      gam_data <- data.frame(y_var = target_binary)
#                    } else {
#                      # Ensure control vars are in the data frame
#                      gam_data <- data.frame(y_var = target_binary, 
#                                             meta.data[, control_vars, drop=FALSE])
#                      # Create formula string *from names in gam_data*
#                      formula_base <- paste(" +", paste(colnames(covariate_mat_model), collapse=" + "))
#                      # Add the covariates matrix to gam_data
#                      gam_data <- cbind(gam_data, covariate_mat_model)
#                    }
                   
#                    # [FIX-GAM] Use all genes remaining in expr_mat (1:nrow)
#                    for (i in 1:nrow(expr_mat)) {
#                      gam_data$gene_expr <- expr_mat[i, ]
                     
#                      if (n_groups == 2) {
#                        gam_data$y_var_numeric <- as.numeric(gam_data$y_var) - 1
#                        full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
#                        gam_fit <- try(mgcv::gam(as.formula(full_formula_str), 
#                                             data = gam_data, family="binomial", ...), silent=TRUE)
#                      } else {
#                        gam_data$y_var_numeric <- as.numeric(gam_data$y_var)
#                        full_formula_str <- paste("y_var_numeric ~ s(gene_expr)", formula_base)
#                        gam_fit <- try(mgcv::gam(as.formula(full_formula_str), 
#                                             data = gam_data, ...), silent=TRUE)
#                      }
                     
#                      if(inherits(gam_fit, "try-error")) {
#                         deviance_explained[i] <- 0
#                      } else {
#                         deviance_explained[i] <- summary(gam_fit)$dev.expl
#                      }
#                    }
                   
#                    weights_magnitude <- deviance_explained
#                    top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0, na.rm=TRUE))])
                   
#                    if (length(top_genes) == 0) {
#                      warning("GAM method found no genes with deviance > 0.")
#                      return(list(genes=character(0), weights=numeric(0), scores=numeric(0), performance=list()))
#                    }
                   
#                    weights_magnitude <- weights_magnitude[top_genes]

#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]
#                      g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      effect_size <- mean_g2 - mean_g1
#                      weights <- weights_magnitude * sign(effect_size)
#                    } else {
#                      warning("gam: n_groups > 2. Score represents magnitude (importance), not direction.")
#                      weights <- weights_magnitude
#                    }
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(deviance_explained = weights_magnitude)
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = NULL)
#                  },
                 
#                  pca_loadings = {
#                    # [V4.1] control_vars가 있었다면 3.8에서 expr_mat이 이미 보정됨.
#                    # PCA는 스케일링된 데이터를 사용해야 함 (3.x에서 처리됨)
                   
#                    pca_res <- prcomp(X, center=FALSE, scale.=FALSE) 
                   
#                    n_pcs_to_test <- min(50, ncol(pca_res$x))
#                    pc_cors <- numeric(n_pcs_to_test)
                   
#                    for (k in 1:n_pcs_to_test) {
#                      if (n_groups == 2) {
#                        pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
#                      } else {
#                        aov_res <- try(summary(aov(pca_res$x[, k] ~ target_binary)), silent=TRUE)
#                        if(inherits(aov_res, "try-error")) {
#                            pc_cors[k] <- 0
#                        } else {
#                            pc_cors[k] <- aov_res[[1]][1, "F value"]
#                        }
#                      }
#                    }
                   
#                    best_pc <- which.max(pc_cors)
#                    weights_all <- pca_res$rotation[, best_pc]
                   
#                    top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
#                    weights <- weights_all[top_genes]
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(PC = best_pc, correlation = pc_cors[best_pc],
#                                 variance_explained = summary(pca_res)$importance[2, best_pc])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) pca_res else NULL)
#                  }
#   )
  
#   # ===
#   # 5. Return results (v4와 동일)
#   # ===
  
#   result$method <- method
#   result$target_var <- target_var
#   result$n_groups <- n_groups
#   result$n_cells <- ncol(expr_mat)
#   result$formula <- paste(deparse(match.call()), collapse = " ")
  
#   class(result) <- c("gene_signature", "list")
#   return(result)
# }

#' @export
find_gene_signature_v4.2_test_to_delete <- function(data, 
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
fgs_preprocess_data_test_to_delete <- function(data, 
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
find_gene_signature_v5_test_to_delete <- function(data, 
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
fgs_preprocess_data_v5.2_test_to_delete <- function(data, 
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
find_gene_signature_v5.2_test_to_delet <- function(data, 
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


# #' Train a Meta-Learner (L2 Stacking Model)
# #'
# #' Trains multiple L2 candidate models on L1 scores and selects the best
# #' one based on cross-validation performance.
# #'
# #' @param l1_signatures A named list of trained L1 signatures (outputs from
# #'                      find_gene_signature). e.g., list(lasso=sig1, rf=sig2)
# #' @param holdout_data A Seurat object or matrix (genes x cells) NOT used
# #'                      for training L1 models.
# #' @param target_var The target variable column name in holdout_data@meta.data.
# #' @param l2_methods A vector of model methods supported by caret 
# #'                   (e.g., c("glm", "ranger", "xgbTree", "svmRadial")).
# #' @param k_folds Number of folds for cross-validation to select the best L2 model.
# #' @param metric The metric to optimize (e.g., "AUC", "Accuracy"). 
# #'               (Note: For "AUC", target must be binary factor).
# #' @param fgs_seed Seed for reproducibility.
# #'
# #' @return A list containing:
# #'    - $best_model: The final L2 model (a caret 'train' object) 
# #'                     trained on all l2_train_df.
# #'    - $model_comparison: Results from caret::resamples comparing L2 candidates.
# #'    - $l2_train_df: The dataframe used for L2 training (L1 scores + target).
# #'    - $l1_signatures: The provided L1 signatures (for reference).
# #' @export
# train_meta_learner_v1 <- function(l1_signatures,
#                                holdout_data,
#                                target_var,
#                                l2_methods = c("glm", "ranger", "xgbTree"),
#                                k_folds = 5,
#                                metric = "AUC",
#                                fgs_seed = 42) {
  
#   if (!requireNamespace("caret", quietly = TRUE)) {
#     stop("caret package required. Install with: install.packages('caret')")
#   }
#   set.seed(fgs_seed)

#   # === 1. L2 Feature 생성 ===
#   message("--- Generating L2 features from holdout_data ---")
  
#   l2_features_list <- lapply(l1_signatures, function(sig) {
#     score_signature(expr_data = holdout_data,
#                     signature = sig,
#                     normalize = TRUE) # 스케일 통일을 위해 정규화
#   })
  
#   l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  
#   # L1 모델 이름이 없으면 기본 이름 할당
#   if (is.null(names(l1_signatures))) {
#     names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures))
#   }
#   colnames(l2_train_df) <- names(l1_signatures)

#   # === 2. L2 훈련 데이터 준비 ===
  
#   # Seurat 객체에서 메타데이터 추출
#   if (inherits(holdout_data, "Seurat")) {
#     meta.data <- holdout_data@meta.data
#     # holdout_data의 셀 순서와 l2_train_df의 행 순서가 일치하는지 확인
#     if(!all(rownames(l2_train_df) == rownames(meta.data))) {
#         meta.data <- meta.data[rownames(l2_train_df), ]
#     }
#     l2_target <- meta.data[[target_var]]
#   } else {
#     # 행렬인 경우, target_var는 별도의 벡터로 제공되어야 함 (이 예제에서는 단순화)
#     stop("holdout_data is not a Seurat object. Please provide target_var as a separate vector matching cell order.")
#   }

#   # caret은 factor형 타겟 변수를 선호 (특히 분류 및 AUC 계산 시)
#   if (!is.factor(l2_target)) {
#     l2_target <- factor(l2_target)
#   }
  
#   # caret이 R 변수명으로 부적합한 이름을 처리하도록 함 (예: "limma-fold")
#   colnames(l2_train_df) <- make.names(colnames(l2_train_df))
  
#   # === 3. L2 모델 교차 검증 및 선택 ===
#   message("--- Training and comparing L2 models via k-fold CV ---")

#   # CV 설정
#   train_control <- caret::trainControl(
#     method = "cv",
#     number = k_folds,
#     summaryFunction = caret::twoClassSummary, # AUC, Sens, Spec 계산
#     classProbs = TRUE, # AUC 계산을 위해 확률 예측 필요
#     savePredictions = "final",
#     allowParallel = TRUE
#   )
  
#   # metric이 AUC가 아니면 기본 summary 사용
#   if (metric != "AUC") {
#     train_control$summaryFunction <- caret::defaultSummary
#   }
  
#   # 여러 L2 모델을 훈련시키기 위한 리스트
#   model_list <- list()
  
#   for (method in l2_methods) {
#     tryCatch({
#       message(sprintf("Training L2 candidate: %s", method))
#       model_list[[method]] <- caret::train(
#         x = l2_train_df,
#         y = l2_target,
#         method = method,
#         trControl = train_control,
#         metric = metric,
#         tuneLength = 5 # 각 모델의 기본 하이퍼파라미터 튜닝
#       )
#     }, error = function(e) {
#       warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
#     })
#   }

#   if (length(model_list) == 0) {
#     stop("No L2 models were successfully trained.")
#   }

#   # 훈련된 모델들의 성능 비교
#   model_comparison <- caret::resamples(model_list)
#   print(summary(model_comparison))

#   # === 4. Best 모델 선택 및 재훈련 (caret이 자동으로 처리) ===
#   # caret::train은 CV를 사용해 최적의 하이퍼파라미터를 찾고,
#   # 그 파라미터로 "전체 L2 데이터"를 재훈련시킨 모델($finalModel)을 반환함.
#   # 우리는 이 훈련된 모델들 중 'Best'를 선택하기만 하면 됨.
  
#   best_model_name <- names(which.max(
#     sapply(model_list, function(m) max(m$results[, metric]))
#   ))
  
#   best_model <- model_list[[best_model_name]]
  
#   message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
#                   best_model_name, metric, max(best_model$results[, metric])))
  
#   # === 5. 결과 반환 ===
#   return(list(
#     best_model = best_model,
#     best_model_name = best_model_name,
#     model_comparison = model_comparison,
#     l2_train_df = data.frame(l2_train_df, target = l2_target),
#     l1_signatures = l1_signatures
#   ))
# }


# # ---
# # 2. train_meta_learner (수정본)
# #    - (FIX 1): GetAssayData 1회 호출
# #    - (FIX 2): NA target 명시적 제거
# # ---
# # (caret, ranger, xgboost 등 필요 패키지 로드 가정)
# #' Train a Meta-Learner (L2 Stacking Model)
# #'
# #' Trains multiple L2 candidate models on L1 scores and selects the best
# #' one based on cross-validation performance.
# #'
# #' @param l1_signatures A named list of trained L1 signatures (outputs from
# #'                      find_gene_signature). e.g., list(lasso=sig1, rf=sig2)
# #' @param holdout_data A Seurat object or matrix (genes x cells) NOT used
# #'                      for training L1 models.
# #' @param target_var The target variable column name in holdout_data@meta.data.
# #' @param l2_methods A vector of model methods supported by caret 
# #'                   (e.g., c("glm", "ranger", "xgbTree", "svmRadial")).
# #' @param k_folds Number of folds for cross-validation to select the best L2 model.
# #' @param metric The metric to optimize (e.g., "AUC", "Accuracy"). 
# #'               (Note: For "AUC", target must be binary factor).
# #' @param fgs_seed Seed for reproducibility.
# #'
# #' @return A list containing:
# #'    - $best_model: The final L2 model (a caret 'train' object) 
# #'                     trained on all l2_train_df.
# #'    - $model_comparison: Results from caret::resamples comparing L2 candidates.
# #'    - $l2_train_df: The dataframe used for L2 training (L1 scores + target).
# #'    - $l1_signatures: The provided L1 signatures (for reference).
# #' @export
# train_meta_learner_v2 <- function(l1_signatures,
#                                holdout_data,
#                                target_var,
#                                l2_methods = c("glm", "ranger", "xgbTree"),
#                                k_folds = 5,
#                                metric = "AUC",
#                                fgs_seed = 42) {
  
#   if (!requireNamespace("caret", quietly = TRUE)) {
#     stop("caret package required.")
#   }
#   set.seed(fgs_seed)

#   # === 1. 데이터 추출 (PERFORMANCE FIX 1) ===
#   message("--- Preparing holdout data (Extracting once) ---")
#   if (inherits(holdout_data, "Seurat")) {
#     # GetAssayData를 루프 밖에서 단 한 번만 호출 (희소 행렬로)
#     expr_mat <- Seurat::GetAssayData(holdout_data, layer="data")
#     meta_data <- holdout_data@meta.data
#   } else {
#     # Seurat 객체가 아니면, holdout_data가 이미 (희소) 행렬이라고 가정
#     expr_mat <- holdout_data
#     # 이 경우 target_var는 반드시 메타데이터 프레임이거나 벡터여야 함 (로직 단순화)
#     if (!is.vector(target_var) && !is.factor(target_var)) {
#         stop("If holdout_data is a matrix, target_var must be the actual target vector.")
#     }
#     meta_data <- NULL # 메타데이터 객체가 없음
#   }

#   # === 2. L2 Feature 생성 ===
#   message("--- Generating L2 features from holdout data ---")
  
#   l2_features_list <- lapply(l1_signatures, function(sig) {
#     # Seurat 객체 대신 추출된 'expr_mat'를 전달
#     score_signature(expr_data = expr_mat, 
#                     signature = sig,
#                     normalize = TRUE) 
#   })
  
#   l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  
#   if (is.null(names(l1_signatures))) {
#     names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures))
#   }
#   # caret이 인식하도록 R 변수명으로 정리
#   colnames(l2_train_df) <- make.names(names(l1_signatures))

#   # === 3. L2 타겟 변수 준비 ===
#   if (!is.null(meta_data)) {
#     # l2_train_df의 행 순서(cell ID)와 meta_data의 행 순서가 일치하는지 확인/정렬
#     if (!all(rownames(l2_train_df) == rownames(meta_data))) {
#         meta_data <- meta_data[rownames(l2_train_df), , drop = FALSE]
#     }
#     l2_target <- meta_data[[target_var]]
#   } else {
#     l2_target <- target_var # 이미 벡터라고 가정
#     if (length(l2_target) != nrow(l2_train_df)) {
#       stop("Target vector length does not match expression data cell count.")
#     }
#   }

#   if (!is.factor(l2_target)) {
#     l2_target <- factor(l2_target)
#   }
  
#   # === 4. NA 타겟 제거 (BUG FIX 2) ===
#   message(sprintf("Initial L2 target length: %d", length(l2_target)))
#   na_indices <- is.na(l2_target)
  
#   if (any(na_indices)) {
#     n_na <- sum(na_indices)
#     message(sprintf("Found and removing %d cells with NA target.", n_na))
    
#     l2_target_filtered <- l2_target[!na_indices]
#     l2_train_df_filtered <- l2_train_df[!na_indices, , drop = FALSE]
    
#   } else {
#     l2_target_filtered <- l2_target
#     l2_train_df_filtered <- l2_train_df
#   }
  
#   message(sprintf("Final L2 training set size: %d cells", length(l2_target_filtered)))

#   if (length(l2_target_filtered) == 0 || nrow(l2_train_df_filtered) == 0) {
#       stop("No valid data remaining after removing NAs from target variable.")
#   }

#   # === 5. L2 모델 교차 검증 및 선택 ===
#   message("--- Training and comparing L2 models via k-fold CV ---")

#   # 타겟 유형(이진/다중)에 따라 CV 컨트롤 분기
#   if (length(levels(l2_target_filtered)) == 2) {
#     message("Binary target detected. Using twoClassSummary (for AUC).")
#     train_control <- caret::trainControl(
#       method = "cv", number = k_folds,
#       summaryFunction = caret::twoClassSummary,
#       classProbs = TRUE,
#       savePredictions = "final",
#       allowParallel = TRUE
#     )
#     # 이진 분류 시 metric이 AUC가 아니면 AUC로 강제
#     if (metric != "AUC") {
#       message("NOTE: Target is binary, forcing metric to 'AUC'.")
#       metric <- "AUC"
#     }
#   } else {
#     message("Multi-class target detected. Using defaultSummary.")
#     train_control <- caret::trainControl(
#       method = "cv", number = k_folds,
#       summaryFunction = caret::defaultSummary,
#       savePredictions = "final",
#       allowParallel = TRUE
#     )
#     # 다중 분류 시 AUC 사용 불가
#     if (metric == "AUC") {
#       message("WARNING: Target is multi-class, 'AUC' metric is invalid. Defaulting to 'Accuracy'.")
#       metric <- "Accuracy"
#     }
#   }
  
#   model_list <- list()
  
#   for (method in l2_methods) {
#     tryCatch({
#       message(sprintf("Training L2 candidate: %s", method))
#       model_list[[method]] <- caret::train(
#         x = l2_train_df_filtered, # <-- 필터링된 데이터 사용
#         y = l2_target_filtered,   # <-- 필터링된 타겟 사용
#         method = method,
#         trControl = train_control,
#         metric = metric,
#         tuneLength = 5 # 기본 하이퍼파라미터 튜닝
#       )
#     }, error = function(e) {
#       warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
#     })
#   }

#   if (length(model_list) == 0) {
#     stop("No L2 models were successfully trained.")
#   }

#   # 훈련된 모델들의 성능 비교
#   model_comparison <- caret::resamples(model_list)
#   print(summary(model_comparison))

#   # Best 모델 선택
#   best_metric_values <- sapply(model_list, function(m) max(m$results[, metric], na.rm = TRUE))
#   best_model_name <- names(which.max(best_metric_values))
#   best_model <- model_list[[best_model_name]]
  
#   message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
#                   best_model_name, metric, max(best_model$results[, metric], na.rm = TRUE)))
  
#   # === 6. 결과 반환 ===
#   return(list(
#     best_model = best_model,
#     best_model_name = best_model_name,
#     model_comparison = model_comparison,
#     l2_train_df_filtered = data.frame(l2_train_df_filtered, target = l2_target_filtered),
#     l1_signatures = l1_signatures
#   ))
# }

# #' Train a Meta-Learner (L2 Stacking Model)
# #'
# #' Trains multiple L2 candidate models on L1 scores and selects the best
# #' one based on cross-validation performance.
# #'
# #' @param l1_signatures A named list of trained L1 signatures (outputs from
# #'                      find_gene_signature). e.g., list(lasso=sig1, rf=sig2)
# #' @param holdout_data A Seurat object or matrix (genes x cells) NOT used
# #'                      for training L1 models.
# #' @param target_var The target variable column name in holdout_data@meta.data.
# #' @param l2_methods A vector of model methods supported by caret 
# #'                   (e.g., c("glm", "ranger", "xgbTree", "svmRadial")).
# #' @param k_folds Number of folds for cross-validation to select the best L2 model.
# #' @param metric The metric to optimize (e.g., "AUC", "Accuracy"). 
# #'               (Note: For "AUC", target must be binary factor).
# #' @param fgs_seed Seed for reproducibility.
# #'
# #' @return A list containing:
# #'    - $best_model: The final L2 model (a caret 'train' object) 
# #'                     trained on all l2_train_df.
# #'    - $model_comparison: Results from caret::resamples comparing L2 candidates.
# #'    - $l2_train_df: The dataframe used for L2 training (L1 scores + target).
# #'    - $l1_signatures: The provided L1 signatures (for reference).
# #' @export
# train_meta_learner_v3 <- function(l1_signatures,
#                                holdout_data,
#                                target_var,
#                                l2_methods = c("glm", "ranger", "xgbTree"),
#                                k_folds = 5,
#                                metric = "AUC",
#                                fgs_seed = 42) {
  
#   if (!requireNamespace("caret", quietly = TRUE)) {
#     stop("caret package required.")
#   }
#   set.seed(fgs_seed)

#   # === 1. 데이터 추출 (PERFORMANCE FIX 1) ===
#   message("--- Preparing holdout data (Extracting once) ---")
#   if (inherits(holdout_data, "Seurat")) {
#     expr_mat <- Seurat::GetAssayData(holdout_data, layer="data")
#     meta_data <- holdout_data@meta.data
#   } else {
#     expr_mat <- holdout_data
#     if (!is.vector(target_var) && !is.factor(target_var)) {
#         stop("If holdout_data is a matrix, target_var must be the actual target vector.")
#     }
#     meta_data <- NULL 
#   }

#   # === 2. L2 Feature 생성 ===
#   message("--- Generating L2 features from holdout data ---")
  
#   l2_features_list <- lapply(l1_signatures, function(sig) {
#     score_signature(expr_data = expr_mat, 
#                     signature = sig,
#                     normalize = TRUE) 
#   })
  
#   l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  
#   if (is.null(names(l1_signatures))) {
#     names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures))
#   }
#   colnames(l2_train_df) <- make.names(names(l1_signatures))

#   # === 3. L2 타겟 변수 준비 ===
#   if (!is.null(meta_data)) {
#     if (!all(rownames(l2_train_df) == rownames(meta_data))) {
#         meta_data <- meta_data[rownames(l2_train_df), , drop = FALSE]
#     }
#     l2_target <- meta_data[[target_var]]
#   } else {
#     l2_target <- target_var 
#   }

#   if (!is.factor(l2_target)) {
#     l2_target <- factor(l2_target)
#   }
  
#   # === (BUG FIX 3): caret을 위한 Factor Level 이름 소독 ===
#   original_levels <- levels(l2_target)
#   sanitized_levels <- make.names(original_levels)
  
#   if (!all(original_levels == sanitized_levels)) {
#     message(sprintf("Sanitizing target levels for caret: '%s' -> '%s'", 
#                     paste(original_levels, collapse="' / '"), 
#                     paste(sanitized_levels, collapse="' / '")))
#     levels(l2_target) <- sanitized_levels
#   }
#   # === End of fix ===

#   # === 4. NA 타겟 제거 (BUG FIX 2) ===
#   message(sprintf("Initial L2 target length: %d", length(l2_target)))
#   na_indices <- is.na(l2_target)
  
#   if (any(na_indices)) {
#     n_na <- sum(na_indices)
#     message(sprintf("Found and removing %d cells with NA target.", n_na))
    
#     l2_target_filtered <- l2_target[!na_indices]
#     l2_train_df_filtered <- l2_train_df[!na_indices, , drop = FALSE]
    
#   } else {
#     l2_target_filtered <- l2_target
#     l2_train_df_filtered <- l2_train_df
#   }
  
#   message(sprintf("Final L2 training set size: %d cells", length(l2_target_filtered)))

#   if (length(l2_target_filtered) == 0 || nrow(l2_train_df_filtered) == 0) {
#       stop("No valid data remaining after removing NAs from target variable.")
#   }

#   # === 5. L2 모델 교차 검증 및 선택 ===
#   message("--- Training and comparing L2 models via k-fold CV ---")

#   if (length(levels(l2_target_filtered)) == 2) {
#     message("Binary target detected. Using twoClassSummary (for AUC).")
#     train_control <- caret::trainControl(
#       method = "cv", number = k_folds,
#       summaryFunction = caret::twoClassSummary,
#       classProbs = TRUE,
#       savePredictions = "final",
#       allowParallel = TRUE
#     )
#     if (metric != "AUC") {
#       message("NOTE: Target is binary, forcing metric to 'AUC'.")
#       metric <- "AUC"
#     }
#   } else {
#     message("Multi-class target detected. Using defaultSummary.")
#     train_control <- caret::trainControl(
#       method = "cv", number = k_folds,
#       summaryFunction = caret::defaultSummary,
#       savePredictions = "final",
#       allowParallel = TRUE
#     )
#     if (metric == "AUC") {
#       message("WARNING: Target is multi-class, 'AUC' metric is invalid. Defaulting to 'Accuracy'.")
#       metric <- "Accuracy"
#     }
#   }
  
#   model_list <- list()
  
#   for (method in l2_methods) {
#     tryCatch({
#       message(sprintf("Training L2 candidate: %s", method))
#       model_list[[method]] <- caret::train(
#         x = l2_train_df_filtered, 
#         y = l2_target_filtered,   
#         method = method,
#         trControl = train_control,
#         metric = metric,
#         tuneLength = 5 
#       )
#     }, error = function(e) {
#       warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
#     })
#   }

#   if (length(model_list) == 0) {
#     stop("No L2 models were successfully trained.")
#   }

#   model_comparison <- caret::resamples(model_list)
#   print(summary(model_comparison))

#   best_metric_values <- sapply(model_list, function(m) max(m$results[, metric], na.rm = TRUE))
#   best_model_name <- names(which.max(best_metric_values))
#   best_model <- model_list[[best_model_name]]
  
#   message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
#                   best_model_name, metric, max(best_model$results[, metric], na.rm = TRUE)))
  
#   # === 6. 결과 반환 ===
#   return(list(
#     best_model = best_model,
#     best_model_name = best_model_name,
#     model_comparison = model_comparison,
#     l2_train_df_filtered = data.frame(l2_train_df_filtered, target = l2_target_filtered),
#     l1_signatures = l1_signatures
#   ))
# }

# # (score_signature 함수는 이전의 'Robust Fix' 버전을 사용한다고 가정)
# # (caret, ranger, xgboost 등 필요 패키지 로드 가정)

# train_meta_learner_v4 <- function(l1_signatures,
#                                holdout_data,
#                                target_var,
#                                l2_methods = c("glm", "ranger", "xgbTree"),
#                                k_folds = 5,
#                                metric = "AUC",
#                                fgs_seed = 42) {
  
#   if (!requireNamespace("caret", quietly = TRUE)) {
#     stop("caret package required.")
#   }
#   set.seed(fgs_seed)

#   # === 1. 데이터 추출 ===
#   message("--- Preparing holdout data (Extracting once) ---")
#   if (inherits(holdout_data, "Seurat")) {
#     expr_mat <- Seurat::GetAssayData(holdout_data, layer="data")
#     meta_data <- holdout_data@meta.data
#   } else {
#     expr_mat <- holdout_data
#     if (!is.vector(target_var) && !is.factor(target_var)) {
#         stop("If holdout_data is a matrix, target_var must be the actual target vector.")
#     }
#     meta_data <- NULL 
#   }

#   # === 2. L2 Feature 생성 ===
#   message("--- Generating L2 features from holdout data ---")
#   l2_features_list <- lapply(l1_signatures, function(sig) {
#     score_signature(expr_data = expr_mat, signature = sig, normalize = TRUE) 
#   })
#   l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
#   if (is.null(names(l1_signatures))) { names(l1_signatures) <- paste0("L1_model_", seq_along(l1_signatures)) }
#   colnames(l2_train_df) <- make.names(names(l1_signatures))

#   # === 3. L2 타겟 변수 준비 및 소독 (BUG FIX 3) ===
#   if (!is.null(meta_data)) {
#     if (!all(rownames(l2_train_df) == rownames(meta_data))) {
#         meta_data <- meta_data[rownames(l2_train_df), , drop = FALSE]
#     }
#     l2_target <- meta_data[[target_var]]
#   } else { l2_target <- target_var }
#   if (!is.factor(l2_target)) { l2_target <- factor(l2_target) }
  
#   original_levels <- levels(l2_target)
#   sanitized_levels <- make.names(original_levels)
#   if (!all(original_levels == sanitized_levels)) {
#     message(sprintf("Sanitizing target levels for caret: '%s' -> '%s'", 
#                     paste(original_levels, collapse="' / '"), 
#                     paste(sanitized_levels, collapse="' / '")))
#     levels(l2_target) <- sanitized_levels
#   }

#   # === 4. NA 타겟 제거 (BUG FIX 2) ===
#   message(sprintf("Initial L2 target length: %d", length(l2_target)))
#   na_target_indices <- is.na(l2_target)
#   if (any(na_target_indices)) {
#     n_na <- sum(na_target_indices)
#     message(sprintf("Found and removing %d cells with NA target.", n_na))
#     l2_target <- l2_target[!na_target_indices]
#     l2_train_df <- l2_train_df[!na_target_indices, , drop = FALSE]
#   }
  
#   # === 5. NA/NaN/Inf Feature 제거 (BUG FIX 4) ===
#   message("Checking for NA/NaN/Inf in L2 features (scores)...")
#   complete_feature_indices <- complete.cases(l2_train_df)
#   if (!all(complete_feature_indices)) {
#     n_na_features <- sum(!complete_feature_indices)
#     message(sprintf("Found and removing %d cells with NA/NaN/Inf features.", n_na_features))
#     l2_target_filtered <- l2_target[complete_feature_indices]
#     l2_train_df_filtered <- l2_train_df[complete_feature_indices, , drop = FALSE]
#   } else {
#     l2_target_filtered <- l2_target
#     l2_train_df_filtered <- l2_train_df
#   }

#   # === (BUG FIX 5): Zero-Variance Predictor 제거 ===
#   message("Checking for zero-variance predictors (constant features)...")
#   nzv_indices <- caret::nearZeroVar(l2_train_df_filtered, saveMetrics = FALSE)
  
#   if (length(nzv_indices) > 0) {
#     n_nzv <- length(nzv_indices)
#     nzv_names <- colnames(l2_train_df_filtered)[nzv_indices]
#     message(sprintf("Found and removing %d constant features: %s",
#                     n_nzv, paste(nzv_names, collapse=", ")))
    
#     l2_train_df_filtered <- l2_train_df_filtered[, -nzv_indices, drop = FALSE]
#   }
#   # === End of fix ===

#   message(sprintf("Final L2 training set size: %d cells, %d features", 
#                   length(l2_target_filtered), ncol(l2_train_df_filtered)))
                  
#   if (length(l2_target_filtered) == 0 || nrow(l2_train_df_filtered) == 0) {
#       stop("No valid data remaining after removing NAs from target and/or features.")
#   }
#   if (ncol(l2_train_df_filtered) == 0) {
#       stop("No valid features remaining after removing zero-variance predictors.")
#   }

#   # === 6. L2 모델 교차 검증 및 선택 ===
#   message("--- Training and comparing L2 models via k-fold CV ---")
  
#   if (length(levels(l2_target_filtered)) == 2) {
#     message("Binary target detected. Using twoClassSummary (for AUC).")
#     train_control <- caret::trainControl(
#       method = "cv", number = k_folds,
#       summaryFunction = caret::twoClassSummary,
#       classProbs = TRUE,
#       savePredictions = "final",
#       allowParallel = TRUE
#     )
#     if (metric != "AUC") { message("NOTE: Target is binary, forcing metric to 'AUC'."); metric <- "AUC" }
#   } else {
#     message("Multi-class target detected. Using defaultSummary.")
#     train_control <- caret::trainControl(method = "cv", number = k_folds, summaryFunction = caret::defaultSummary, savePredictions = "final", allowParallel = TRUE)
#     if (metric == "AUC") { message("WARNING: Target is multi-class, 'AUC' metric is invalid. Defaulting to 'Accuracy'."); metric <- "Accuracy" }
#   }
  
#   model_list <- list()
#   for (method in l2_methods) {
#     tryCatch({
#       message(sprintf("Training L2 candidate: %s", method))
#       model_list[[method]] <- caret::train(
#         x = l2_train_df_filtered,
#         y = l2_target_filtered,
#         method = method,
#         trControl = train_control,
#         metric = metric,
#         tuneLength = 5
#         # 'na.action' 인수 제거됨
#       )
#     }, error = function(e) {
#       warning(sprintf("Failed to train L2 model: %s. Error: %s", method, e$message))
#     })
#   }

#   if (length(model_list) < 2) {
#     if (length(model_list) == 0) {
#         stop("No L2 models were successfully trained.")
#     }
#     warning("Only one L2 model was successfully trained. Comparison skipped.")
#     model_comparison <- NULL
#     best_model <- model_list[[1]]
#     best_model_name <- names(model_list)[1]
    
#   } else {
#     model_comparison <- caret::resamples(model_list)
#     print(summary(model_comparison))
#     best_metric_values <- sapply(model_list, function(m) max(m$results[, metric], na.rm = TRUE))
#     best_model_name <- names(which.max(best_metric_values))
#     best_model <- model_list[[best_model_name]]
#   }
  
#   # 최종 에러 방지 (best_model이 메트릭 생성에 실패한 경우)
#   best_metric_val <- tryCatch({
#     max(best_model$results[, metric], na.rm = TRUE)
#   }, error = function(e) {
#     warning(sprintf("Could not extract metric '%s' from best model '%s'.", metric, best_model_name))
#     return(NA_real_)
#   })

#   message(sprintf("--- Best L2 model selected: %s (CV %s: %.4f) ---",
#                   best_model_name, metric, best_metric_val))
  
#   # === 7. 결과 반환 ===
#   return(list(
#     best_model = best_model,
#     best_model_name = best_model_name,
#     model_comparison = model_comparison,
#     l2_train_df_filtered = data.frame(l2_train_df_filtered, target = l2_target_filtered),
#     l1_signatures = l1_signatures
#   ))
# }

# #' @export
# train_meta_learner_v5 <- function(...) {
#   .Deprecated("TML6", package = "myR")
#   TML6(...)
# }

# #' @export
# train_meta_learner_v6 <- function(...) {
#   .Deprecated("TML6", package = "myR")
#   TML6(...)
# }

# scatter_smooth_colored3_in_test <- function(object,
#                                    feature,
#                                    group.by   = "sample_no",
#                                    x_var      = "nih_change",
#                                    transpose  = FALSE,
#                                    color_by   = NULL,
#                                    split.by   = NULL,
#                                    palette    = NULL,
#                                    transparency    = TRUE,
#                                    transparency_desc = FALSE,
#                                    fitted_line = c("linear", "loess", "lasso", NULL)) {
#     fitted_line <- match.arg(fitted_line)
#     stopifnot(is.character(feature), length(feature) == 1)
    
#     # 1. Build per‑cell tibble ------------------------------------------------
    
#     if (inherits(object, "Seurat")) {
#         expr_vec <- Seurat::FetchData(object, vars = feature)[, 1]
#         meta_df  <- tibble::as_tibble(object@meta.data)
#         cell_df  <- dplyr::mutate(meta_df, !!feature := expr_vec)
#     } else {
#         cell_df <- tibble::as_tibble(object)
#         if (!feature %in% names(cell_df))
#             stop("In data.frame mode the column '", feature, "' must exist.")
#     }
    
#     for (col in c(group.by, x_var, color_by, split.by)) {
#         if (!is.null(col) && !col %in% names(cell_df))
#             stop("Column '", col, "' not found in data.")
#     }
    
#     if (!is.null(split.by)) {
#         if (is.numeric(cell_df[[split.by]])) {
#             stop("`split.by` column '", split.by, "' must be categorical (character or factor), not numeric.")
#         }
#         cell_df[[split.by]] <- as.factor(cell_df[[split.by]])
#     }
    
    
#     # 2. Aggregate by group ---------------------------------------------------
    
#     is_color_numeric <- if (!is.null(color_by)) {
#         is.numeric(cell_df[[color_by]])
#     } else {
#         FALSE
#     }
    
#     agg_df <- cell_df %>%
#         dplyr::group_by(.data[[group.by]]) %>%
#         dplyr::summarise(
#             avg_expr = mean(.data[[feature]], na.rm = TRUE),
#             x_val    = mean(.data[[x_var]], na.rm = TRUE),
#             colour   = if (!is.null(color_by)) {
#                 if (is_color_numeric) {
#                     mean(.data[[color_by]], na.rm = TRUE)
#                 } else {
#                     dplyr::first(.data[[color_by]], na_rm = TRUE)
#                 }
#             } else {
#                 NA
#             },
#             split_col = if (!is.null(split.by)) {
#                 dplyr::first(.data[[split.by]], na_rm = TRUE)
#             } else {
#                 NA
#             },
#             .groups  = "drop"
#         )
    
    
#     # 3. Aesthetics -----------------------------------------------------------
    
#     x_col <- if (transpose) "avg_expr" else "x_val"
#     y_col <- if (transpose) "x_val"   else "avg_expr"
    
#     p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
    
#     if (!is.null(split.by)) {
#         if (!is.null(color_by)) {
#             warning("`color_by` argument is ignored when `split.by` is provided.")
#         }
#         p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data[["split_col"]]), size = 3)
        
#         pal <- palette %||% RColorBrewer::brewer.pal(max(3, length(unique(agg_df$split_col))), "Set1")
#         p <- p + ggplot2::scale_colour_manual(values = pal, name = split.by)
        
#     } else if (!is.null(color_by)) {
#         if (is.numeric(cell_df[[color_by]])) {
#             p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour,
#                                                       alpha  = colour), size = 3)
#             alpha_range <- if (transparency_desc) c(1, 0.2) else c(0.2, 1)
#             if (transparency) {
#                 p <- p + ggplot2::scale_alpha(range = alpha_range, guide = "none")
#             } else {
#                 p <- p + ggplot2::guides(alpha = "none")
#             }
#             pal <- if (is.null(palette)) viridisLite::viridis(256) else palette
#             p <- p + ggplot2::scale_colour_gradientn(colours = pal, name = color_by)
#         } else {
#             p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour), size = 3)
#             pal <- palette %||% RColorBrewer::brewer.pal(max(3, length(unique(agg_df$colour))), "Set1")
#             p <- p + ggplot2::scale_colour_manual(values = pal, name = color_by)
#         }
#     } else {
#         p <- p + ggplot2::geom_point(size = 3)
#     }
    
    
#     # 4. Smoothing line -------------------------------------------------------
    
#     if (!is.null(fitted_line)) {
        
#         if (!is.null(split.by)) {
            
#             # *** 수정된 부분 시작 ***
#             # `fitted_line` 인수를 `geom_smooth`가 이해하는 `method` 문자열로 변환
            
#             method_val <- NULL # 초기화
            
#             if (fitted_line == "linear") {
#                 method_val <- "lm"
#             } else if (fitted_line == "loess") {
#                 method_val <- "loess"
#             } else if (fitted_line == "lasso") {
#                 warning("ggplot2::geom_smooth does not support 'lasso'. Falling back to 'lm' method.")
#                 method_val <- "lm"
#             }
#             # *** 수정된 부분 끝 ***
            
#             warning("Custom statistical annotations (p-value, equation) are disabled when `split.by` is used.")
            
#             p <- p + ggplot2::geom_smooth(ggplot2::aes(colour = .data[["split_col"]]),
#                                           method = method_val, # "linear" 대신 "lm"이 전달됨
#                                           se = TRUE,
#                                           show.legend = FALSE) 
            
#         } else {
#             # (원본 로직: split.by가 NULL일 때만 실행)
#             if (fitted_line == "linear") {
#                 p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black")
#                 fit <- stats::lm(agg_df[[y_col]] ~ agg_df[[x_col]])
#                 coef <- round(stats::coef(fit), 3)
#                 pval <- signif(summary(fit)$coefficients[2, 4], 3)
#                 annot <- paste0("y = ", coef[1], " + ", coef[2], " * x\np = ", pval)
#                 p <- p + ggplot2::annotate("text", x = min(agg_df[[x_col]], na.rm = TRUE),
#                                            y = max(agg_df[[y_col]], na.rm = TRUE),
#                                            label = annot, hjust = 0, vjust = 1, size = 4)
#             } else if (fitted_line == "loess") {
#                 p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "black")
#             } else if (fitted_line == "lasso") {
#                 if (!requireNamespace("glmnet", quietly = TRUE)) {
#                     warning("glmnet not installed; falling back to linear fit.")
#                     p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE)
#                 } else {
#                     xmat <- as.matrix(agg_df[[x_col]])
#                     fit  <- glmnet::cv.glmnet(xmat, agg_df[[y_col]], alpha = 1)
#                     preds <- as.numeric(glmnet::predict.glmnet(fit$glmnet.fit, newx = xmat,
#                                                                s = fit$lambda.min))
#                     pred_df <- agg_df %>% dplyr::mutate(pred = preds)
#                     p <- p + ggplot2::geom_line(data = pred_df[order(pred_df[[x_col]]), ],
#                                                 ggplot2::aes(x = .data[[x_col]], y = pred),
#                                                 colour = "red", linewidth = 1)
#                 }
#             }
#         }
#     }
    
    
#     # 5. Labels & theme -------------------------------------------------------
    
#     p <- p + ggplot2::theme_bw() +
#         ggplot2::labs(x = if (transpose) paste("Average", feature, "expression") else x_var,
#                       y = if (transpose) x_var else paste("Average", feature, "expression"))
    
#     p
#     return(p)
# }

# plot_volcano_in_test <- function(lmm_summary,
#                          x_col="estimate",
#                          y_col="p_value",
#                          filter_col="term",
#                         filter_pattern = NULL,
#                         title = "Volcano Plot",
#                         effect_threshold = 0.5,
#                         p_threshold = 0.05,
#                         resize_x=TRUE,
#                         resize_y=TRUE) {
#   #ver2: generalized
#   plot_data <- lmm_summary
  
#   if (!is.null(filter_pattern)) {
#     plot_data <- plot_data %>%
#       filter(grepl(filter_pattern, term))
#   }
#   if (!is.null(x_col)){plot_data=plot_data%>%mutate(effect_size= !!sym(x_col))}
#   if (!is.null(y_col)){plot_data=plot_data%>%mutate(p_value= !!sym(y_col))}

#   # plot label definition
#   x_label="Effect Size"
#   y_label="Significance"
  
#   # plot resizing
#   if(resize_x){xlim_val=quantile(abs(plot_data$effect_size),0.95,na.rm=TRUE)}
#   if(resize_y){ylim_val=quantile(abs(plot_data$effect_size),0.95,na.rm=TRUE)}
  
#   # gene significance label arrange
#   plot_data <- plot_data %>%
#     mutate(
#       log_p = -log10(p_value),
#       category = case_when(
#         abs(effect_size) > effect_threshold & p_value < p_threshold ~ "Significant",
#         abs(effect_size) > effect_threshold ~ "Large effect",
#         p_value  < p_threshold ~ "Small effect",
#         TRUE ~ "Not significant"
#       )
#     )
#   ggplot(plot_data, aes(x = effect_size, y = log_p, color = category)) +
#     geom_point(alpha = 0.6) +
#     geom_hline(yintercept = -log10(p_threshold), 
#                linetype = "dashed", color = "gray") +
#     geom_vline(xintercept = c(-effect_threshold, effect_threshold), 
#                linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("Significant" = "red",
#                                  "Large effect" = "orange", 
#                                  "Small effect" = "blue",
#                                  "Not significant" = "gray")) +
#     labs(title = title,
#          x = x_label,
#          y = y_label) +
#     theme_bw()+
#     coord_cartesian(xlim = c(-xlim_val, xlim_val))
# }


# # validation tools (not using) -----

# #' Validate Seurat Object
# #'
# #' Checks if the input is a valid Seurat object with expected components.
# #'
# #' @param obj Object to validate
# #' @param assay Optional assay name to check for
# #' @param reduction Optional reduction name to check for
# #' @param min_cells Minimum number of cells required
# #' @param min_features Minimum number of features required
# #'
# #' @return TRUE if valid, stops with error message otherwise
# #' @keywords internal
# #' @export
# validate_seurat <- function(obj, 
#                             assay = NULL, 
#                             reduction = NULL,
#                             min_cells = 0,
#                             min_features = 0) {
  
#   # Check if it's a Seurat object
#   if (!inherits(obj, "Seurat")) {
#     stop("Input must be a Seurat object. Current class: ", 
#          paste(class(obj), collapse = ", "))
#   }
  
#   # Check cell count
#   n_cells <- ncol(obj)
#   if (n_cells < min_cells) {
#     stop("Seurat object has ", n_cells, " cells, but ", min_cells, " required")
#   }
  
#   # Check feature count
#   n_features <- nrow(obj)
#   if (n_features < min_features) {
#     stop("Seurat object has ", n_features, " features, but ", min_features, " required")
#   }
  
#   # Check assay if specified
#   if (!is.null(assay)) {
#     available_assays <- Seurat::Assays(obj)
#     if (!assay %in% available_assays) {
#       stop("Assay '", assay, "' not found. Available assays: ", 
#            paste(available_assays, collapse = ", "))
#     }
#   }
  
#   # Check reduction if specified
#   if (!is.null(reduction)) {
#     available_reductions <- names(obj@reductions)
#     if (!reduction %in% available_reductions) {
#       # Try case-insensitive match
#       matched_reduction <- available_reductions[tolower(available_reductions) == tolower(reduction)]
#       if (length(matched_reduction) == 1) {
#         message("Note: Using case-insensitive match: '", matched_reduction, "' for '", reduction, "'")
#         return(TRUE)
#       }
#       stop("Reduction '", reduction, "' not found. Available reductions: ", 
#            paste(available_reductions, collapse = ", "))
#     }
#   }
  
#   return(TRUE)
# }

# #' Validate Metadata Column
# #'
# #' Checks if a metadata column exists and optionally validates its type.
# #'
# #' @param obj Seurat object or data frame
# #' @param column_name Column name to validate
# #' @param required_type Optional required type ("numeric", "factor", "character")
# #' @param allow_na Whether NA values are allowed
# #'
# #' @return TRUE if valid, stops with error message otherwise
# #' @keywords internal
# #' @export
# validate_metadata_column <- function(obj, 
#                                      column_name, 
#                                      required_type = NULL,
#                                      allow_na = TRUE) {
  
#   # Get metadata
#   if (inherits(obj, "Seurat")) {
#     metadata <- obj@meta.data
#   } else if (is.data.frame(obj)) {
#     metadata <- obj
#   } else {
#     stop("Input must be a Seurat object or data frame")
#   }
  
#   # Check if column exists
#   if (!column_name %in% colnames(metadata)) {
#     stop("Column '", column_name, "' not found in metadata. Available columns: ",
#          paste(head(colnames(metadata), 10), collapse = ", "),
#          if (ncol(metadata) > 10) "..." else "")
#   }
  
#   column_data <- metadata[[column_name]]
  
#   # Check for NA values if not allowed
#   if (!allow_na && any(is.na(column_data))) {
#     n_na <- sum(is.na(column_data))
#     stop("Column '", column_name, "' contains ", n_na, " NA values, but NA not allowed")
#   }
  
#   # Check type if specified
#   if (!is.null(required_type)) {
#     is_correct_type <- switch(
#       required_type,
#       "numeric" = is.numeric(column_data),
#       "factor" = is.factor(column_data),
#       "character" = is.character(column_data),
#       stop("Unknown required_type: ", required_type)
#     )
    
#     if (!is_correct_type) {
#       stop("Column '", column_name, "' must be of type '", required_type, 
#            "', but is: ", paste(class(column_data), collapse = ", "))
#     }
#   }
  
#   return(TRUE)
# }

# #' Validate Gene List
# #'
# #' Checks if genes exist in the object and optionally filters to valid genes.
# #'
# #' @param obj Seurat object or character vector of available genes
# #' @param genes Character vector of genes to validate
# #' @param min_present Minimum number of genes that must be present (default: all)
# #' @param assay Assay to check genes in (for Seurat objects)
# #' @param warn_missing Whether to warn about missing genes
# #'
# #' @return Character vector of valid genes present in the object
# #' @keywords internal
# #' @export
# validate_genes <- function(obj, 
#                           genes, 
#                           min_present = NULL,
#                           assay = NULL,
#                           warn_missing = TRUE) {
  
#   if (!is.character(genes) || length(genes) == 0) {
#     stop("genes must be a non-empty character vector")
#   }
  
#   # Get available genes
#   if (inherits(obj, "Seurat")) {
#     if (is.null(assay)) {
#       assay <- Seurat::DefaultAssay(obj)
#     }
#     available_genes <- rownames(obj[[assay]])
#   } else if (is.character(obj)) {
#     available_genes <- obj
#   } else {
#     stop("obj must be a Seurat object or character vector of gene names")
#   }
  
#   # Find present and missing genes
#   genes_present <- intersect(genes, available_genes)
#   genes_missing <- setdiff(genes, available_genes)
  
#   # Set default minimum
#   if (is.null(min_present)) {
#     min_present <- length(genes)
#   }
  
#   # Check if enough genes are present
#   if (length(genes_present) < min_present) {
#     stop("Only ", length(genes_present), " of ", length(genes), 
#          " genes found, but ", min_present, " required. ",
#          "First missing genes: ", 
#          paste(head(genes_missing, 5), collapse = ", "),
#          if (length(genes_missing) > 5) "..." else "")
#   }
  
#   # Warn about missing genes if requested
#   if (warn_missing && length(genes_missing) > 0) {
#     warning("Genes not found (", length(genes_missing), "/", length(genes), "): ",
#             paste(head(genes_missing, 10), collapse = ", "),
#             if (length(genes_missing) > 10) "..." else "")
#   }
  
#   return(genes_present)
# }

# #' Validate Numeric Range
# #'
# #' Checks if a numeric parameter is within acceptable range.
# #'
# #' @param value Numeric value to validate
# #' @param param_name Name of parameter (for error messages)
# #' @param min Minimum allowed value (inclusive)
# #' @param max Maximum allowed value (inclusive)
# #' @param allow_na Whether NA is allowed
# #'
# #' @return TRUE if valid, stops with error message otherwise
# #' @keywords internal
# #' @export
# validate_numeric_range <- function(value, 
#                                    param_name, 
#                                    min = -Inf, 
#                                    max = Inf,
#                                    allow_na = FALSE) {
  
#   # Check for NA
#   if (is.na(value)) {
#     if (allow_na) {
#       return(TRUE)
#     } else {
#       stop(param_name, " cannot be NA")
#     }
#   }
  
#   # Check if numeric
#   if (!is.numeric(value) || length(value) != 1) {
#     stop(param_name, " must be a single numeric value")
#   }
  
#   # Check range
#   if (value < min || value > max) {
#     stop(param_name, " must be between ", min, " and ", max, 
#          ", but is: ", value)
#   }
  
#   return(TRUE)
# }

# #' Validate Choice
# #'
# #' Validates that a value is one of allowed choices (like match.arg but more informative).
# #'
# #' @param value Value to validate
# #' @param param_name Name of parameter (for error messages)
# #' @param choices Vector of allowed values
# #' @param multiple Whether multiple choices are allowed
# #'
# #' @return The validated value (or values if multiple=TRUE)
# #' @keywords internal
# #' @export
# validate_choice <- function(value, param_name, choices, multiple = FALSE) {
  
#   if (is.null(value)) {
#     stop(param_name, " cannot be NULL")
#   }
  
#   if (multiple) {
#     if (!all(value %in% choices)) {
#       invalid <- setdiff(value, choices)
#       stop(param_name, " contains invalid choices: ", 
#            paste(invalid, collapse = ", "),
#            ". Allowed choices: ", 
#            paste(choices, collapse = ", "))
#     }
#   } else {
#     if (length(value) != 1) {
#       stop(param_name, " must be a single value, not ", length(value))
#     }
#     if (!value %in% choices) {
#       stop(param_name, " must be one of: ", 
#            paste(choices, collapse = ", "),
#            ". Got: ", value)
#     }
#   }
  
#   return(value)
# }

# #' Validate File Path
# #'
# #' Checks if a file path exists and optionally validates extension.
# #'
# #' @param path File path to validate
# #' @param must_exist Whether file must already exist
# #' @param extensions Optional vector of allowed extensions (e.g., c("csv", "txt"))
# #' @param type Type of path ("file" or "directory")
# #'
# #' @return Normalized path if valid, stops with error otherwise
# #' @keywords internal
# #' @export
# validate_path <- function(path, 
#                          must_exist = TRUE, 
#                          extensions = NULL,
#                          type = c("file", "directory")) {
  
#   type <- match.arg(type)
  
#   if (!is.character(path) || length(path) != 1) {
#     stop("path must be a single character string")
#   }
  
#   # Check existence if required
#   if (must_exist) {
#     if (type == "file") {
#       if (!file.exists(path)) {
#         stop("File not found: ", path)
#       }
#       if (dir.exists(path)) {
#         stop("Expected a file but got a directory: ", path)
#       }
#     } else if (type == "directory") {
#       if (!dir.exists(path)) {
#         stop("Directory not found: ", path)
#       }
#     }
#   }
  
#   # Check extension if specified
#   if (!is.null(extensions) && type == "file") {
#     file_ext <- tools::file_ext(path)
#     if (!file_ext %in% extensions) {
#       stop("File extension must be one of: ", 
#            paste(extensions, collapse = ", "),
#            ". Got: ", file_ext)
#     }
#   }
  
#   # Return normalized path
#   return(normalizePath(path, mustWork = must_exist))
# }

# #' Create Informative Error Message
# #'
# #' Helper to create consistent, informative error messages.
# #'
# #' @param context Context where error occurred (e.g., function name)
# #' @param message Main error message
# #' @param suggestion Optional suggestion for fixing the error
# #'
# #' @return Formatted error message
# #' @keywords internal
# #' @export
# create_error_message <- function(context, message, suggestion = NULL) {
#   msg <- paste0("[", context, "] ", message)
#   if (!is.null(suggestion)) {
#     msg <- paste0(msg, "\nSuggestion: ", suggestion)
#   }
#   return(msg)
# }

# #' Check Package Dependencies
# #'
# #' Checks if required packages are installed and optionally loads them.
# #'
# #' @param packages Character vector of package names
# #' @param load Whether to load the packages (default: FALSE)
# #'
# #' @return TRUE if all packages available, stops with error otherwise
# #' @keywords internal
# #' @export
# check_packages <- function(packages, load = FALSE) {
  
#   missing <- character(0)
  
#   for (pkg in packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       missing <- c(missing, pkg)
#     } else if (load) {
#       library(pkg, character.only = TRUE)
#     }
#   }
  
#   if (length(missing) > 0) {
#     stop("Required packages not installed: ", 
#          paste(missing, collapse = ", "),
#          "\nInstall with: install.packages(c('", 
#          paste(missing, collapse = "', '"), "'))")
#   }
  
#   return(TRUE)
# }


# ============================================================================
# DEPRECATED: 이 파일의 모든 함수들은 test_analysis.R로 통합되었습니다.
# 최신 버전은 test_analysis.R의 다음 함수들을 사용하세요:
# - runNEBULA_v1 (또는 runNEBULA2_v1)
# - runMAST_v1
# - runMUSCAT_v5 (또는 runMUSCAT2_v1)
# ============================================================================

# #' @title NEBULA 파이프라인 함수 (DEPRECATED)
# #' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runNEBULA_v1을 사용하세요.
# #' @export 
# runNEBULA <- function(...) {
#   .Deprecated("runNEBULA_v1", package = "myR", msg = "runNEBULA is deprecated. Use runNEBULA_v1 from test_analysis.R instead.")
#   if (exists("runNEBULA_v1", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runNEBULA_v1", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runNEBULA_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

# #' @title MAST 파이프라인 함수 (DEPRECATED)
# #' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMAST_v1을 사용하세요.
# #' @export
# runMAST_v1 <- function(...) {
#   .Deprecated("runMAST_v1", package = "myR", msg = "runMAST_v1 in test_cursor.R is deprecated. Use runMAST_v1 from test_analysis.R instead.")
#   if (exists("runMAST_v1", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runMAST_v1", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runMAST_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

# #' @title Muscat 파이프라인 함수 v1 (DEPRECATED)
# #' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMUSCAT_v5 또는 runMUSCAT2_v1을 사용하세요.
# #' @export
# runMUSCAT_v1 <- function(...) {
#   .Deprecated("runMUSCAT_v5", package = "myR", msg = "runMUSCAT_v1 is deprecated. Use runMUSCAT_v5 or runMUSCAT2_v1 from test_analysis.R instead.")
#   if (exists("runMUSCAT_v5", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runMUSCAT_v5", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runMUSCAT_v5 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

# #' @title MAST 파이프라인 함수 (DEPRECATED)
# #' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMAST_v1을 사용하세요.
# #' @export
# runMAST <- function(...) {
#   .Deprecated("runMAST_v1", package = "myR", msg = "runMAST is deprecated. Use runMAST_v1 from test_analysis.R instead.")
#   if (exists("runMAST_v1", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runMAST_v1", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runMAST_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

# #' @title Muscat 파이프라인 함수 (DEPRECATED)
# #' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMUSCAT_v5 또는 runMUSCAT2_v1을 사용하세요.
# #' @export
# runMUSCAT <- function(...) {
#   .Deprecated("runMUSCAT_v5", package = "myR", msg = "runMUSCAT is deprecated. Use runMUSCAT_v5 or runMUSCAT2_v1 from test_analysis.R instead.")
#   if (exists("runMUSCAT_v5", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runMUSCAT_v5", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runMUSCAT_v5 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

#=========================




# from test_working.R
# #' @title NEBULA 파이프라인 함수
# #' @description Seurat 객체를 받아 NA 처리, 데이터 정렬, NEBULA GLMM을 실행.
# #'              Modeling에서, y~NB(mu,phi)이며, log(mu)=log(o)+BX+random effect 에서 "o"는 offset, random effect는 patient_col, X는 fixed effects와 covar_effects.
# #'
# #' @param sobj (Seurat) Seurat 객체
# #' @param layer (character) Raw count가 저장된 layer (예: "counts")
# #' @param fixed_effects (character vector) 고정 효과로 사용할 meta.data 컬럼명
# #'                    (예: c("target_col", "celltype_col"))
# #' @param covar_effects (character vector) 임의 효과처럼 보정하고 싶은 공변량.
# #'                    NEBULA의 한계로 인해 '고정 효과'로 처리됩니다.
# #'                    (예: c("batch_col", "percent.mt"))
# #' @param patient_col (character) 환자 ID 컬럼. NEBULA의 'id' (유일한 임의 효과)로 사용됨.
# #' @param offset (character) Offset으로 사용할 'numeric' 컬럼 (예: "nCount_RNA")
# #' @param min_count (numeric) 유전자 필터링 기준 (최소 발현 세포 수)
# #'
# #' @return (nebula) NEBULA 실행 결과 객체
# #'
# #' @export 
# runNEBULA <- function(...) {
#   .Deprecated("runNEBULA_v1", package = "myR", msg = "runNEBULA in test_working.R is deprecated. Use runNEBULA_v1 from test_analysis.R instead.")
#   if (exists("runNEBULA_v1", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runNEBULA_v1", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runNEBULA_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

# #' @title Muscat (Pseudo-bulking) 파이프라인 함수 (V5 - prepData 사용)
# #'
# #' @description prepData를 사용하여 SCE 객체를 포맷하고, pbDS로 GLM(edgeR/DESeq2) 실행
# #' @note 이 함수는 'pbDS'의 한계로 인해 혼합 효과(Random Effects/block)를 지원하지 않습니다.
# #'
# #' @param sobj (Seurat) Seurat 객체
# #' @param cluster_id (character) 세포 유형 컬럼명
# #' @param sample_id (character) 환자/샘플 ID 컬럼명
# #' @param group_id (character) 비교할 그룹 컬럼명
# #' @param formula_str (character) 포뮬러 (예: "~ 0 + group + set")
# #'                 'group' 키워드는 'group_id'로 자동 변환됩니다.
# #' @param contrast (character) Contrast (예: "groupA-groupB")
# #' @param method (character) "edgeR", "DESeq2", "limma-trend", "limma-voom" (block 없음)
# #'
# #' @export
# runMUSCAT <- function(...) {
#   .Deprecated("runMUSCAT_v5", package = "myR", msg = "runMUSCAT in test_working.R is deprecated. Use runMUSCAT_v5 or runMUSCAT2_v1 from test_analysis.R instead.")
#   if (exists("runMUSCAT_v5", envir = asNamespace("myR"), inherits = FALSE)) {
#     fun <- get("runMUSCAT_v5", envir = asNamespace("myR"))
#     return(fun(...))
#   } else {
#     stop("runMUSCAT_v5 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
#   }
# }

# # Original runMUSCAT function body removed - use runMUSCAT_v5 from test_analysis.R instead
# # The following code block is kept for reference but will not execute:
# if (FALSE) {
#   runMUSCAT_original <- function(
#   sobj,
#   cluster_id = "seurat_clusters",
#   sample_id  = "hos_no",
#   group_id   = "type",
#   batch_id   = NULL,                 # ex) "exp_batch"
#   contrast   = NULL,                 # ex) "IS - SAH"
#   method     = "edgeR",
#   pb_min_cells = 3,
#   filter_genes = c("none","genes","both","edgeR"),
#   keep_clusters = NULL,
#   cluster_label_map = NULL
# ){
#   if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")
#   filter_genes <- match.arg(filter_genes)

#   # deps
#   req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr")
#   miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
#   if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

#   # 1) Seurat -> SCE, prepSCE
#   sce <- Seurat::as.SingleCellExperiment(sobj)
#   sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

#   # factor 보장
#   sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
#   sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
#   sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
#   if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
#     sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
#   }

#   # 2) Pseudobulk
#   pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

#   # (선택) 특정 클러스터만
#   if (!is.null(keep_clusters)) {
#     keep_clusters <- as.character(keep_clusters)
#     pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
#     if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
#   }

#   # 2-1) pb 메타 보강 (sample_id / group_id / batch)
#   pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))

#   # sample_id 없으면 assay의 colnames로 복구
#   if (!"sample_id" %in% names(pb_meta)) {
#     first_assay <- names(SummarizedExperiment::assays(pb))[1]
#     sid_guess <- colnames(SummarizedExperiment::assays(pb)[[first_assay]])
#     if (is.null(sid_guess)) stop("pb에 sample_id가 없습니다.")
#     pb_meta$sample_id <- sid_guess
#     rownames(pb_meta) <- pb_meta$sample_id
#     SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)
#   }

#   # sce에서 (sample_id -> group_id / batch) map
#   sce_meta <- as.data.frame(SummarizedExperiment::colData(sce))
#   map_cols <- c("sample_id","group_id")
#   if (!is.null(batch_id) && batch_id %in% names(sce_meta)) map_cols <- c(map_cols, batch_id)
#   sce_map <- unique(sce_meta[, map_cols, drop=FALSE])

#   # pb에 group_id / batch 보강
#   pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
#   need_fix <- (!"group_id" %in% names(pb_meta)) ||
#               (length(unique(pb_meta$group_id)) < 2) ||
#               (all(unique(pb_meta$group_id) %in% c("type","group","group_id", NA, "")))
#   if (need_fix || (!is.null(batch_id) && !batch_id %in% names(pb_meta))) {
#     pb_meta2 <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
#     if ("group_id.x" %in% names(pb_meta2) && "group_id.y" %in% names(pb_meta2)) {
#       pb_meta2$group_id <- ifelse(is.na(pb_meta2$group_id.y), pb_meta2$group_id.x, pb_meta2$group_id.y)
#       pb_meta2$group_id.x <- NULL; pb_meta2$group_id.y <- NULL
#     }
#     rownames(pb_meta2) <- rownames(pb_meta)
#     SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta2)
#   }

#   # factor화
#   pb$sample_id <- droplevels(factor(SummarizedExperiment::colData(pb)$sample_id))
#   pb$group_id  <- droplevels(factor(SummarizedExperiment::colData(pb)$group_id))
#   if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb))) {
#     pb[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[batch_id]]))
#   }

#   # 3) contrast 그룹만 자동 subset
#   extract_groups <- function(contrast_str, levels_available){
#     z <- gsub("\\s+", "", contrast_str)
#     toks <- unique(gsub("^group(_id)?", "", unlist(strsplit(z, "[^A-Za-z0-9_]+"))))
#     toks <- toks[nchar(toks) > 0]
#     keep <- intersect(toks, levels_available)
#     if (length(keep) < 1) {
#       g2 <- levels_available[vapply(levels_available, function(g) grepl(g, z), logical(1))]
#       keep <- unique(g2)
#     }
#     keep
#   }
#   grp_lvls <- levels(SummarizedExperiment::colData(pb)$group_id)
#   tg <- extract_groups(contrast, grp_lvls)
#   if (length(tg) < 2) stop(sprintf("contrast에서 추출한 그룹이 부족합니다. contrast='%s', 사용가능레벨=%s",
#                                    contrast, paste(grp_lvls, collapse=", ")))

#   keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
#   pb_sub <- pb[, keep_idx]
#   pb_sub$group_id <- droplevels(factor(SummarizedExperiment::colData(pb_sub)$group_id))

#   # **sce도 동일 기준으로 subset (resDS용 필수)**
#   sce_sub <- sce[, sce$sample_id %in% SummarizedExperiment::colData(pb_sub)$sample_id &
#                     sce$group_id  %in% tg]
#   sce_sub$cluster_id <- droplevels(factor(sce_sub$cluster_id))
#   sce_sub$sample_id  <- droplevels(factor(sce_sub$sample_id))
#   sce_sub$group_id   <- droplevels(factor(sce_sub$group_id))
#   if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce_sub))) {
#     sce_sub[[batch_id]] <- droplevels(factor(sce_sub[[batch_id]]))
#   }

#   # 4) design/contrast (batch는 'batch'로 복사해서 사용)
#   pb_sub$group <- pb_sub$group_id
#   if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb_sub))) {
#     pb_sub$batch <- droplevels(factor(SummarizedExperiment::colData(pb_sub)[[batch_id]]))
#     design <- stats::model.matrix(~ 0 + group + batch,
#                                   data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
#   } else {
#     design <- stats::model.matrix(~ 0 + group,
#                                   data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
#   }

#   fix_contrast <- function(contrast_str, design_cols){
#     z <- gsub("\\s+", "", contrast_str)
#     toks <- unlist(strsplit(z, "([+\\-])", perl=TRUE))
#     ops  <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl=TRUE)))
#     rebuild <- function(tok){
#       tok <- gsub("^group(_id)?", "group", tok)
#       if (!grepl("^group", tok)) tok <- paste0("group", tok)
#       tok
#     }
#     toks2 <- vapply(toks, rebuild, character(1))
#     out <- toks2[1]; if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i+1])
#     out
#   }
#   contrast_fixed <- fix_contrast(contrast, colnames(design))
#   contrast_matrix <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)

#   # 5) pbDS
#   res <- muscat::pbDS(
#     pb_sub,
#     design    = design,
#     method    = method,
#     contrast  = contrast_matrix,
#     min_cells = pb_min_cells,
#     filter    = filter_genes,
#     verbose   = TRUE
#   )

#   # 6) 결과 평탄화: **sce_sub가 먼저, res가 다음**
#   combined <- muscat::resDS(sce_sub, res)

#   # cluster_id 정리 + 라벨
#   if ("cluster" %in% names(combined) && !"cluster_id" %in% names(combined)) {
#     combined$cluster_id <- combined$cluster
#   }
#   if (!"cluster_id" %in% names(combined)) stop("resDS 결과에 'cluster_id'가 없습니다.")
#   if (!is.null(cluster_label_map)) {
#     combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
#     combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
#   } else {
#     combined$cluster_label <- as.character(combined$cluster_id)
#   }

#   # 7) 요약 헬퍼
#   .pick_cols_for <- function(tab, contrast_fixed) {
#     cstr <- gsub("\\s+", "", contrast_fixed)
#     patt <- paste0("(^|\\.)", gsub("([+\\-])", "\\\\\\1", cstr), "$")

#     # 1) 접미사 있는 형태 먼저 시도
#     logfc <- grep("^logFC(\\.|_)?", names(tab), value=TRUE); logfc <- logfc[grep(patt, logfc)]
#     padj  <- grep("^(p_adj(\\.loc|\\.glb)?|FDR)(\\.|_)?", names(tab), value=TRUE); padj <- padj[grep(patt, padj)]
#     pcol  <- grep("^(p_val|PValue)(\\.|_)?", names(tab), value=TRUE); pcol <- pcol[grep(patt, pcol)]

#     # 2) 접미사 매칭 실패 시 기본 컬럼으로 폴백
#     if (!length(logfc) && "logFC" %in% names(tab)) logfc <- "logFC"
#     if (!length(pcol)  && "p_val" %in% names(tab)) pcol  <- "p_val"
#     if (!length(padj)) {
#       if ("p_adj.loc" %in% names(tab))      padj <- "p_adj.loc"
#       else if ("p_adj.glb" %in% names(tab)) padj <- "p_adj.glb"
#       else if ("FDR" %in% names(tab))       padj <- "FDR"
#     }

#     if (!length(logfc) || !length(padj)) {
#       stop("결과 테이블에서 logFC/adj.p 컬럼을 찾지 못했습니다. resDS 출력 컬럼명을 확인하세요.")
#     }
#     list(logfc = logfc[1], padj = padj[1], p = if (length(pcol)) pcol[1] else NULL)
#   }

#   # --- [REPLACE] 클러스터별 Top-N
#   top_by_cluster <- function(n=25){
#     tab_use <- combined
#     # contrast 컬럼이 있으면 현재 contrast만 사용
#     if ("contrast" %in% names(tab_use)) {
#       target <- gsub("\\s+", "", contrast_fixed)
#       tab_use <- tab_use[gsub("\\s+", "", tab_use$contrast) == target, , drop=FALSE]
#     }
#     cols <- .pick_cols_for(tab_use, contrast_fixed)

#     tab_use |>
#       dplyr::mutate(
#         logFC_view = .data[[cols$logfc]],
#         padj_view  = .data[[cols$padj]],
#         p_view     = if (!is.null(cols$p)) .data[[cols$p]] else NA_real_
#       ) |>
#       dplyr::arrange(padj_view) |>
#       dplyr::group_by(cluster_label) |>
#       dplyr::slice_head(n = n) |>
#       dplyr::ungroup()
#   }

#   # --- [REPLACE] 전체 Top-N
#   top_overall <- function(n=100){
#     tab_use <- combined
#     if ("contrast" %in% names(tab_use)) {
#       target <- gsub("\\s+", "", contrast_fixed)
#       tab_use <- tab_use[gsub("\\s+", "", tab_use$contrast) == target, , drop=FALSE]
#     }
#     cols <- .pick_cols_for(tab_use, contrast_fixed)

#     tab_use |>
#       dplyr::mutate(
#         logFC_view = .data[[cols$logfc]],
#         padj_view  = .data[[cols$padj]],
#         p_view     = if (!is.null(cols$p)) .data[[cols$p]] else NA_real_
#       ) |>
#       dplyr::arrange(padj_view) |>
#       dplyr::slice_head(n = n)
#   }

#   list(
#     sce        = sce,
#     sce_sub    = sce_sub,
#     pb         = pb,
#     pb_sub     = pb_sub,
#     design     = design,
#     contrast_fixed = contrast_fixed,
#     res_raw    = res,
#     combined   = combined,
#     top_by_cluster = top_by_cluster,
#     top_overall    = top_overall
#   )
# }
# } # End of if (FALSE) block


# # from signature.R
# #' Find Gene Signatures that Best Separate Target Variable
# #'
# #' @param data Seurat object, count matrix, or data.frame with genes as rows/columns
# #' @param meta.data Optional metadata data.frame (required if data is not Seurat)
# #' @param target_var Column name in metadata representing the target variable
# #' @param target_group For numeric targets: quantile cutoff (0-1) or list(low=0.25, high=0.75)
# #'                     For factor: specific levels to compare (default: all vs all)
# #' @param method One of: "tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings"
# #' @param n_features Number of top features to return (default: 50)
# #' @param preprocess Logical, whether to normalize/scale data (default: TRUE)
# #' @param min_cells Minimum cells expressing gene (default: 10)
# #' @param min_pct Minimum percentage of cells expressing gene (default: 0.01)
# #' @param return_model Logical, return full model object (default: FALSE)
# #' @param seed Random seed for reproducibility (default: 42)
# #' @param ... Additional method-specific parameters
# #'
# #' @return List containing:
# #'   - genes: Character vector of selected genes
# #'   - weights: Named numeric vector of gene weights/importance
# #'   - scores: Per-cell signature scores (if applicable)
# #'   - performance: Separation metrics (AUC, accuracy, etc.)
# #'   - method: Method used
# #'   - model: Full model object (if return_model=TRUE)
# #'
# #' @examples
# #' # Random Forest for binary classification
# #' result <- find_gene_signature(seurat_obj, target_var="cell_type", 
# #'                                target_group=c("TypeA", "TypeB"), 
# #'                                method="tree_based")
# #' 
# #' # LASSO for continuous variable (top vs bottom quartile)
# #' result <- find_gene_signature(seurat_obj, target_var="pseudotime", 
# #'                                target_group=0.25, method="lasso")
# #' 
# #' # NMF for multi-class
# #' result <- find_gene_signature(count_matrix, meta.data=metadata, 
# #'                                target_var="condition", method="nmf")
# #' @export
# find_gene_signature <- function(data, 
#                                 meta.data = NULL,
#                                 target_var,
#                                 target_group = NULL,
#                                 method = c("tree_based", "lasso", "limma", 
#                                            "nmf", "wilcoxon", "gam", "pca_loadings"),
#                                 n_features = 50,
#                                 preprocess = TRUE,
#                                 min_cells = 10,
#                                 min_pct = 0.01,
#                                 return_model = FALSE,
#                                 fgs_seed = 42,
#                                 ...) {
  
#   set.seed(fgs_seed)
#   method <- match.arg(method)
  
#   # ===
#   # 1. Input validation and data extraction
#   # ===
  
#   # Check if Seurat object
#   is_seurat <- inherits(data, "Seurat")
  
#   if (is_seurat) {
#     if (!requireNamespace("Seurat", quietly = TRUE)) {
#       stop("Seurat package required but not installed")
#     }
#     if (is.null(meta.data)) {
#       meta.data <- data@meta.data
#     }
#     # Extract normalized data or raw counts
#     if ("data" %in% names(data@assays[[1]])) {
#       expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
#     } else {
#       expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
#     }
#   } else {
#     # Assume matrix or data.frame
#     if (is.null(meta.data)) {
#       stop("meta.data must be provided when data is not a Seurat object")
#     }
#     expr_mat <- as.matrix(data)
#     # Transpose if cells are rows
#     if (nrow(expr_mat) == nrow(meta.data)) {
#       expr_mat <- t(expr_mat)
#     }
#   }
  
#   # Validate target_var exists
#   if (!target_var %in% colnames(meta.data)) {
#     stop(sprintf("target_var '%s' not found in metadata columns: %s", 
#                  target_var, paste(colnames(meta.data), collapse=", ")))
#   }
  
#   # Align cells
#   common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
#   if (length(common_cells) == 0) {
#     stop("No common cells found between expression matrix and metadata")
#   }
#   expr_mat <- expr_mat[, common_cells]
#   meta.data <- meta.data[common_cells, ]
  
#   # ===
#   # 2. Process target variable
#   # ===
  
#   target_values <- meta.data[[target_var]]
#   target_type <- class(target_values)[1]
  
#   if (is.numeric(target_values)) {
#     # Continuous variable - create binary groups
#     if (is.null(target_group)) {
#       target_group <- 0.25  # default: bottom 25% vs top 25%
#     }
    
#     if (is.list(target_group)) {
#       low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
#       high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
#       group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                              ifelse(target_values >= high_cutoff, "High", NA))
#     } else if (length(target_group) == 1 && target_group < 1) {
#       low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
#       high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
#       group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                              ifelse(target_values >= high_cutoff, "High", NA))
#     } else {
#       group_labels <- ifelse(target_values < target_group, "Low", "High")
#     }
    
#     # Remove NAs (middle group)
#     keep_cells <- !is.na(group_labels)
#     expr_mat <- expr_mat[, keep_cells]
#     meta.data <- meta.data[keep_cells, ]
#     target_binary <- factor(group_labels[keep_cells])
    
#   } else {
#     # Categorical variable
#     if (!is.null(target_group)) {
#       # Select specific groups
#       keep_cells <- target_values %in% target_group
#       expr_mat <- expr_mat[, keep_cells]
#       meta.data <- meta.data[keep_cells, ]
#       target_binary <- factor(target_values[keep_cells])
#     } else {
#       target_binary <- factor(target_values)
#     }
#   }
  
#   if (length(unique(target_binary)) < 2) {
#     stop("Target variable must have at least 2 groups after processing")
#   }
  
#   n_groups <- length(unique(target_binary))
  
#   # ===
#   # 3. Filter and preprocess genes
#   # ===
  
#   # Filter low-expressed genes
#   n_cells_expr <- rowSums(expr_mat > 0)
#   pct_cells_expr <- n_cells_expr / ncol(expr_mat)
#   keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
  
#   expr_mat <- expr_mat[keep_genes, ]
  
#   if (nrow(expr_mat) == 0) {
#     stop("No genes pass filtering criteria")
#   }
  
#   # Preprocessing
#   if (preprocess) {
#     # Log-normalize if data appears to be raw counts
#     if (max(expr_mat) > 100) {
#       expr_mat <- log1p(expr_mat)
#     }
#     # Scale genes (z-score)
#     if (method %in% c("lasso", "gam", "pca_loadings")) {
#       gene_means <- rowMeans(expr_mat)
#       gene_sds <- apply(expr_mat, 1, sd)
#       gene_sds[gene_sds == 0] <- 1  # avoid division by zero
#       expr_mat <- (expr_mat - gene_means) / gene_sds
#     }
#   }
  
#   # ===
#   # 4. Method-specific signature discovery
#   # ===
  
#   result <- switch(method,
                   
#                    # ---
#                    # Random Forest / Tree-based methods
#                    # ---
#                    tree_based = {
#                      if (!requireNamespace("randomForest", quietly = TRUE)) {
#                        stop("randomForest package required. Install with: install.packages('randomForest')")
#                      }
                     
#                      # Transpose for modeling (samples as rows)
#                      X <- t(expr_mat)
#                      y <- target_binary
                     
#                      # Subsample genes if too many (RF can be slow)
#                      if (nrow(expr_mat) > 2000) {
#                        # Pre-filter by variance
#                        gene_vars <- apply(expr_mat, 1, var)
#                        top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
#                        X <- X[, top_var_genes]
#                      }
                     
#                      # Train Random Forest
#                      rf_model <- randomForest::randomForest(
#                        x = X, y = y,
#                        ntree = 500,
#                        importance = TRUE,
#                        ...
#                      )
                     
#                      # Extract feature importance
#                      importance_scores <- randomForest::importance(rf_model)
#                      if (n_groups == 2) {
#                        weights <- importance_scores[, "MeanDecreaseGini"]
#                      } else {
#                        weights <- rowMeans(importance_scores[, grep("MeanDecreaseGini", 
#                                                                     colnames(importance_scores))])
#                      }
                     
#                      # Select top features
#                      top_genes <- names(sort(weights, decreasing=TRUE)[1:min(n_features, length(weights))])
#                      weights <- weights[top_genes]
                     
#                      # Calculate signature scores
#                      scores <- as.numeric(X[, top_genes] %*% weights)
#                      names(scores) <- rownames(X)
                     
#                      # Performance metrics
#                      pred <- rf_model$predicted
#                      if (n_groups == 2) {
#                        if (requireNamespace("pROC", quietly = TRUE)) {
#                          roc_obj <- pROC::roc(y, scores, quiet=TRUE)
#                          auc <- as.numeric(pROC::auc(roc_obj))
#                        } else {
#                          auc <- NA
#                        }
#                        acc <- mean(pred == y)
#                        perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
#                      } else {
#                        acc <- mean(pred == y)
#                        perf <- list(accuracy = acc, confusion = table(pred, y))
#                      }
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = if(return_model) rf_model else NULL)
#                    },
                   
#                    # ---
#                    # LASSO regression
#                    # ---
#                    lasso = {
#                      if (!requireNamespace("glmnet", quietly = TRUE)) {
#                        stop("glmnet package required. Install with: install.packages('glmnet')")
#                      }
                     
#                      X <- t(expr_mat)
#                      y <- target_binary
                     
#                      # Fit LASSO with cross-validation
#                      if (n_groups == 2) {
#                        cv_fit <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1, ...)
#                      } else {
#                        cv_fit <- glmnet::cv.glmnet(X, y, family="multinomial", alpha=1, ...)
#                      }
                     
#                      # Extract coefficients at lambda.1se (more regularized)
#                      coefs <- coef(cv_fit, s="lambda.1se")
                     
#                      if (n_groups == 2) {
#                        weights <- as.numeric(coefs[-1])  # remove intercept
#                        names(weights) <- rownames(coefs)[-1]
#                      } else {
#                        # Average coefficients across classes
#                        coef_list <- lapply(coefs, function(x) as.numeric(x[-1]))
#                        weights <- rowMeans(do.call(cbind, coef_list))
#                        names(weights) <- rownames(coefs[[1]])[-1]
#                      }
                     
#                      # Select non-zero coefficients
#                      nonzero_genes <- names(weights)[weights != 0]
#                      if (length(nonzero_genes) == 0) {
#                        warning("No genes selected by LASSO. Returning top genes by coefficient magnitude.")
#                        top_genes <- names(sort(abs(weights), decreasing=TRUE)[1:min(n_features, length(weights))])
#                      } else {
#                        top_genes <- nonzero_genes[1:min(n_features, length(nonzero_genes))]
#                      }
#                      weights <- weights[top_genes]
                     
#                      # Signature scores
#                      scores <- as.numeric(X[, top_genes] %*% weights)
#                      names(scores) <- rownames(X)
                     
#                      # Performance
#                      pred_probs <- predict(cv_fit, newx=X, s="lambda.1se", type="response")
#                      if (n_groups == 2) {
#                        pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
#                        if (requireNamespace("pROC", quietly = TRUE)) {
#                          roc_obj <- pROC::roc(y, as.numeric(pred_probs), quiet=TRUE)
#                          auc <- as.numeric(pROC::auc(roc_obj))
#                        } else {
#                          auc <- NA
#                        }
#                        acc <- mean(pred == y)
#                        perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
#                      } else {
#                        pred <- levels(y)[apply(pred_probs, 1, which.max)]
#                        acc <- mean(pred == y)
#                        perf <- list(accuracy = acc, confusion = table(pred, y))
#                      }
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = if(return_model) cv_fit else NULL)
#                    },
                   
#                    # ---
#                    # Differential Expression (limma/wilcoxon)
#                    # ---
#                    limma = {
#                      if (!requireNamespace("limma", quietly = TRUE)) {
#                        stop("limma package required. Install with: BiocManager::install('limma')")
#                      }
                     
#                      design <- model.matrix(~0 + target_binary)
#                      colnames(design) <- levels(target_binary)
                     
#                      # Fit linear model
#                      fit <- limma::lmFit(expr_mat, design)
                     
#                      # Define contrasts
#                      if (n_groups == 2) {
#                        contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
#                      } else {
#                        # All pairwise comparisons
#                        contrast_str <- limma::makeContrasts(
#                          contrasts = combn(levels(target_binary), 2, function(x) paste(x, collapse="-")),
#                          levels = design
#                        )
#                      }
                     
#                      contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
#                      fit2 <- limma::contrasts.fit(fit, contrast_mat)
#                      fit2 <- limma::eBayes(fit2)
                     
#                      # Extract results
#                      top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
                     
#                      # Weights as moderated t-statistics
#                      weights <- top_table$t
#                      names(weights) <- rownames(top_table)
                     
#                      top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
#                      weights <- weights[top_genes]
                     
#                      # Signature scores
#                      scores <- colSums(expr_mat[top_genes, ] * weights)
                     
#                      # Performance (simple)
#                      score_cutoff <- median(scores)
#                      pred <- factor(ifelse(scores > score_cutoff, levels(target_binary)[2], 
#                                            levels(target_binary)[1]), levels=levels(target_binary))
#                      acc <- mean(pred == target_binary)
#                      perf <- list(accuracy = acc, top_table = top_table[top_genes, ])
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = if(return_model) fit2 else NULL)
#                    },
                   
#                    wilcoxon = {
#                      # Fast Wilcoxon rank-sum test
#                      pvals <- numeric(nrow(expr_mat))
#                      effect_sizes <- numeric(nrow(expr_mat))
                     
#                      for (i in 1:nrow(expr_mat)) {
#                        if (n_groups == 2) {
#                          group1 <- expr_mat[i, target_binary == levels(target_binary)[1]]
#                          group2 <- expr_mat[i, target_binary == levels(target_binary)[2]]
#                          test <- wilcox.test(group1, group2)
#                          pvals[i] <- test$p.value
#                          effect_sizes[i] <- median(group2) - median(group1)
#                        } else {
#                          test <- kruskal.test(expr_mat[i, ] ~ target_binary)
#                          pvals[i] <- test$p.value
#                          effect_sizes[i] <- var(tapply(expr_mat[i, ], target_binary, median))
#                        }
#                      }
                     
#                      names(pvals) <- rownames(expr_mat)
#                      names(effect_sizes) <- rownames(expr_mat)
                     
#                      # Adjust p-values
#                      padj <- p.adjust(pvals, method="BH")
                     
#                      # Select by p-value and effect size
#                      ranks <- rank(pvals) + rank(-abs(effect_sizes))
#                      top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
                     
#                      weights <- effect_sizes[top_genes]
#                      scores <- colSums(expr_mat[top_genes, ] * weights)
                     
#                      perf <- list(pvalues = pvals[top_genes], padj = padj[top_genes],
#                                   effect_sizes = effect_sizes[top_genes])
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = NULL)
#                    },
                   
#                    # ---
#                    # NMF
#                    # ---
#                    nmf = {
#                      if (!requireNamespace("NMF", quietly = TRUE)) {
#                        stop("NMF package required. Install with: install.packages('NMF')")
#                      }
                     
#                      # Ensure non-negative data
#                      expr_mat_pos <- expr_mat - min(expr_mat) + 0.01
                     
#                      # Run NMF
#                      rank <- min(n_groups + 2, 10)  # adaptive rank
                     
#                      # The 'seed' argument from find_gene_signature conflicts with an NMF internal generic.
#                      # We remove it and rely on the set.seed() call at the start of the function.
#                      dots <- list(...)
#                      dots$seed <- NULL
#                      nmf_res <- do.call(NMF::nmf, c(list(x = expr_mat_pos, rank = rank), dots))

#                      W <- NMF::basis(nmf_res)  # gene loadings
#                      H <- NMF::coef(nmf_res)   # cell scores
                     
#                      # Find component most associated with target
#                      component_cors <- numeric(rank)
#                      for (k in 1:rank) {
#                        if (n_groups == 2) {
#                          component_cors[k] <- abs(cor(H[k, ], as.numeric(target_binary)))
#                        } else {
#                          # ANOVA F-statistic
#                          component_cors[k] <- summary(aov(H[k, ] ~ target_binary))[[1]][1, "F value"]
#                        }
#                      }
                     
#                      best_component <- which.max(component_cors)
#                      weights <- W[, best_component]
#                      names(weights) <- rownames(expr_mat)
                     
#                      top_genes <- names(sort(weights, decreasing=TRUE)[1:min(n_features, length(weights))])
#                      weights <- weights[top_genes]
                     
#                      scores <- H[best_component, ]
#                      names(scores) <- colnames(expr_mat)
                     
#                      perf <- list(component = best_component, correlation = component_cors[best_component])
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = if(return_model) nmf_res else NULL)
#                    },
                   
#                    # ---
#                    # GAM (Generalized Additive Model)
#                    # ---
#                    gam = {
#                      if (!requireNamespace("mgcv", quietly = TRUE)) {
#                        stop("mgcv package required. Install with: install.packages('mgcv')")
#                      }
                     
#                      # Fit GAM for each gene
#                      deviance_explained <- numeric(nrow(expr_mat))
                     
#                      for (i in 1:min(500, nrow(expr_mat))) {  # limit for speed
#                        gene_expr <- expr_mat[i, ]
#                        if (n_groups == 2) {
#                          gam_fit <- mgcv::gam(as.numeric(target_binary) ~ s(gene_expr), 
#                                               family="binomial")
#                        } else {
#                          gam_fit <- mgcv::gam(as.numeric(target_binary) ~ s(gene_expr))
#                        }
#                        deviance_explained[i] <- summary(gam_fit)$dev.expl
#                      }
                     
#                      names(deviance_explained) <- rownames(expr_mat)[1:length(deviance_explained)]
                     
#                      top_genes <- names(sort(deviance_explained, decreasing=TRUE)[1:min(n_features, 
#                                                                                         sum(deviance_explained > 0))])
#                      weights <- deviance_explained[top_genes]
                     
#                      scores <- colSums(expr_mat[top_genes, ] * weights)
                     
#                      perf <- list(deviance_explained = deviance_explained[top_genes])
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = NULL)
#                    },
                   
#                    # ---
#                    # PCA loadings
#                    # ---
#                    pca_loadings = {
#                      # Run PCA
#                      pca_res <- prcomp(t(expr_mat), center=FALSE, scale.=FALSE)
                     
#                      # Find PC most correlated with target
#                      pc_cors <- numeric(min(50, ncol(pca_res$x)))
#                      for (k in 1:length(pc_cors)) {
#                        if (n_groups == 2) {
#                          pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
#                        } else {
#                          pc_cors[k] <- summary(aov(pca_res$x[, k] ~ target_binary))[[1]][1, "F value"]
#                        }
#                      }
                     
#                      best_pc <- which.max(pc_cors)
#                      weights <- pca_res$rotation[, best_pc]
                     
#                      top_genes <- names(sort(abs(weights), decreasing=TRUE)[1:min(n_features, length(weights))])
#                      weights <- weights[top_genes]
                     
#                      scores <- pca_res$x[, best_pc]
                     
#                      perf <- list(PC = best_pc, correlation = pc_cors[best_pc],
#                                   variance_explained = summary(pca_res)$importance[2, best_pc])
                     
#                      list(genes = top_genes, weights = weights, scores = scores,
#                           performance = perf, model = if(return_model) pca_res else NULL)
#                    }
#   )
  
#   # ===
#   # 5. Return results
#   # ===
  
#   result$method <- method
#   result$target_var <- target_var
#   result$n_groups <- n_groups
#   result$n_cells <- ncol(expr_mat)
  
#   class(result) <- c("gene_signature", "list")
#   return(result)
# }
# #' Find Gene Signatures that Best Separate Target Variable
# #'
# #' @param data Seurat object, count matrix, or data.frame with genes as rows/columns
# #' @param meta.data Optional metadata data.frame (required if data is not Seurat)
# #' @param target_var Column name in metadata representing the target variable
# #' @param target_group For numeric targets: quantile cutoff (0-1) or list(low=0.25, high=0.75)
# #'                     For factor: specific levels to compare (default: all vs all)
# #' @param method One of: "tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings"
# #' @param n_features Number of top features to return (default: 50)
# #' @param preprocess Logical, whether to normalize/scale data (default: TRUE)
# #' @param min_cells Minimum cells expressing gene (default: 10)
# #' @param min_pct Minimum percentage of cells expressing gene (default: 0.01)
# #' @param return_model Logical, return full model object (default: FALSE)
# #' @param seed Random seed for reproducibility (default: 42)
# #' @param ... Additional method-specific parameters
# #'
# #' @return List containing:
# #'   - genes: Character vector of selected genes
# #'   - weights: Named numeric vector of gene weights/importance
# #'   - scores: Per-cell signature scores (if applicable)
# #'   - performance: Separation metrics (AUC, accuracy, etc.)
# #'   - method: Method used
# #'   - model: Full model object (if return_model=TRUE)
# #'
# #' @examples
# #' # Random Forest for binary classification
# #' result <- find_gene_signature(seurat_obj, target_var="cell_type", 
# #'                                target_group=c("TypeA", "TypeB"), 
# #'                                method="tree_based")
# #' 
# #' # LASSO for continuous variable (top vs bottom quartile)
# #' result <- find_gene_signature(seurat_obj, target_var="pseudotime", 
# #'                                target_group=0.25, method="lasso")
# #' 
# #' # NMF for multi-class
# #' result <- find_gene_signature(count_matrix, meta.data=metadata, 
# #'                                target_var="condition", method="nmf")
# #' @export
# find_gene_signature_v3 <- function(data, 
#                                  meta.data = NULL,
#                                  target_var,
#                                  target_group = NULL,
#                                  method = c("tree_based", "lasso", "limma", 
#                                             "nmf", "wilcoxon", "gam", "pca_loadings"),
#                                  n_features = 50,
#                                  preprocess = TRUE,
#                                  min_cells = 10,
#                                  min_pct = 0.01,
#                                  return_model = FALSE,
#                                  fgs_seed = 42,
#                                  lambda_selection = "lambda.1se",
#                                  ...) {
  
#   # ===
#   # 0. Batch Processing (v2와 동일)
#   # ===
  
#   all_methods <- c("tree_based", "lasso", "limma", "nmf", "wilcoxon", "gam", "pca_loadings")
  
#   if (is.null(method)) {
#     method <- all_methods
#   }
  
#   if (length(method) > 1) {
#     current_call <- match.call()
    
#     results_list <- lapply(method, function(m) {
#       single_call <- current_call
#       single_call$method <- m
#       tryCatch({
#         eval(single_call)
#       }, error = function(e) {
#         warning(sprintf("Method '%s' failed with error: %s", m, e$message))
#         return(list(method = m, error = e$message))
#       })
#     })
    
#     names(results_list) <- method
#     return(results_list)
#   }
  
#   if (!method %in% all_methods) {
#     stop(sprintf("Invalid method '%s'. Choose from: %s", 
#                  method, paste(all_methods, collapse=", ")))
#   }
  
#   set.seed(fgs_seed)
  
#   # ===
#   # 1. Input validation and data extraction (v2와 동일)
#   # ===
  
#   is_seurat <- inherits(data, "Seurat")
  
#   if (is_seurat) {
#     if (!requireNamespace("Seurat", quietly = TRUE)) {
#       stop("Seurat package required but not installed")
#     }
#     if (is.null(meta.data)) {
#       meta.data <- data@meta.data
#     }
#     if ("data" %in% names(data@assays[[1]])) {
#       expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "data"))
#     } else {
#       expr_mat <- as.matrix(Seurat::GetAssayData(data, layer = "counts"))
#     }
#   } else {
#     if (is.null(meta.data)) {
#       stop("meta.data must be provided when data is not a Seurat object")
#     }
#     expr_mat <- as.matrix(data)
#     if (nrow(expr_mat) == nrow(meta.data)) {
#       expr_mat <- t(expr_mat)
#     }
#   }
  
#   if (!target_var %in% colnames(meta.data)) {
#     stop(sprintf("target_var '%s' not found in metadata columns: %s", 
#                  target_var, paste(colnames(meta.data), collapse=", ")))
#   }
  
#   common_cells <- intersect(colnames(expr_mat), rownames(meta.data))
#   if (length(common_cells) == 0) {
#     stop("No common cells found between expression matrix and metadata")
#   }
#   expr_mat <- expr_mat[, common_cells]
#   meta.data <- meta.data[common_cells, ]
  
#   # ===
#   # 2. Process target variable (v2와 동일)
#   # ===
  
#   target_values <- meta.data[[target_var]]
#   target_type <- class(target_values)[1]
  
#   if (is.numeric(target_values)) {
#     if (is.null(target_group)) {
#       target_group <- 0.25
#     }
#     if (is.list(target_group)) {
#       low_cutoff <- quantile(target_values, target_group$low, na.rm=TRUE)
#       high_cutoff <- quantile(target_values, target_group$high, na.rm=TRUE)
#       group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                              ifelse(target_values >= high_cutoff, "High", NA))
#     } else if (length(target_group) == 1 && target_group < 1) {
#       low_cutoff <- quantile(target_values, target_group, na.rm=TRUE)
#       high_cutoff <- quantile(target_values, 1 - target_group, na.rm=TRUE)
#       group_labels <- ifelse(target_values <= low_cutoff, "Low",
#                              ifelse(target_values >= high_cutoff, "High", NA))
#     } else {
#       group_labels <- ifelse(target_values < target_group, "Low", "High")
#     }
#     keep_cells <- !is.na(group_labels)
#     expr_mat <- expr_mat[, keep_cells]
#     meta.data <- meta.data[keep_cells, ]
#     target_binary <- factor(group_labels[keep_cells])
#   } else {
#     if (!is.null(target_group)) {
#       keep_cells <- target_values %in% target_group
#       expr_mat <- expr_mat[, keep_cells]
#       meta.data <- meta.data[keep_cells, ]
#       target_binary <- factor(target_values[keep_cells])
#     } else {
#       target_binary <- factor(target_values)
#     }
#   }
  
#   if (length(unique(target_binary)) < 2) {
#     stop("Target variable must have at least 2 groups after processing")
#   }
  
#   n_groups <- length(unique(target_binary))
  
#   # ===
#   # 3. Filter and preprocess genes (v2와 동일)
#   # ===
  
#   n_cells_expr <- rowSums(expr_mat > 0)
#   pct_cells_expr <- n_cells_expr / ncol(expr_mat)
#   keep_genes <- (n_cells_expr >= min_cells) & (pct_cells_expr >= min_pct)
#   expr_mat <- expr_mat[keep_genes, ]
  
#   if (nrow(expr_mat) == 0) {
#     stop("No genes pass filtering criteria")
#   }
  
#   if (preprocess) {
#     if (max(expr_mat) > 100) {
#       expr_mat <- log1p(expr_mat)
#     }
#     if (method %in% c("lasso", "gam", "pca_loadings")) {
#       gene_means <- rowMeans(expr_mat)
#       gene_sds <- apply(expr_mat, 1, sd)
#       gene_sds[gene_sds == 0] <- 1
#       expr_mat <- (expr_mat - gene_means) / gene_sds
#     }
#   }
  
#   # ===
#   # 4. Method-specific signature discovery
#   # ===
  
#   # X (samples x genes)는 여러 메서드에서 공통으로 사용
#   X <- t(expr_mat)
#   y <- target_binary
  
#   result <- switch(method,
                 
#                  tree_based = {
#                    if (!requireNamespace("randomForest", quietly = TRUE)) {
#                      stop("randomForest package required. Install with: install.packages('randomForest')")
#                    }
                   
#                    X_rf <- X
#                    if (nrow(expr_mat) > 2000) {
#                      gene_vars <- apply(expr_mat, 1, var)
#                      top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:2000])
#                      X_rf <- X[, top_var_genes]
#                    }
                   
#                    rf_model <- randomForest::randomForest(x = X_rf, y = y, ntree = 500, importance = TRUE, ...)
                   
#                    importance_scores <- randomForest::importance(rf_model)
#                    if (n_groups == 2) {
#                      weights_magnitude <- importance_scores[, "MeanDecreaseGini"]
#                    } else {
#                      weights_magnitude <- rowMeans(importance_scores[, grep("MeanDecreaseGini", colnames(importance_scores))])
#                    }
                   
#                    top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, length(weights_magnitude))])
#                    weights_magnitude <- weights_magnitude[top_genes]
                   
#                    # === v3: 방향성 보정 ===
#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]
#                      g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      # effect_size: group2 (예: NR) - group1 (예: R)
#                      effect_size <- mean_g2 - mean_g1 
#                      weights <- weights_magnitude * sign(effect_size) # 중요도 * 방향
#                    } else {
#                      warning("tree_based: n_groups > 2. Score represents magnitude (importance), not direction.")
#                      weights <- weights_magnitude
#                    }
#                    # === v3 끝 ===
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    pred <- rf_model$predicted
#                    if (n_groups == 2) {
#                      if (requireNamespace("pROC", quietly = TRUE)) {
#                        roc_obj <- pROC::roc(y, scores, quiet=TRUE)
#                        auc <- as.numeric(pROC::auc(roc_obj))
#                      } else { auc <- NA }
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
#                    } else {
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, confusion = table(pred, y))
#                    }
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) rf_model else NULL)
#                  },
                 
#                  lasso = {
#                    if (!requireNamespace("glmnet", quietly = TRUE)) {
#                      stop("glmnet package required. Install with: install.packages('glmnet')")
#                    }
                   
#                    if (n_groups == 2) {
#                      cv_fit <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1, ...)
#                    } else {
#                      cv_fit <- glmnet::cv.glmnet(X, y, family="multinomial", alpha=1, ...)
#                    }
                   
#                    coefs <- coef(cv_fit, s = lambda_selection) 
                   
#                    if (n_groups == 2) {
#                      weights_all <- as.numeric(coefs[-1])
#                      names(weights_all) <- rownames(coefs)[-1]
#                    } else {
#                      coef_list <- lapply(coefs, function(x) as.numeric(x[-1]))
#                      weights_all <- rowMeans(do.call(cbind, coef_list))
#                      names(weights_all) <- rownames(coefs[[1]])[-1]
#                    }
                   
#                    nonzero_genes <- names(weights_all)[weights_all != 0]
#                    if (length(nonzero_genes) == 0) {
#                      warning(sprintf("No genes selected by LASSO with s='%s'. Returning top genes by coefficient magnitude.", lambda_selection))
#                      top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
#                    } else {
#                      # v3 수정: non-zero 유전자 중에서도 절대값 크기 순으로 정렬
#                      nonzero_weights <- weights_all[nonzero_genes]
#                      top_genes <- names(sort(abs(nonzero_weights), decreasing=TRUE)[1:min(n_features, length(nonzero_weights))])
#                    }
                   
#                    weights <- weights_all[top_genes] # 최종 가중치 (부호 있음)
                   
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    pred_probs <- predict(cv_fit, newx=X, s = lambda_selection, type="response")
#                    if (n_groups == 2) {
#                      pred <- factor(ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]), levels=levels(y))
#                      if (requireNamespace("pROC", quietly = TRUE)) {
#                        # v3 수정: AUC 계산 시 점수(scores) 사용
#                        roc_obj <- pROC::roc(y, scores, quiet=TRUE) 
#                        auc <- as.numeric(pROC::auc(roc_obj))
#                      } else { auc <- NA }
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, auc = auc, confusion = table(pred, y))
#                    } else {
#                      pred <- levels(y)[apply(pred_probs, 1, which.max)]
#                      acc <- mean(pred == y)
#                      perf <- list(accuracy = acc, confusion = table(pred, y))
#                    }
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) cv_fit else NULL)
#                  },
                 
#                  limma = {
#                    if (!requireNamespace("limma", quietly = TRUE)) {
#                      stop("limma package required. Install with: BiocManager::install('limma')")
#                    }
#                    design <- model.matrix(~0 + target_binary)
#                    colnames(design) <- levels(target_binary)
#                    fit <- limma::lmFit(expr_mat, design)
#                    if (n_groups == 2) {
#                      contrast_str <- paste(levels(target_binary)[2], levels(target_binary)[1], sep="-")
#                    } else {
#                      contrast_str <- limma::makeContrasts(
#                        contrasts = combn(levels(target_binary), 2, function(x) paste(x, collapse="-")),
#                        levels = design
#                      )
#                    }
#                    contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
#                    fit2 <- limma::contrasts.fit(fit, contrast_mat)
#                    fit2 <- limma::eBayes(fit2)
#                    top_table <- limma::topTable(fit2, number=Inf, sort.by="B")
                   
#                    weights_all <- top_table$t
#                    names(weights_all) <- rownames(top_table)
                   
#                    top_genes <- rownames(top_table)[1:min(n_features, nrow(top_table))]
#                    weights <- weights_all[top_genes] # 최종 가중치 (부호 있음)
                   
#                    # v3 수정: 점수 계산 방식 통일
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    # ... (limma의 performance 부분은 v2와 동일하게 유지)
#                    perf <- list(top_table = top_table[top_genes, ])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) fit2 else NULL)
#                  },
                 
#                  wilcoxon = {
#                    pvals <- numeric(nrow(expr_mat))
#                    effect_sizes <- numeric(nrow(expr_mat))
                   
#                    for (i in 1:nrow(expr_mat)) {
#                      if (n_groups == 2) {
#                        group1 <- expr_mat[i, target_binary == levels(target_binary)[1]]
#                        group2 <- expr_mat[i, target_binary == levels(target_binary)[2]]
#                        test <- wilcox.test(group1, group2)
#                        pvals[i] <- test$p.value
#                        effect_sizes[i] <- median(group2) - median(group1) # 부호 있음
#                      } else {
#                        test <- kruskal.test(expr_mat[i, ] ~ target_binary)
#                        pvals[i] <- test$p.value
#                        effect_sizes[i] <- var(tapply(expr_mat[i, ], target_binary, median)) # 부호 없음
#                      }
#                    }
                   
#                    names(pvals) <- rownames(expr_mat)
#                    names(effect_sizes) <- rownames(expr_mat)
                   
#                    padj <- p.adjust(pvals, method="BH")
                   
#                    # v3 수정: 랭킹 방식을 p-value와 effect-size의 절대값으로
#                    ranks <- rank(pvals) + rank(-abs(effect_sizes))
#                    top_genes <- names(sort(ranks)[1:min(n_features, length(ranks))])
                   
#                    weights <- effect_sizes[top_genes] # 최종 가중치 (부호 있음)
                   
#                    # v3 수정: 점수 계산 방식 통일
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(pvalues = pvals[top_genes], padj = padj[top_genes],
#                                 effect_sizes = effect_sizes[top_genes])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = NULL)
#                  },
                 
#                  nmf = {
#                    if (!requireNamespace("NMF", quietly = TRUE)) {
#                      stop("NMF package required. Install with: install.packages('NMF')")
#                    }
#                    expr_mat_pos <- expr_mat - min(expr_mat) + 0.01
#                    rank <- min(n_groups + 2, 10)
#                    dots <- list(...)
#                    dots$seed <- NULL 
#                    nmf_res <- do.call(NMF::nmf, c(list(x = expr_mat_pos, rank = rank), dots))
#                    W <- NMF::basis(nmf_res)
#                    H <- NMF::coef(nmf_res)
                   
#                    component_cors <- numeric(rank)
#                    for (k in 1:rank) {
#                      if (n_groups == 2) {
#                        component_cors[k] <- abs(cor(H[k, ], as.numeric(target_binary)))
#                      } else {
#                        component_cors[k] <- summary(aov(H[k, ] ~ target_binary))[[1]][1, "F value"]
#                      }
#                    }
                   
#                    best_component <- which.max(component_cors)
#                    weights_magnitude <- W[, best_component] # 항상 양수
#                    names(weights_magnitude) <- rownames(expr_mat)
                   
#                    top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, length(weights_magnitude))])
#                    weights_magnitude <- weights_magnitude[top_genes]
                   
#                    # === v3: 방향성 보정 ===
#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]
#                      g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      effect_size <- mean_g2 - mean_g1
#                      weights <- weights_magnitude * sign(effect_size) # 중요도 * 방향
#                    } else {
#                      warning("nmf: n_groups > 2. Score represents magnitude (importance), not direction.")
#                      weights <- weights_magnitude
#                    }
#                    # === v3 끝 ===
                   
#                    # v3 수정: 점수 계산 방식 통일
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(component = best_component, correlation = component_cors[best_component])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) nmf_res else NULL)
#                  },
                 
#                  gam = {
#                    if (!requireNamespace("mgcv", quietly = TRUE)) {
#                      stop("mgcv package required. Install with: install.packages('mgcv')")
#                    }
                   
#                    # v3 수정: 더 많은 유전자를 테스트하기 위해 상한을 높임 (속도 주의)
#                    n_test_genes <- min(2000, nrow(expr_mat)) 
#                    test_genes_idx <- 1:n_test_genes
                   
#                    # 상위 2000개 유전자가 아닌, 분산이 높은 2000개 유전자를 테스트
#                    if(nrow(expr_mat) > 2000) {
#                        gene_vars <- apply(expr_mat, 1, var)
#                        test_genes_idx <- order(gene_vars, decreasing = TRUE)[1:2000]
#                    }
                   
#                    deviance_explained <- numeric(nrow(expr_mat))
#                    names(deviance_explained) <- rownames(expr_mat)
                   
#                    for (i in test_genes_idx) {
#                      gene_expr <- expr_mat[i, ]
#                      if (n_groups == 2) {
#                        y_numeric_gam <- as.numeric(target_binary) - 1 
#                        gam_fit <- mgcv::gam(y_numeric_gam ~ s(gene_expr), family="binomial", ...)
#                      } else {
#                        gam_fit <- mgcv::gam(as.numeric(target_binary) ~ s(gene_expr), ...)
#                      }
#                      deviance_explained[i] <- summary(gam_fit)$dev.expl
#                    }
                   
#                    weights_magnitude <- deviance_explained
#                    top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0))])
#                    weights_magnitude <- weights_magnitude[top_genes]

#                    # === v3: 방향성 보정 ===
#                    if (n_groups == 2) {
#                      g1_cells <- y == levels(y)[1]
#                      g2_cells <- y == levels(y)[2]
#                      mean_g1 <- colMeans(X[g1_cells, top_genes, drop=FALSE])
#                      mean_g2 <- colMeans(X[g2_cells, top_genes, drop=FALSE])
#                      effect_size <- mean_g2 - mean_g1
#                      weights <- weights_magnitude * sign(effect_size) # 중요도 * 방향
#                    } else {
#                      warning("gam: n_groups > 2. Score represents magnitude (importance), not direction.")
#                      weights <- weights_magnitude
#                    }
#                    # === v3 끝 ===
                   
#                    # v3 수정: 점수 계산 방식 통일
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(deviance_explained = weights_magnitude) # 방향성 보정 전의 중요도
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = NULL)
#                  },
                 
#                  pca_loadings = {
#                    # v3 수정: PCA는 스케일링된 데이터를 사용해야 함
#                    # 함수 상단의 전처리 블록에서 method="pca_loadings"일 때 이미 스케일링됨.
#                    pca_res <- prcomp(X, center=FALSE, scale.=FALSE) 
                   
#                    pc_cors <- numeric(min(50, ncol(pca_res$x)))
#                    for (k in 1:length(pc_cors)) {
#                      if (n_groups == 2) {
#                        pc_cors[k] <- abs(cor(pca_res$x[, k], as.numeric(target_binary)))
#                      } else {
#                        pc_cors[k] <- summary(aov(pca_res$x[, k] ~ target_binary))[[1]][1, "F value"]
#                      }
#                    }
                   
#                    best_pc <- which.max(pc_cors)
#                    weights_all <- pca_res$rotation[, best_pc] # 부호 있음
                   
#                    top_genes <- names(sort(abs(weights_all), decreasing=TRUE)[1:min(n_features, length(weights_all))])
#                    weights <- weights_all[top_genes] # 최종 가중치 (부호 있음)
                   
#                    # v3 수정: 점수 계산 방식 통일
#                    # (v2의 scores = pca_res$x[, best_pc]는 '모든' 유전자를 사용한 점수였음)
#                    scores <- as.numeric(X[, top_genes] %*% weights)
#                    names(scores) <- rownames(X)
                   
#                    perf <- list(PC = best_pc, correlation = pc_cors[best_pc],
#                                 variance_explained = summary(pca_res)$importance[2, best_pc])
                   
#                    list(genes = top_genes, weights = weights, scores = scores,
#                         performance = perf, model = if(return_model) pca_res else NULL)
#                  }
#   )
  
#   # ===
#   # 5. Return results (v2와 동일)
#   # ===
  
#   result$method <- method
#   result$target_var <- target_var
#   result$n_groups <- n_groups
#   result$n_cells <- ncol(expr_mat)
#   result$formula <- paste(deparse(match.call()), collapse = " ")
  
#   class(result) <- c("gene_signature", "list")
#   return(result)
# }


#' @export
find_gene_signature_v5.3_test_to_delete <- function(data, 
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
                                 lambda_selection = "lambda.1se",
                                 enet.alpha = 0.5,
                                 pca.n_pcs = 1,
                                 gam.min_unique = 15,
                                 gam.k = NULL,  # v5.3: NULL이면 동적 k 사용
                                 gam.k_dynamic_factor = 5,  # v5.3: 동적 k 계산 시 사용 (n_unique_vals / factor)
                                 ...) {
  
  # v5.3은 v5.2와 동일하지만 gam 부분만 동적 k를 사용
  # v5.2를 호출하되, gam 메서드일 때만 동적 k를 적용
  
  all_methods <- c(
    "random_forest", "random_forest_ranger", "xgboost",
    "lasso", "ridge", "elastic_net",
    "pca_loadings", "nmf_loadings",
    "gam", "limma", "wilcoxon"
  )
  
  if (is.null(method)) {
    method <- all_methods
  }
  
  # gam 메서드가 포함되어 있고 gam.k가 NULL이면 동적 k를 사용
  use_dynamic_k <- "gam" %in% method && is.null(gam.k)
  
  if (use_dynamic_k) {
    # v5.3: gam에 동적 k 적용을 위해 직접 구현
    # v5.2의 전처리 함수 재사용
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
    
    for (m in method) {
      if (!m %in% all_methods) {
        warning(sprintf("Invalid method '%s'. Skipping.", m))
        next
      }
      
      message(sprintf("--- Running Method: %s ---", m))
      
      tryCatch({
        if (m %in% c("limma", "wilcoxon", "nmf_loadings", "random_forest", "random_forest_ranger")) {
          expr_mat_method <- expr_mat_base
        } else if (m %in% c("lasso", "ridge", "elastic_net", "gam", "pca_loadings", "xgboost")) {
          expr_mat_method <- if(is.null(preprocessed_data$expr_mat_scaled)) expr_mat_base else preprocessed_data$expr_mat_scaled
        }
        
        if (m %in% c("wilcoxon", "pca_loadings")) {
          if (!is.null(preprocessed_data$expr_mat_corrected)) {
            expr_mat_method <- preprocessed_data$expr_mat_corrected
          }
        }
        
        X <- t(expr_mat_method)
        y <- target_binary
        
        if (m == "gam") {
          # v5.3: 동적 k 적용
          if (!requireNamespace("mgcv", quietly = TRUE)) {
            stop("mgcv package required.")
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
            
            # v5.3: unique value threshold 먼저 적용
            if (n_unique_vals < gam.min_unique) {
              deviance_explained[i] <- 0
              genes_skipped <- genes_skipped + 1
              next
            }
            
            # v5.3: 동적 k 계산
            k_dynamic <- max(3, min(10, floor(n_unique_vals / gam.k_dynamic_factor)))
            if (k_dynamic != 10) genes_with_dynamic_k <- genes_with_dynamic_k + 1
            
            full_formula_str <- paste("y_var_numeric ~ s(gene_expr, k=", k_dynamic, ", bs='cr')", formula_base)
            
            fit_result <- tryCatch({
              if (n_groups == 2) {
                mgcv::bam(as.formula(full_formula_str), data = gam_data, family="binomial", ...)
              } else {
                mgcv::bam(as.formula(full_formula_str), data = gam_data, ...)
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
          if (genes_with_dynamic_k > 0) {
            message(sprintf("GAM: Applied dynamic k for %d genes (k < 10).", genes_with_dynamic_k))
          }
          
          weights_magnitude <- deviance_explained
          top_genes <- names(sort(weights_magnitude, decreasing=TRUE)[1:min(n_features, sum(weights_magnitude > 0, na.rm=TRUE))])
          
          if (length(top_genes) == 0) {
            warning("GAM method found no genes with deviance > 0.")
            result <- list(genes=character(0), weights=numeric(0), scores=numeric(0), performance=list())
          } else {
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
            
            result <- list(genes = top_genes, weights = weights, scores = scores,
                         performance = perf, model = NULL)
          }
        } else {
          # gam이 아닌 다른 메서드는 v5.2를 호출
          v52_call <- match.call()
          v52_call[[1]] <- as.name("find_gene_signature_v5.2")
          v52_call$gam.k <- if (is.null(gam.k)) 10 else gam.k
          v52_call$gam.k_dynamic_factor <- NULL
          v52_call$method <- m
          result <- eval(v52_call)
          result <- result[[m]]
        }
        
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
    }
    
    return(results_list)
  } else {
    # gam.k가 지정되어 있으면 v5.2를 그대로 사용
    v52_call <- match.call()
    v52_call[[1]] <- as.name("find_gene_signature_v5.2")
    v52_call$gam.k_dynamic_factor <- NULL
    return(eval(v52_call))
  }
}

# # ============================================================================
# # DEPRECATED: FGS_v5.2 and FGS_v5.3 - Moved to test_to_delete.R
# # Use FGS() instead, which is an alias for find_gene_signature_v5.3
# # ============================================================================

# #' Find Gene Signature v5.2 (FGS_v5.2) - DEPRECATED
# #'
# #' @description
# #' DEPRECATED: Use \code{FGS()} instead, which is an alias for \code{find_gene_signature_v5.3}.
# #' This function is kept for backward compatibility only.
# #'
# #' @param ... All arguments passed to \code{find_gene_signature_v5.2}
# #' @return See \code{find_gene_signature_v5.2}
# #' @seealso \code{\link{FGS}}, \code{\link{find_gene_signature_v5.2}}
# #' @export
# FGS_v5.2 <- function(...) {
#   .Deprecated("FGS", msg = "FGS_v5.2 is deprecated. Use FGS() instead, which uses find_gene_signature_v5.3.")
#   find_gene_signature_v5.2(...)
# }

# #' Find Gene Signature v5.3 (FGS_v5.3) - DEPRECATED
# #'
# #' @description
# #' DEPRECATED: Use \code{FGS()} instead, which is an alias for \code{find_gene_signature_v5.3}.
# #' This function is kept for backward compatibility only.
# #'
# #' @param ... All arguments passed to \code{find_gene_signature_v5.3}
# #' @return See \code{find_gene_signature_v5.3}
# #' @seealso \code{\link{FGS}}, \code{\link{find_gene_signature_v5.3}}
# #' @export
# FGS_v5.3 <- function(...) {
#   .Deprecated("FGS", msg = "FGS_v5.3 is deprecated. Use FGS() instead.")
#   find_gene_signature_v5.3(...)
# }



# test.R


# analysis -----
NEBULA_test=function(){}


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
#' @title MAST 파이프라인 함수 (DEPRECATED)
#' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMAST_v1을 사용하세요.
#' @export
runMAST_test <- function(...) {
  .Deprecated("runMAST_v1", package = "myR", msg = "runMAST in test.R is deprecated. Use runMAST_v1 from test_analysis.R instead.")
  if (exists("runMAST_v1", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runMAST_v1", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runMAST_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

# [V4] find_gene_signature_v4
#
# 변경점:
# 1. control_vars = NULL 인자 추가 (교란 변수 보정)
# 2. test_n = NULL 인자 추가 (속도 향상을 위한 사전 필터링)
# 3. limma, lasso, tree_based, gam: Native 방식 보정
# 4. wilcoxon, nmf, pca_loadings: limma::removeBatchEffect 방식 보정

#' @export
find_gene_signature_v4_test <- function(data, 
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
find_gene_signature_v4.1_test <- function(data, 
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
find_gene_signature_v4.2_test <- function(data, 
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
fgs_preprocess_data_test <- function(data, 
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
find_gene_signature_v5_test <- function(data, 
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
  covariate_names <- if (is.null(covariate_mat_model)) character(0) else colnames(covariate_mat_model)

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
      gene_feature_names <- colnames(X)
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
      result$fgs_version <- method_impl
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
fgs_preprocess_data_v5.2_test <- function(data, 
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
find_gene_signature_v5.2_test <- function(data, 
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
                                 method_impl = c("v5.2","v5.4"),
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
          
          gene_names_in_model <- colnames(X_rf)[!colnames(X_rf) %in% covariate_names]
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
          build_ranger_model <- function(importance_mode = "impurity") {
            ranger::ranger(
              y ~ .,
              data = ranger_data,
              num.trees = 500,
              importance = importance_mode,
              respect.unordered.factors = if (method_impl == "v5.4") "order" else NULL,
              ...
            )
          }
          
          importance_mode <- if (method_impl == "v5.4") "permutation" else "impurity"
          rf_model <- build_ranger_model(importance_mode)
          
          weights_magnitude_all <- ranger::importance(rf_model)
          if (method_impl == "v5.4" && (is.null(weights_magnitude_all) || all(!is.finite(weights_magnitude_all)))) {
            warning("ranger: permutation importance returned no usable values; falling back to impurity-based importance.")
            rf_model <- build_ranger_model("impurity")
            weights_magnitude_all <- ranger::importance(rf_model)
          }
          weights_magnitude_all <- weights_magnitude_all[is.finite(weights_magnitude_all)]
          
          gene_names_in_model <- names(weights_magnitude_all)[!names(weights_magnitude_all) %in% covariate_names]
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
          
          gene_names_in_model <- names(weights_magnitude_all)[!names(weights_magnitude_all) %in% covariate_names]
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
          if (method_impl == "v5.4") {
            X_model <- as.matrix(X_model)
          }
          
          # alpha=0 (Ridge)
          y_glmnet <- if (method_impl == "v5.4" && n_groups == 2) as.numeric(y == levels(y)[2]) else y
          if (n_groups == 2) {
            cv_fit <- glmnet::cv.glmnet(
              X_model, y_glmnet,
              family = "binomial",
              alpha = 0,
              penalty.factor = penalty_vec,
              ...
            )
          } else {
            cv_fit <- glmnet::cv.glmnet(
              X_model, y,
              family = "multinomial",
              alpha = 0,
              penalty.factor = penalty_vec,
              ...
            )
          }
          
          coefs <- coef(cv_fit, s = lambda_selection)
          if (n_groups == 2) {
            coef_mat <- as.matrix(coefs)
            weights_all <- coef_mat[-1, 1]
            names(weights_all) <- rownames(coef_mat)[-1]
          } else {
            coef_list <- lapply(coefs, function(x) as.matrix(x)[-1, 1])
            weights_all <- rowMeans(do.call(cbind, coef_list))
            names(weights_all) <- rownames(as.matrix(coefs[[1]]))[-1]
          }
          if (!is.null(control_vars)) {
            weights_all <- weights_all[names(weights_all) %in% gene_feature_names]
          }
          weights_all <- weights_all[is.finite(weights_all)]
          
          if (length(weights_all) == 0) {
            stop("Ridge returned no gene coefficients.")
          }
          
          top_genes <- names(sort(abs(weights_all), decreasing = TRUE)[
            1:min(n_features, length(weights_all))
          ])
          weights <- weights_all[top_genes]
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          pred_probs <- predict(cv_fit, newx = X_model, s = lambda_selection, type = "response")
          if (n_groups == 2) {
            pred_probs <- as.vector(pred_probs)
          }
          if (n_groups == 2) {
            pred <- factor(
              ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]),
              levels = levels(y)
            )
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
        },
        
        elastic_net = { # [Req 4] 신규
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
          if (method_impl == "v5.4") {
            X_model <- as.matrix(X_model)
          }
          
          # alpha=enet.alpha (기본 0.5)
          y_glmnet <- if (method_impl == "v5.4" && n_groups == 2) as.numeric(y == levels(y)[2]) else y
          if (n_groups == 2) {
            cv_fit <- glmnet::cv.glmnet(
              X_model, y_glmnet,
              family = "binomial",
              alpha = enet.alpha,
              penalty.factor = penalty_vec,
              ...
            )
          } else {
            cv_fit <- glmnet::cv.glmnet(
              X_model, y,
              family = "multinomial",
              alpha = enet.alpha,
              penalty.factor = penalty_vec,
              ...
            )
          }
          
          coefs <- coef(cv_fit, s = lambda_selection)
          
          if (n_groups == 2) {
            coef_mat <- as.matrix(coefs)
            weights_all <- coef_mat[-1, 1]
            names(weights_all) <- rownames(coef_mat)[-1]
          } else {
            coef_list <- lapply(coefs, function(x) as.matrix(x)[-1, 1])
            weights_all <- rowMeans(do.call(cbind, coef_list))
            names(weights_all) <- rownames(as.matrix(coefs[[1]]))[-1]
          }
          if (!is.null(control_vars)) {
            weights_all <- weights_all[names(weights_all) %in% gene_feature_names]
          }
          weights_all <- weights_all[is.finite(weights_all)]
          
          if (length(weights_all) == 0) {
            stop("Elastic net returned no gene coefficients.")
          }
          
          top_genes <- names(sort(abs(weights_all), decreasing = TRUE)[
            1:min(n_features, length(weights_all))
          ])
          weights <- weights_all[top_genes]
          
          scores <- as.numeric(X[, top_genes] %*% weights)
          names(scores) <- rownames(X)
          
          pred_probs <- predict(cv_fit, newx = X_model, s = lambda_selection, type = "response")
          if (n_groups == 2) {
            pred_probs <- as.vector(pred_probs)
          }
          if (n_groups == 2) {
            pred <- factor(
              ifelse(pred_probs > 0.5, levels(y)[2], levels(y)[1]),
              levels = levels(y)
            )
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
          if (method_impl == "v5.4") {
            expr_mat_nmf <- as.matrix(expr_mat_nmf)
          }
          
          if (!is.null(control_vars)) {
            if (!requireNamespace("limma", quietly = TRUE)) {
              stop("limma package required for NMF batch correction.")
            }
            warning(sprintf("Method NMF: Applying limma::removeBatchEffect for: %s...", m))
            covariate_mat <- model.matrix(~ . - 1, data = meta.data_clean[, control_vars, drop=FALSE])
            expr_mat_nmf <- limma::removeBatchEffect(expr_mat_nmf, covariates = covariate_mat)
          }
          if (method_impl == "v5.4") {
            min_val <- suppressWarnings(min(expr_mat_nmf, na.rm = TRUE))
            shift <- if (is.finite(min_val) && min_val < 0) abs(min_val) + 1e-6 else 1e-6
            expr_mat_pos <- expr_mat_nmf + shift
          } else {
            expr_mat_pos <- expr_mat_nmf - min(expr_mat_nmf) + 0.01 
          }
          
          max_rank_allowed <- max(2, min(ncol(expr_mat_pos) - 1L, 10))
          rank <- max(2, min(n_groups + 2, max_rank_allowed))
          if (!is.finite(rank) || rank < 2) {
            stop("NMF: Unable to determine a valid rank; check input matrix dimensions.")
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
          
          nmf_res <- tryCatch(do.call(NMF::nmf, nmf_args), error = function(e) {
            warning(sprintf("NMF failed to converge: %s", e$message))
            return(NULL)
          })
          if (is.null(nmf_res)) {
            return(list(
              genes = character(0),
              weights = numeric(0),
              scores = numeric(0),
              performance = list(error = "NMF convergence failure"),
              model = NULL
            ))
          }
          
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

#' @export
train_meta_learner_v1_test <- function(l1_signatures,
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
train_meta_learner_v2_test <- function(l1_signatures,
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
train_meta_learner_v3_test <- function(l1_signatures,
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

train_meta_learner_v4_test <- function(l1_signatures,
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
train_meta_learner_v5_test <- function(...) {
  .Deprecated("TML7", package = "myR")
  TML7(...)
}

#' @export
train_meta_learner_v6_test <- function(...) {
  .Deprecated("TML7", package = "myR")
  TML7(...)
}