#!/usr/bin/env Rscript
# Test individual L2 methods for TML7
# Usage: 
#   Rscript test_l2_method_individual.R <method_name>
#   or set METHOD environment variable
#   or edit METHOD variable in this script
#
# Example:
#   taskset -c 0-7 Rscript test_l2_method_individual.R glm
#   taskset -c 0-7 Rscript test_l2_method_individual.R xgbTree

# FGS 환경 초기화 (start.R 사용, 병렬 처리 비활성화)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')

# 필요한 패키지 로드
library(qs)
library(pROC)

# FGS 함수 로드
devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE)

# ===== 설정 =====
# 명령줄 인자 파싱
# Usage: Rscript test_l2_method_individual.R <method> [target_var] [data_path] [cv_group_var]
# Example: Rscript test_l2_method_individual.R glm response /data/user3/sobj/data_seurat_251104.qs hos_no
# Example: Rscript test_l2_method_individual.R glm g3 /data/user3/sobj/is5.qs emrid

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # 환경 변수나 스크립트 내 변수 확인
  if (exists("METHOD") && !is.null(METHOD)) {
    METHOD <- METHOD
  } else {
    METHOD <- Sys.getenv("METHOD", unset = NA)
    if (is.na(METHOD)) {
      stop("Method name required. Usage: Rscript test_l2_method_individual.R <method> [target_var] [data_path] [cv_group_var]\n",
           "  Example: Rscript test_l2_method_individual.R glm response\n",
           "  Example: Rscript test_l2_method_individual.R glm g3 /data/user3/sobj/is5.qs emrid\n",
           "\nSupported methods: glm, ranger, xgbTree, glmnet, svmRadial, mlp, mlpKerasDropout, nnet, earth")
    }
  }
  TARGET_VAR <- Sys.getenv("TARGET_VAR", unset = "response")
  DATA_PATH <- Sys.getenv("DATA_PATH", unset = "/data/user3/sobj/data_seurat_251104.qs")
  CV_GROUP_VAR <- Sys.getenv("CV_GROUP_VAR", unset = "hos_no")
} else {
  METHOD <- args[1]
  TARGET_VAR <- if (length(args) >= 2) args[2] else "response"
  DATA_PATH <- if (length(args) >= 3) args[3] else "/data/user3/sobj/data_seurat_251104.qs"
  CV_GROUP_VAR <- if (length(args) >= 4) args[4] else "hos_no"
}

# 지원되는 methods 확인
supported_methods <- c("glm", "ranger", "xgbTree", "glmnet", "svmRadial", 
                       "mlp", "mlpKerasDropout", "nnet", "earth")

if (!METHOD %in% supported_methods) {
  stop(sprintf("Unsupported method: %s\nSupported methods: %s", 
               METHOD, paste(supported_methods, collapse = ", ")))
}

message("========================================")
message(sprintf("Testing L2 method: %s", METHOD))
message("========================================\n")

# CPU 설정 (taskset 사용 시 이 설정들은 무시될 수 있음)
Sys.setenv(
  FGS_MAX_CPU_CORES = "1",
  FGS_BLAS_THREADS = "1",
  FGS_DISABLE_PARALLEL = "FALSE"
)

# 데이터 로드
message("Loading data...")
fgs2 <- qs::qread("/data/user3/sobj/fgs/fgs2.qs")
message(sprintf("Loading Seurat data from: %s", DATA_PATH))
data_seurat <- qs::qread(DATA_PATH)
message("Data loaded.\n")

# target_var가 데이터에 존재하는지 확인
if (!TARGET_VAR %in% colnames(data_seurat@meta.data)) {
  # 자동 감지 시도: response 또는 g3
  if ("response" %in% colnames(data_seurat@meta.data)) {
    message(sprintf("Warning: '%s' not found, using 'response' instead", TARGET_VAR))
    TARGET_VAR <- "response"
  } else if ("g3" %in% colnames(data_seurat@meta.data)) {
    message(sprintf("Warning: '%s' not found, using 'g3' instead", TARGET_VAR))
    TARGET_VAR <- "g3"
  } else {
    stop(sprintf("Target variable '%s' not found in meta.data. Available columns: %s",
                 TARGET_VAR, paste(colnames(data_seurat@meta.data)[1:min(10, ncol(data_seurat@meta.data))], collapse = ", ")))
  }
}

# TML7 실행
message(sprintf("Running TML7 with method: %s", METHOD))
message(sprintf("  - target_var: %s", TARGET_VAR))
message(sprintf("  - cv_group_var: %s", CV_GROUP_VAR))
message(sprintf("  - data_path: %s", DATA_PATH))
message("  - k_folds: 5")
message("  - metric: AUC\n")

start_time <- Sys.time()

result <- tryCatch({
  tml_result <- TML7(
    l1_signatures = fgs2,
    holdout_data = data_seurat,
    target_var = TARGET_VAR,
    l2_methods = METHOD,  # 단일 method만 테스트
    cv_group_var = CV_GROUP_VAR,
    fgs_seed = 42,
    metric = "AUC",
    k_folds = 5
  )
  
  # 성능 계산
  sig_cols <- setdiff(colnames(tml_result$l2_train), ".target")
  score_mat <- tml_result$l2_train[, sig_cols, drop = FALSE]
  prob_mat <- predict(tml_result$best_model, newdata = score_mat, type = "prob")
  pos_class <- tml_result$positive_class
  prob_vec <- prob_mat[, pos_class]
  prob_vec <- pmin(pmax(prob_vec, 1e-6), 1 - 1e-6)
  roc_obj <- pROC::roc(response = tml_result$l2_train$.target, predictor = prob_vec, quiet = TRUE)
  auc_val <- as.numeric(roc_obj$auc)
  
  list(
    success = TRUE,
    tml = tml_result,
    auc = auc_val,
    best_model = tml_result$best_model_name,
    cv_metric = if (!is.null(tml_result$best_model$results)) {
      max(tml_result$best_model$results$ROC, na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}, error = function(e) {
  list(
    success = FALSE,
    error = conditionMessage(e),
    auc = NA_real_,
    best_model = NA_character_,
    cv_metric = NA_real_
  )
})

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

# 결과 출력
message("\n========================================")
message("Results")
message("========================================")

if (result$success) {
  message(sprintf("✓ SUCCESS: %s", METHOD))
  message(sprintf("  Elapsed time: %.2f seconds (%.2f minutes)", elapsed, elapsed / 60))
  message(sprintf("  Best model: %s", result$best_model))
  message(sprintf("  CV ROC (best): %.4f", result$cv_metric))
  message(sprintf("  Training AUC: %.4f", result$auc))
  
  # 결과 저장
  output_file <- sprintf("/data/user3/sobj/tml2_v7_%s.qs", METHOD)
  message(sprintf("\nSaving result to: %s", output_file))
  qs::qsave(result$tml, output_file)
  message("✓ Saved.")
  
} else {
  message(sprintf("✗ FAILED: %s", METHOD))
  message(sprintf("  Elapsed time: %.2f seconds", elapsed))
  message(sprintf("  Error: %s", result$error))
  
  # 패키지 의존성 문제인 경우 안내
  if (grepl("package required|package not available|No usable models", result$error, ignore.case = TRUE)) {
    method_deps <- list(
      glm = "base R (no extra package)",
      ranger = "install.packages('ranger')",
      xgbTree = "install.packages('xgboost')",
      glmnet = "install.packages('glmnet')",
      svmRadial = "install.packages('kernlab')",
      mlp = "install.packages('RSNNS')",
      mlpKerasDropout = "install.packages('keras')",
      nnet = "install.packages('nnet')",
      earth = "install.packages('earth')"
    )
    
    if (METHOD %in% names(method_deps)) {
      message(sprintf("\n  To install required package, run:\n    %s", method_deps[[METHOD]]))
    }
  }
}

message("\n========================================\n")

