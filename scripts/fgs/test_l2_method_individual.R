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
# 테스트할 method를 여기서 지정하거나, 명령줄 인자로 전달
# METHOD <- "glm"  # 예시: 이 줄의 주석을 해제하고 method 이름을 지정

# 명령줄 인자에서 method 읽기
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  METHOD <- args[1]
} else if (exists("METHOD") && !is.null(METHOD)) {
  # 스크립트 내에서 METHOD 변수가 설정되어 있으면 사용
} else {
  # 환경 변수에서 읽기
  METHOD <- Sys.getenv("METHOD", unset = NA)
  if (is.na(METHOD)) {
    stop("Method name required. Usage: Rscript test_l2_method_individual.R <method_name>\n",
         "  or set METHOD environment variable\n",
         "  or set METHOD variable in script\n",
         "\nSupported methods: glm, ranger, xgbTree, glmnet, svmRadial, mlp, mlpKerasDropout, nnet, earth")
  }
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
data_seurat <- qs::qread("/data/user3/sobj/data_seurat_251104.qs")
message("Data loaded.\n")

# TML7 실행
message(sprintf("Running TML7 with method: %s", METHOD))
message("  - target_var: response")
message("  - cv_group_var: hos_no")
message("  - k_folds: 5")
message("  - metric: AUC\n")

start_time <- Sys.time()

result <- tryCatch({
  tml_result <- TML7(
    l1_signatures = fgs2,
    holdout_data = data_seurat,
    target_var = 'response',
    l2_methods = METHOD,  # 단일 method만 테스트
    cv_group_var = 'hos_no',
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
}

message("\n========================================\n")

