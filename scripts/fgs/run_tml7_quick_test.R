#!/usr/bin/env Rscript
# Quick TML7 test script with taskset support
# Usage: taskset -c 0-7 Rscript run_tml7_quick_test.R [method1] [method2] ...
# Example: taskset -c 0-7 Rscript run_tml7_quick_test.R glm ranger xgbTree

# FGS 환경 초기화 (start.R 사용, 병렬 처리 비활성화)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')

# 필요한 패키지 로드
library(qs)
library(pROC)

# FGS 함수 로드
devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE)

# 데이터 로드
message("Loading data...")
fgs2 <- qs::qread("/data/user3/sobj/fgs/fgs2.qs")
data_seurat <- qs::qread("/data/user3/sobj/data_seurat_251104.qs")
message("Data loaded.\n")

# 명령줄 인자에서 methods 읽기 (없으면 기본값 사용)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  l2_methods <- args
} else {
  # 기본 methods
  l2_methods <- c('glm', 'ranger', 'xgbTree')
}

message("========================================")
message("TML7 Quick Test")
message("========================================")
message(sprintf("Methods: %s", paste(l2_methods, collapse = ", ")))
message("target_var: response")
message("cv_group_var: hos_no")
message("k_folds: 5")
message("metric: AUC")
message("========================================\n")

# CPU 설정 (taskset 사용 시 이 설정들은 무시될 수 있음)
Sys.setenv(
  FGS_MAX_CPU_CORES = "1",
  FGS_BLAS_THREADS = "1",
  FGS_DISABLE_PARALLEL = "FALSE"
)

# TML7 실행
start_time <- Sys.time()

tml_result <- TML7(
  l1_signatures = fgs2,
  holdout_data = data_seurat,
  target_var = 'response',
  l2_methods = l2_methods,
  cv_group_var = 'hos_no',
  fgs_seed = 42,
  metric = "AUC",
  k_folds = 5
)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

message("\n========================================")
message("Results")
message("========================================")
message(sprintf("Elapsed time: %.2f seconds (%.2f minutes)", elapsed, elapsed / 60))
message(sprintf("Best model: %s", tml_result$best_model_name))

# 성능 계산
sig_cols <- setdiff(colnames(tml_result$l2_train), ".target")
score_mat <- tml_result$l2_train[, sig_cols, drop = FALSE]
prob_mat <- predict(tml_result$best_model, newdata = score_mat, type = "prob")
pos_class <- tml_result$positive_class
prob_vec <- prob_mat[, pos_class]
prob_vec <- pmin(pmax(prob_vec, 1e-6), 1 - 1e-6)
roc_obj <- pROC::roc(response = tml_result$l2_train$.target, predictor = prob_vec, quiet = TRUE)
auc_val <- as.numeric(roc_obj$auc)

message(sprintf("Training AUC: %.4f", auc_val))
if (!is.null(tml_result$best_model$results)) {
  cv_roc <- max(tml_result$best_model$results$ROC, na.rm = TRUE)
  message(sprintf("CV ROC (best): %.4f", cv_roc))
}

# 결과 저장
output_file <- "/data/user3/sobj/tml2_v7_quick_test.qs"
message(sprintf("\nSaving result to: %s", output_file))
qs::qsave(tml_result, output_file)
message("✓ Saved.")

message("\n========================================\n")

