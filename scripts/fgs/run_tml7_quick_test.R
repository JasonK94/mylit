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

# 명령줄 인자 파싱
# Usage: Rscript run_tml7_quick_test.R [method1] [method2] ... [target_var] [data_path] [cv_group_var]
# Example: Rscript run_tml7_quick_test.R glm ranger xgbTree
# Example: Rscript run_tml7_quick_test.R glm ranger g3 /data/user3/sobj/is5.qs emrid

args <- commandArgs(trailingOnly = TRUE)

# 기본값
l2_methods <- c('glm', 'ranger', 'xgbTree')
TARGET_VAR <- "response"
DATA_PATH <- "/data/user3/sobj/data_seurat_251104.qs"
CV_GROUP_VAR <- "hos_no"

# 인자 파싱: methods는 알려진 method 이름들, 나머지는 설정
if (length(args) > 0) {
  supported_methods <- c("glm", "ranger", "xgbTree", "glmnet", "svmRadial", 
                         "mlp", "mlpKerasDropout", "nnet", "earth")
  
  # methods와 설정 구분
  method_indices <- which(args %in% supported_methods)
  if (length(method_indices) > 0) {
    l2_methods <- args[method_indices]
    args <- args[-method_indices]
  }
  
  # 남은 인자 처리 (target_var, data_path, cv_group_var 순서)
  if (length(args) >= 1) TARGET_VAR <- args[1]
  if (length(args) >= 2) DATA_PATH <- args[2]
  if (length(args) >= 3) CV_GROUP_VAR <- args[3]
}

# 데이터 로드
message("Loading data...")
fgs2 <- qs::qread("/data/user3/sobj/fgs/fgs2.qs")
message(sprintf("Loading Seurat data from: %s", DATA_PATH))
data_seurat <- qs::qread(DATA_PATH)
message("Data loaded.\n")

# target_var 확인
if (!TARGET_VAR %in% colnames(data_seurat@meta.data)) {
  if ("response" %in% colnames(data_seurat@meta.data)) {
    message(sprintf("Warning: '%s' not found, using 'response' instead", TARGET_VAR))
    TARGET_VAR <- "response"
  } else if ("g3" %in% colnames(data_seurat@meta.data)) {
    message(sprintf("Warning: '%s' not found, using 'g3' instead", TARGET_VAR))
    TARGET_VAR <- "g3"
  } else {
    stop(sprintf("Target variable '%s' not found in meta.data.", TARGET_VAR))
  }
}

message("========================================")
message("TML7 Quick Test")
message("========================================")
message(sprintf("Methods: %s", paste(l2_methods, collapse = ", ")))
message(sprintf("target_var: %s", TARGET_VAR))
message(sprintf("cv_group_var: %s", CV_GROUP_VAR))
message(sprintf("data_path: %s", DATA_PATH))
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
  target_var = TARGET_VAR,
  l2_methods = l2_methods,
  cv_group_var = CV_GROUP_VAR,
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

