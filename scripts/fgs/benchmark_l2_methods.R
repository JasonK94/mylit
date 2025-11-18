#!/usr/bin/env Rscript
# Benchmark individual L2 methods for TML7
# Measures execution time for each method separately

# Use FGS-specific lightweight environment initialization
# This avoids loading all packages from start.R and prevents parallel worker issues
suppressPackageStartupMessages({
  source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")
})

library(qs)
library(pROC)
suppressPackageStartupMessages(devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE))

# Output directory
fgs_root <- '/home/user3/data_user3/git_repo/_wt/fgs'
output_dir <- file.path(fgs_root, 'outputs', 'benchmark_l2')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load pre-computed FGS results
fgs_path <- '/data/user3/sobj/fgs/fgs2.qs'
message("Loading pre-computed FGS results from: ", fgs_path)
fgs_result <- qs::qread(fgs_path)

# Load Seurat data
data_path <- '/data/user3/sobj/data_seurat_251104.qs'
message("Loading Seurat object from: ", data_path)
data_seurat <- qs::qread(data_path)

# L2 methods to benchmark (one at a time)
l2_methods_list <- list(
  glm = "glm",
  glmnet = "glmnet",
  ranger = "ranger",
  xgbTree = "xgbTree",
  svmRadial = "svmRadial",
  mlp = "mlp",
  mlpKerasDropout = "mlpKerasDropout",
  nnet = "nnet",
  earth = "earth"
)

# Results storage
benchmark_results <- data.frame(
  method = character(),
  elapsed_time_sec = numeric(),
  success = logical(),
  error_message = character(),
  auc = numeric(),
  stringsAsFactors = FALSE
)

message("\n=== Starting L2 Method Benchmark ===")
message("Total methods to test: ", length(l2_methods_list))
message("Using cv_group_var: emrid")
message("k_folds: 5\n")

for (i in seq_along(l2_methods_list)) {
  method_name <- names(l2_methods_list)[i]
  method_value <- l2_methods_list[[i]]
  
  message(sprintf("[%d/%d] Testing method: %s", i, length(l2_methods_list), method_name))
  message("----------------------------------------")
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    tml <- TML7(
      l1_signatures = fgs_result,
      holdout_data = data_seurat,
      target_var = "response",
      l2_methods = method_value,
      metric = "AUC",
      cv_group_var = "emrid",
      fgs_seed = 42,
      allow_parallel = FALSE,
      k_folds = 5
    )
    
    # Calculate AUC on training data
    sig_cols <- setdiff(colnames(tml$l2_train), ".target")
    score_mat <- tml$l2_train[, sig_cols, drop = FALSE]
    prob_mat <- predict(tml$best_model, newdata = score_mat, type = "prob")
    pos_class <- tml$positive_class
    prob_vec <- prob_mat[, pos_class]
    prob_vec <- pmin(pmax(prob_vec, 1e-6), 1 - 1e-6)
    roc_obj <- pROC::roc(response = tml$l2_train$.target, predictor = prob_vec, quiet = TRUE)
    auc_val <- as.numeric(roc_obj$auc)
    
    list(
      success = TRUE,
      auc = auc_val,
      best_model = tml$best_model_name,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      success = FALSE,
      auc = NA_real_,
      best_model = NA_character_,
      error = conditionMessage(e)
    )
  })
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  benchmark_results <- rbind(benchmark_results, data.frame(
    method = method_name,
    elapsed_time_sec = elapsed,
    success = result$success,
    error_message = ifelse(is.na(result$error), "", result$error),
    auc = ifelse(is.na(result$auc), NA_real_, result$auc),
    stringsAsFactors = FALSE
  ))
  
  if (result$success) {
    message(sprintf("✓ Success: %.2f seconds, AUC=%.4f, model=%s", 
                    elapsed, result$auc, result$best_model))
  } else {
    message(sprintf("✗ Failed: %.2f seconds, error: %s", elapsed, result$error))
  }
  message("")
}

# Save results
timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
results_file <- file.path(output_dir, paste0('l2_benchmark_', timestamp, '.csv'))
utils::write.csv(benchmark_results, results_file, row.names = FALSE)

# Summary
message("\n=== Benchmark Summary ===")
message("Results saved to: ", results_file)
message("\nMethod Performance:")
print(benchmark_results[, c("method", "elapsed_time_sec", "success", "auc")])

successful <- sum(benchmark_results$success)
message(sprintf("\nSuccessful: %d/%d", successful, nrow(benchmark_results)))
message(sprintf("Total time: %.2f seconds (%.2f minutes)", 
                sum(benchmark_results$elapsed_time_sec),
                sum(benchmark_results$elapsed_time_sec) / 60))

