#!/usr/bin/env Rscript
# Full TML7 pipeline on IS6 dataset (is5s)
# Uses expanded L2 methods with CPU core limiting

# Use FGS-specific lightweight environment initialization
# This avoids loading all packages from start.R and prevents parallel worker issues
suppressPackageStartupMessages({
  source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")
})

library(qs)
library(pROC)
suppressPackageStartupMessages(devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE))

# Configuration
fgs_root <- '/home/user3/data_user3/git_repo/_wt/fgs'
output_dir <- file.path(fgs_root, 'outputs', 'fgs_is5s')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj_path <- '/data/user3/sobj/IS6_sex_added_0.1x_251110.qs'
message('Loading is5s from ', obj_path)
is5s <- qs::qread(obj_path)

# Derive response column if needed
if (!'response' %in% colnames(is5s@meta.data)) {
  if (!'g3' %in% colnames(is5s@meta.data)) {
    stop("Neither 'response' nor 'g3' columns exist in meta.data")
  }
  response_factor <- ifelse(
    is.na(is5s@meta.data$g3),
    NA_character_,
    ifelse(is5s@meta.data$g3 == 2, 'R', 'NR')
  )
  is5s@meta.data$response <- factor(response_factor, levels = c('NR', 'R'))
  message('Derived response column from g3 (NR=1, R=2).')
}

# Step 1: Find gene signatures
methods_to_run <- c(
  'random_forest', 'random_forest_ranger', 'lasso', 'ridge', 'elastic_net',
  'pca_loadings', 'nmf_loadings', 'gam', 'limma', 'wilcoxon', 'xgboost'
)

message('Running find_gene_signature_v5.4 (n_features=200, control_vars=hos_no) ...')
fgs_result <- find_gene_signature_v5.4(
  data = is5s,
  target_var = 'response',
  control_vars = 'hos_no',
  n_features = 200,
  method = methods_to_run,
  fgs_seed = 42
)
message('Collected ', length(fgs_result), ' signatures: ', paste(names(fgs_result), collapse = ', '))

# Step 2: Train TML7 with expanded methods
# Start with fast methods only for initial testing
l2_methods <- c('glm', 'glmnet', 'ranger', 'xgbTree', 'nnet')
message('Training TML7 with group_var=hos_no and methods: ', paste(l2_methods, collapse = ', '))

tml <- TML7(
  l1_signatures = fgs_result,
  holdout_data = is5s,
  target_var = 'response',
  l2_methods = l2_methods,
  metric = 'AUC',
  cv_group_var = 'hos_no',
  fgs_seed = 42,
  allow_parallel = FALSE,
  k_folds = 5
)

message('Best model: ', tml$best_model_name)

# Step 3: Evaluate performance
sig_cols <- setdiff(colnames(tml$l2_train), '.target')
score_mat <- tml$l2_train[, sig_cols, drop = FALSE]
prob_mat <- predict(tml$best_model, newdata = score_mat, type = 'prob')
pos_class <- tml$positive_class
prob_vec <- prob_mat[, pos_class]
prob_vec <- pmin(pmax(prob_vec, 1e-6), 1 - 1e-6)
roc_obj <- pROC::roc(response = tml$l2_train$.target, predictor = prob_vec, quiet = TRUE)
auc_val <- as.numeric(roc_obj$auc)
message(sprintf('Observed ROC AUC: %.4f', auc_val))

# Step 4: Compute meta gene importance
cmgi <- compute_meta_gene_importance(tml)

# Step 5: Save outputs
stamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
fgs_path <- file.path(output_dir, paste0('fgs_is5s_v54_n200_', stamp, '.qs'))
tml_path <- file.path(output_dir, paste0('tml7_is5s_groupcv_', stamp, '.qs'))
cmgi_path <- file.path(output_dir, paste0('cmgi_is5s_groupcv_', stamp, '.qs'))

qs::qsave(fgs_result, fgs_path)
qs::qsave(tml, tml_path)
qs::qsave(cmgi, cmgi_path)

# Step 6: Write log
log_path <- file.path(output_dir, paste0('tml7_run_', stamp, '.log'))
log_lines <- c(
  sprintf('timestamp=%s', Sys.time()),
  sprintf('fgs_file=%s', fgs_path),
  sprintf('tml_file=%s', tml_path),
  sprintf('cmgi_file=%s', cmgi_path),
  sprintf('best_model=%s', tml$best_model_name),
  sprintf('positive_class=%s', tml$positive_class),
  sprintf('auc=%.4f', auc_val),
  sprintf('l2_methods=%s', paste(l2_methods, collapse = ',')),
  sprintf('cv_group_var=%s', 'hos_no'),
  sprintf('n_l1_signatures=%d', length(fgs_result))
)
writeLines(log_lines, log_path)
cat(paste(log_lines, collapse = '\n'), '\n')

message('\n=== Pipeline Complete ===')
message('Outputs saved to: ', output_dir)

