#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  source("/home/user3/GJC_KDW_250721/start.R")
})

library(qs)
library(pROC)
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")

is5s_path <- "/data/user3/sobj/IS6_sex_added_0.1x_251110.qs"
output_dir <- "/home/user3/data_user3/git_repo/_wt/fgs/outputs/fgs_is5s"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading Seurat object: ", is5s_path)
is5s <- qs::qread(is5s_path)
if (!"response" %in% colnames(is5s@meta.data)) {
  if (!"g3" %in% colnames(is5s@meta.data)) {
    stop("Metadata lacks both 'response' and 'g3' columns.")
  }
  response_factor <- ifelse(
    is.na(is5s@meta.data$g3),
    NA_character_,
    ifelse(is5s@meta.data$g3 == 2, "R", "NR")
  )
  is5s@meta.data$response <- factor(response_factor, levels = c("NR", "R"))
}

methods_to_run <- c(
  "random_forest", "random_forest_ranger", "lasso", "ridge", "elastic_net",
  "pca_loadings", "nmf_loadings", "gam", "limma", "wilcoxon", "xgboost"
)

message("Running find_gene_signature_v5.4 across methods: ", paste(methods_to_run, collapse = ", "))
fgs <- find_gene_signature_v5.4(
  data = is5s,
  target_var = "response",
  control_vars = "hos_no",
  n_features = 200,
  method = methods_to_run,
  fgs_seed = 42
)

message("Training TML7 meta-learner...")
tml <- TML7(
  l1_signatures = fgs,
  holdout_data = is5s,
  target_var = "response",
  metric = "AUC",
  cv_group_var = "hos_no",
  fgs_seed = 42
)

message("Computing meta gene importance and meta signature scores...")
cmgi <- compute_meta_gene_importance(tml)
is5s_with_meta <- add_meta_signature_score(
  seurat_obj = is5s,
  gene_importance_result = cmgi,
  signature_name = "meta_signature_score"
)

positive_class <- cmgi$positive_class
pred_probs <- predict(tml$best_model, tml$l2_train, type = "prob")
prob_vec <- pred_probs[, positive_class]
prob_vec <- pmin(pmax(prob_vec, 1e-6), 1 - 1e-6)
logit_vec <- stats::qlogis(prob_vec)

meta_scores <- is5s_with_meta@meta.data[rownames(tml$l2_train), "meta_signature_score"]

score_df <- data.frame(
  cell = rownames(tml$l2_train),
  target = tml$l2_train$.target,
  prob = prob_vec,
  logit = logit_vec,
  meta_signature_score = meta_scores,
  stringsAsFactors = FALSE
)

auc <- as.numeric(pROC::auc(pROC::roc(score_df$target, score_df$prob, quiet = TRUE)))
correlation <- suppressWarnings(stats::cor(score_df$logit, score_df$meta_signature_score, use = "complete.obs"))

fgs_path <- file.path(output_dir, "fgs_is5s.qs")
tml_path <- file.path(output_dir, "tml_is5s.qs")
scores_path <- file.path(output_dir, "l2_predictions_with_meta.csv")
log_path <- file.path(output_dir, "run_log.txt")

qs::qsave(fgs, fgs_path)
qs::qsave(tml, tml_path)
utils::write.csv(score_df, scores_path, row.names = FALSE)

log_lines <- c(
  sprintf("fgs_methods=%s", paste(names(fgs), collapse = ",")),
  sprintf("positive_class=%s", positive_class),
  sprintf("AUC=%.4f", auc),
  sprintf("logit_vs_meta_cor=%.4f", correlation),
  sprintf("fgs_saved=%s", fgs_path),
  sprintf("tml_saved=%s", tml_path),
  sprintf("scores_saved=%s", scores_path),
  sprintf("timestamp=%s", Sys.time())
)

writeLines(log_lines, log_path)
cat(paste(log_lines, collapse = "\n"), "\n")

