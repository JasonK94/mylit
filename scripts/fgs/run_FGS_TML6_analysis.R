# ============================================================================
# FGS and TML6 Analysis Script
# Data: IS6_sex_added_251110.qs
# ============================================================================

# Load required packages
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs package required. Install with: install.packages('qs')")
}
library(qs)

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat package required. Install with: install.packages('Seurat')")
}
library(Seurat)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("dplyr package required. Install with: install.packages('dplyr')")
}
library(dplyr)

# Load the package
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools package required. Install with: install.packages('devtools')")
}
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")

# ============================================================================
# 1. Data Loading and Preparation
# ============================================================================

message("Loading data...")
is6 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# Check metadata
message("Checking metadata...")
colnames(is6@meta.data)
table(is6@meta.data$hos_no, useNA = "ifany") # sample ids
table(is6@meta.data$sex, useNA = "ifany") # biological cofactor
table(is6@meta.data$g3, useNA = "ifany") # target group variable; values are "NA", 1, 2
table(is6@meta.data$anno3.scvi, useNA = "ifany") # Cluster_annotated (scvi integrated)
table(is6@meta.data$GEM, useNA = "ifany") # batch variable

# Key variables
sample_key <- "hos_no"
batch_key <- "GEM"
group_key <- "g3"
covar_key <- "sex"
cluster_key <- "anno3.scvi"

# Check confounding relationships
meta_clinical <- is6@meta.data %>% distinct(hos_no, .keep_all = TRUE)
table(meta_clinical$hos_no, meta_clinical$GEM) # patients are in a specific GEM
table(meta_clinical$hos_no, meta_clinical$g3) # patients are in a specific group
table(meta_clinical$GEM, meta_clinical$g3) # some are perfectly separated

# ============================================================================
# 2. Data Cleaning: Remove NA from g3
# ============================================================================

message("Removing cells with NA in g3...")
is6_clean <- is6[, !is.na(is6@meta.data$g3)]
message(sprintf(
  "Original cells: %d, After removing NA: %d",
  ncol(is6), ncol(is6_clean)
))

# Convert g3 to factor to avoid numeric confusion
is6_clean@meta.data$g3 <- factor(is6_clean@meta.data$g3, levels = c("1", "2"))
message("g3 levels:", paste(levels(is6_clean@meta.data$g3), collapse = ", "))

# ============================================================================
# 3. FGS Analysis with Multiple Methods
# ============================================================================

message("\n=== Running FGS Analysis ===")

# Define control variables
# Note: FGS supports fixed effects only, so we use:
# - g3: main fixed effect
# - sex: covariate for sex differences
# - anno3.scvi: covariate for cluster-specific differences
# - GEM: batch effect (as fixed effect, since random effects not supported)

control_vars <- c("sex", "anno3.scvi", "GEM")

# Run FGS with multiple methods
# Using a subset of methods for faster execution
# Full list: c("random_forest", "random_forest_ranger", "lasso", "ridge",
#              "elastic_net", "pca_loadings", "nmf_loadings", "gam",
#              "limma", "wilcoxon", "xgboost")

fgs_methods <- c("lasso", "random_forest", "limma", "gam")

message("Running FGS with methods: ", paste(fgs_methods, collapse = ", "))

# Use FGS() which is an alias for find_gene_signature_v5.3
fgs_results <- FGS(
  data = is6_clean,
  target_var = "g3",
  control_vars = control_vars,
  method = fgs_methods,
  n_features = 100, # Number of top genes per method
  preprocess = TRUE,
  min_cells = 10,
  min_pct = 0.01,
  fgs_seed = 42,
  gam.min_unique = 15,
  gam.k = NULL, # Use dynamic k adjustment (v5.3 feature)
  gam.k_dynamic_factor = 5
)

# Save FGS results
message("\nSaving FGS results...")
qs::qsave(fgs_results, "/data/user3/sobj/FGS_results_IS6_251110.qs")
message("FGS results saved to: /data/user3/sobj/FGS_results_IS6_251110.qs")

# ============================================================================
# 4. Extract Signatures for TML6
# ============================================================================

message("\n=== Preparing Signatures for TML6 ===")

# Debug: Check FGS results structure
message("FGS results structure:")
message(sprintf("  Type: %s", class(fgs_results)))
message(sprintf("  Length: %d", length(fgs_results)))
message(sprintf("  Names: %s", paste(names(fgs_results), collapse = ", ")))

# Check first result structure
if (length(fgs_results) > 0) {
  first_result <- fgs_results[[1]]
  message(sprintf("  First result type: %s", class(first_result)))
  message(sprintf("  First result names: %s", paste(names(first_result), collapse = ", ")))
  if (!is.null(first_result$error)) {
    message(sprintf("  First result has error: %s", first_result$error))
  }
}

# Extract signatures from FGS results
# Convert to format expected by TML6 (named numeric vectors)
l1_signatures <- list()

for (method_name in names(fgs_results)) {
  result_item <- fgs_results[[method_name]]

  # Check if result has error
  if (!is.null(result_item$error)) {
    warning(sprintf("Method %s failed with error: %s", method_name, result_item$error))
    next
  }

  # Check if result has genes and weights
  if (!is.null(result_item$genes) && !is.null(result_item$weights)) {
    # Get genes and weights
    genes <- result_item$genes
    weights <- result_item$weights

    # Validate
    if (length(genes) == 0 || length(weights) == 0) {
      warning(sprintf("Method %s: Empty genes or weights", method_name))
      next
    }

    if (length(genes) != length(weights)) {
      warning(sprintf(
        "Method %s: Mismatch between genes (%d) and weights (%d)",
        method_name, length(genes), length(weights)
      ))
      next
    }

    # Create named numeric vector
    sig_vec <- weights
    names(sig_vec) <- genes

    l1_signatures[[method_name]] <- sig_vec

    message(sprintf("  %s: %d genes", method_name, length(genes)))
  } else {
    warning(sprintf(
      "Method %s did not produce valid signature (genes: %s, weights: %s)",
      method_name,
      !is.null(result_item$genes),
      !is.null(result_item$weights)
    ))
    # Debug: print available names
    message(sprintf("    Available names: %s", paste(names(result_item), collapse = ", ")))
  }
}

if (length(l1_signatures) == 0) {
  stop("No valid signatures extracted from FGS results")
}

message(sprintf("Total signatures prepared: %d", length(l1_signatures)))

# ============================================================================
# 5. TML6 Meta-Learner Training
# ============================================================================

message("\n=== Running TML6 Meta-Learner ===")

# Use the same cleaned data as holdout
# Note: In a real scenario, you might want to use a separate holdout set
meta_model <- TML6(
  l1_signatures = l1_signatures,
  holdout_data = is6_clean,
  target_var = "g3",
  l2_methods = c("glm", "ranger", "xgbTree"),
  k_folds = 5,
  metric = "AUC", # For binary classification
  fgs_seed = 42,
  layer = "data",
  allow_parallel = FALSE # Set to TRUE if parallel execution is desired
)

# Save TML6 results
message("\nSaving TML6 results...")
qs::qsave(meta_model, "/data/user3/sobj/TML6_results_IS6_251110.qs")
message("TML6 results saved to: /data/user3/sobj/TML6_results_IS6_251110.qs")

# ============================================================================
# 6. Compute Gene-Level Importance
# ============================================================================

message("\n=== Computing Gene-Level Importance ===")

gene_importance <- compute_meta_gene_importance(meta_model, normalize = TRUE)

# Save gene importance
qs::qsave(gene_importance, "/data/user3/sobj/gene_importance_IS6_251110.qs")
message("Gene importance saved to: /data/user3/sobj/gene_importance_IS6_251110.qs")

# Print summary
message("\n=== Summary ===")
message(sprintf("Best L2 model: %s", meta_model$best_model_name))
message(sprintf(
  "Best metric (%s): %.4f",
  meta_model$best_metric_name,
  max(meta_model$best_model$results[[meta_model$best_metric_name]], na.rm = TRUE)
))
message(sprintf("Top 10 genes by importance:"))
print(head(gene_importance$gene_summary, 10))

message("\n=== Analysis Complete ===")
message("Output files:")
message("  1. /data/user3/sobj/FGS_results_IS6_251110.qs")
message("  2. /data/user3/sobj/TML6_results_IS6_251110.qs")
message("  3. /data/user3/sobj/gene_importance_IS6_251110.qs")
