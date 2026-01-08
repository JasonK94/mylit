#!/usr/bin/env Rscript

# ==============================================================================
# FGS & TML Pipeline CLI Script
#
# Usage:
#   Rscript scripts/fgs/run_fgs_pipeline.R --input <path> --output <path> ...
#
# Description:
#   Runs the full Find Gene Signature (FGS) and Transfer Meta Learning (TML) pipeline.
#   1. FGS: Identifies candidate signatures using multiple L1 methods.
#   2. TML: Ensembles signatures into a meta-learner (L2 model).
#   3. CMGI: Computes gene-level importance from the meta-learner.
#   4. AMSC: adds meta signature scores to the original Seurat object.
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(qs)
    library(Seurat)
    library(dplyr)
})

# --- Option Parsing ---
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Path to input Seurat object (.qs or .rds)", metavar = "file"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output path prefix (e.g. /path/to/results/my_analysis_). Files will be appended with suffixes.", metavar = "path"
    ),

    # FGS Parameters
    make_option(c("-t", "--target_var"),
        type = "character", default = NULL,
        help = "Target variable name in meta.data for classification (e.g. 'response', 'group')", metavar = "string"
    ),
    make_option(c("-c", "--control_vars"),
        type = "character", default = NULL,
        help = "Comma-separated list of control variables (e.g. 'batch,sex')", metavar = "string"
    ),
    make_option(c("-n", "--n_features"),
        type = "integer", default = 50,
        help = "Number of top features to select per method [default: %default]", metavar = "int"
    ),
    make_option(c("--fgs_methods"),
        type = "character", default = "all",
        help = "Comma-separated list of FGS methods or 'all' [default: %default]", metavar = "string"
    ),

    # TML Parameters
    make_option(c("--l2_methods"),
        type = "character", default = "glm,ranger,xgbTree",
        help = "Comma-separated list of L2 meta-learner methods [default: %default]", metavar = "string"
    ),
    make_option(c("--cv_folds"),
        type = "integer", default = 5,
        help = "Number of CV folds for TML [default: %default]", metavar = "int"
    ),
    make_option(c("--cv_group_var"),
        type = "character", default = NULL,
        help = "Variable for group-wise CV (often same as sample id) [default: NULL]", metavar = "string"
    ),

    # Execution Parameters
    make_option(c("--cores"),
        type = "integer", default = 16,
        help = "Number of CPU cores to use [default: %default]", metavar = "int"
    ),
    make_option(c("--name"),
        type = "character", default = "fgs_meta",
        help = "Name prefix for the meta signature score in Seurat object [default: %default]", metavar = "string"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$target_var)) {
    print_help(opt_parser)
    stop("Missing required arguments: --input, --output, and --target_var are mandatory.", call. = FALSE)
}

# --- Environment Setup ---
# Set environment variables for core control
Sys.setenv(
    FGS_MAX_CPU_CORES = as.character(opt$cores),
    FGS_BLAS_THREADS = "1", # Force sequential BLAS to avoid conflicts
    FGS_DISABLE_PARALLEL = "FALSE" # Allow parallel processing via future/parallel
)

# Source core functions
# Assuming the script is run from project root or scripts folder
# We try to locate init_fgs_env.R and signature.R
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))
if (length(script_dir) == 0) script_dir <- getwd() # Fallback for interactive

# Try to find init_fgs_env.R
init_env_path <- file.path(script_dir, "init_fgs_env.R")
if (!file.exists(init_env_path)) {
    # Check if we are in scripts/fgs
    if (file.exists("scripts/fgs/init_fgs_env.R")) {
        init_env_path <- "scripts/fgs/init_fgs_env.R"
    } else {
        warning("Could not find init_fgs_env.R. Assuming environment is already set.")
        init_env_path <- NULL
    }
}

if (!is.null(init_env_path)) {
    suppressPackageStartupMessages(source(init_env_path))
} else {
    # Fallback sourcing if init not found
    if (file.exists("myR/R/signature.R")) source("myR/R/signature.R")
    if (file.exists("myR/R/fgs_core.R")) source("myR/R/fgs_core.R")
    if (file.exists("myR/R/tml_core.R")) source("myR/R/tml_core.R")
    if (file.exists("myR/R/tml_utils.R")) source("myR/R/tml_utils.R")
}

# --- Logging Setup ---
log_dir <- "/data/user3/log/fgs"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, paste0("fgs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

# Print log location immediately
message("Logging to: ", log_file)

# Divert output to file (and keep console if interactive, but Rscript usually needs sinking)
# Note: R does not support split=TRUE for 'message' (stderr).
# So: stdout -> Console & File
#     stderr -> File only (Console will be silent for warnings/messages)
con <- file(log_file, open = "wt")
sink(con, type = "output", split = TRUE)
sink(con, type = "message")

# Log command line arguments
message("\n=== Execution Details ===")
message("Start Time: ", Sys.time())
message("Command Path: ", commandArgs(trailingOnly = FALSE)[1])
message("Script Path: ", if (exists("script_dir")) script_dir else "unknown")
message("Arguments:")
print(opt)
message("=========================\n")

# Ensure signature.R functions are loaded if init didn't load them (init loads start.R which might load local packages)
# But we should be explicit
if (file.exists("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/fgs_core.R")) {
    source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/fgs_core.R")
    source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/tml_core.R")
    source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/tml_utils.R")
    source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/utils_fgs.R") # Should be empty or utility only
}

# --- Main Pipeline ---

message("\n=== Starting FGS Pipeline ===")
message("Input: ", opt$input)
message("Output Prefix: ", opt$output)
message("Target: ", opt$target_var)

# 1. Load Data
message("\n--- Loading Data ---")
if (grepl("\\.qs$", opt$input)) {
    sobj <- qs::qread(opt$input)
} else {
    sobj <- readRDS(opt$input)
}
message("Dimensions: ", paste(dim(sobj), collapse = " x "))

# Prepare variables
control_vars <- if (!is.null(opt$control_vars)) trimws(strsplit(opt$control_vars, ",")[[1]]) else NULL
fgs_methods <- if (opt$fgs_methods == "all") {
    c(
        "random_forest_ranger", "lasso", "ridge", "elastic_net",
        "pca_loadings", "nmf_loadings", "gam", "limma", "wilcoxon", "xgboost"
    )
} else {
    trimws(strsplit(opt$fgs_methods, ",")[[1]])
}
l2_methods <- trimws(strsplit(opt$l2_methods, ",")[[1]])

# 2. Run FGS
message("\n--- 1. Running FGS (Find Gene Signature) ---")
message("Methods: ", paste(fgs_methods, collapse = ", "))
message("N Features: ", opt$n_features)

fgs_results <- find_gene_signature(
    data = sobj,
    target_var = opt$target_var,
    control_vars = control_vars,
    method = fgs_methods,
    n_features = opt$n_features,
    preprocess = TRUE
)

fgs_out_file <- paste0(opt$output, "_fgs.qs")
if (!dir.exists(dirname(fgs_out_file))) {
    dir.create(dirname(fgs_out_file), recursive = TRUE)
}
qs::qsave(fgs_results, fgs_out_file)
message("✓ FGS results saved to: ", fgs_out_file)

# 3. Run TML7
message("\n--- 2. Running TML7 (Transfer Meta Learning) ---")
message("L2 Methods: ", paste(l2_methods, collapse = ", "))
message("CV Folds: ", opt$cv_folds)
message("CV Group Var: ", opt$cv_group_var)

tml_results <- TML7(
    l1_signatures = fgs_results,
    holdout_data = sobj,
    target_var = opt$target_var,
    l2_methods = l2_methods,
    cv_folds = opt$cv_folds,
    cv_group_var = opt$cv_group_var,
    metric = "AUC" # Defaulting to AUC for now, hardcoded as reasonable default
)

tml_out_file <- paste0(opt$output, "_tml.qs")
qs::qsave(tml_results, tml_out_file)
message("✓ TML results saved to: ", tml_out_file)

# 4. Run CMGI (Compute Meta Gene Importance)
message("\n--- 3. Running CMGI (Compute Meta Gene Importance) ---")
cmgi_results <- tryCatch(
    {
        compute_meta_gene_importance(tml_results)
    },
    error = function(e) {
        warning("CMGI failed: ", e$message)
        return(NULL)
    }
)

if (!is.null(cmgi_results)) {
    cmgi_out_file <- paste0(opt$output, "_cmgi.qs")
    qs::qsave(cmgi_results, cmgi_out_file)

    # Also save CSV for convenience
    cmgi_csv_file <- paste0(opt$output, "_cmgi_genes.csv")
    write.csv(cmgi_results$gene_summary, cmgi_csv_file, row.names = FALSE)
    message("✓ CMGI results saved to: ", cmgi_out_file)
}

# 5. Run AMSC (Add Meta Signature Score)
message("\n--- 4. Running AMSC (Add Meta Signature Score) ---")
if (!is.null(cmgi_results)) {
    sobj_updated <- add_meta_signature_score(
        seurat_obj = sobj,
        gene_importance_result = cmgi_results,
        signature_name = paste0(opt$name, "_score")
    )

    sobj_out_file <- paste0(opt$output, "_sobj.qs")
    qs::qsave(sobj_updated, sobj_out_file)
    message("✓ Updated Seurat object saved to: ", sobj_out_file)
} else {
    warning("Skipping AMSC because CMGI results are missing.")
}

message("\n=== Pipeline Completed Successfully ===")
