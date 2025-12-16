#!/usr/bin/env Rscript
# MASC Plot Generator (from cached results)
#
# - masc_results.qs 파일만 읽어서 플롯 생성
# - 분석 단계(.masc_run_analysis)를 건너뛰고 빠르게 플롯만 생성
#
# Usage example:
#   Rscript scripts/masc/plot_masc.R \
#     --results /data/user3/sobj/masc/stroke_complex/masc_anno3_complex_results.qs \
#     -o /data/user3/sobj/masc/stroke_complex \
#     --prefix masc_anno3_complex \
#     --cluster_var anno3 \
#     --contrast_var g3

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--results"), type = "character", help = "Path to masc_results.qs file"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for plots"),
  make_option(c("--prefix"), type = "character", default = "masc", help = "Output prefix for plot files"),
  make_option(c("--cluster_var"), type = "character", help = "Cluster variable name (for labeling)"),
  make_option(c("--contrast_var"), type = "character", help = "Contrast variable name (for labeling)"),
  make_option(c("--masc_r_path"), type = "character",
              default = "/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/myR/R/masc.R",
              help = "Path to masc.R (defaults to current worktree)"),
  make_option(c("--renv"), type = "character",
              default = "/home/user3/GJC_KDW_250721/renv",
              help = "Optional renv directory path. If specified and exists, activates renv."),
  make_option(c("--start_r_path"), type = "character",
              default = "/home/user3/GJC_KDW_250721/start.R",
              help = "Optional start.R path (loads env/renv). If missing, ignored.")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$results) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("Missing required args: --results, --output_dir", call. = FALSE)
}

if (!file.exists(opt$results)) {
  stop("Results file not found: ", opt$results, call. = FALSE)
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Activate renv if specified and exists
if (!is.null(opt$renv) && opt$renv != "" && dir.exists(opt$renv)) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    cat("Warning: renv package not available, skipping renv activation\n")
  } else {
    tryCatch({
      renv::activate(opt$renv)
      cat(sprintf("Activated renv from: %s\n", opt$renv))
    }, error = function(e) {
      cat(sprintf("Warning: Failed to activate renv from %s: %s\n", opt$renv, e$message))
    })
  }
}

# Also try start.R if it exists (for backward compatibility)
if (file.exists(opt$start_r_path)) {
  suppressWarnings(try(source(opt$start_r_path), silent = TRUE))
}

if (!file.exists(opt$masc_r_path)) stop("masc.R not found: ", opt$masc_r_path, call. = FALSE)
source(opt$masc_r_path)

cat(sprintf("Loading MASC results from: %s\n", opt$results))
masc_results <- qs::qread(opt$results)

if (!inherits(masc_results, "data.frame")) {
  stop("Results file does not contain a data.frame. Found: ", class(masc_results)[1], call. = FALSE)
}

if (!"cluster" %in% colnames(masc_results)) {
  stop("Results data.frame does not have 'cluster' column.", call. = FALSE)
}

cat(sprintf("Loaded results for %d clusters\n", nrow(masc_results)))

# Generate plots using the plot function from masc.R
cluster_var <- if (!is.null(opt$cluster_var)) opt$cluster_var else "cluster"
contrast_var <- if (!is.null(opt$contrast_var)) opt$contrast_var else "contrast"

plots_path <- file.path(opt$output_dir, paste0(opt$prefix, "_plots.qs"))

# Extract fixed_effects and random_effects from column names if possible
# (Note: This is a heuristic approach; ideally these should be passed as parameters)
fixed_effects <- NULL
random_effects <- NULL
model_formula <- NULL

plots <- .masc_plot_bundle(
  masc_results = masc_results,
  cluster_var = cluster_var,
  contrast_var = contrast_var,
  save = TRUE,
  save_path = plots_path,
  verbose = TRUE,
  fixed_effects = fixed_effects,
  random_effects = random_effects,
  model_formula = model_formula
)

if (!is.null(plots)) {
  cat("\n=== Plots generated successfully ===\n")
  cat(sprintf("Plot files saved to: %s\n", opt$output_dir))
  if ("or_forest" %in% names(plots)) {
    cat("  - OR forest plot: ", paste0(opt$prefix, "_plots_or_forest.{png,pdf}\n"))
  }
  if ("pvalue_bar" %in% names(plots)) {
    cat("  - P-value bar plot: ", paste0(opt$prefix, "_plots_pvalue_bar.{png,pdf}\n"))
  }
} else {
  cat("\nWarning: No plots were generated.\n")
}

cat("\nDone.\n")

