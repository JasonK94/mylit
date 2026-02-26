#!/usr/bin/env Rscript
# scripts/consensus/run_deg_consensus_cli.R
#
# Flexible CLI script for DEG Consensus Analysis using optparse
#
# Usage Example:
# Rscript run_deg_consensus_cli.R \
#   --input /data/user3/sobj/is2_IS_3_clustered.qs \
#   --output /data/user3/sobj/consensus/MyAnalysis \
#   --cluster anno3big \
#   --group g3 \
#   --contrast "2 - 1" \
#   --covariates "sex,age" \
#   --methods "all" \
#   --cores 8

renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

suppressPackageStartupMessages({
    library(optparse)
    library(qs)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(future)
    library(doFuture)
})

# --- Option Parser ---
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Path to input Seurat object (.qs or .rds) [Required]", metavar = "FILE"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "output",
        help = "Output directory [default= %default]", metavar = "DIR"
    ),
    make_option(c("-c", "--cluster"),
        type = "character", default = "seurat_clusters",
        help = "Cluster column name [default= %default]"
    ),
    make_option(c("-s", "--sampleid"),
        type = "character", default = "hos_no",
        help = "Sample ID column name [default= %default]"
    ),
    make_option(c("-g", "--group"),
        type = "character", default = "g3",
        help = "Group/Condition column name [default= %default]"
    ),
    make_option(c("--contrast"),
        type = "character", default = NULL,
        help = "Contrast string (e.g., '2 - 1' or 'Disease - Control') [Required]"
    ),
    make_option(c("--covariates"),
        type = "character", default = "",
        help = "Comma-separated list of covariates (e.g., 'sex,age'). Empty for none."
    ),
    make_option(c("--methods"),
        type = "character", default = "all",
        help = "Methods to run (comma-separated, e.g., 'edgeR-LRT,muscat-edgeR') or 'all' [default= %default]"
    ),
    make_option(c("--cores"),
        type = "integer", default = 4,
        help = "Number of cores for parallel processing [default= %default]"
    ),
    make_option(c("--renv"),
        type = "character", default = NULL,
        help = "Path to renv library to use (optional)"
    ),
    make_option(c("--src"),
        type = "character", default = "/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus",
        help = "Directory containing R source files [default= %default]"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# --- Validation ---
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file argument (--input) is required.")
}
if (is.null(opt$contrast)) {
    print_help(opt_parser)
    stop("Contrast argument (--contrast) is required.")
}

# --- Setup Environment ---
if (!is.null(opt$renv) && dir.exists(opt$renv)) {
    .libPaths(c(opt$renv, .libPaths()))
    message("Added renv library: ", opt$renv)
}

# Source Functions
if (dir.exists(opt$src)) {
    src_files <- list.files(opt$src, pattern = "\\.R$", full.names = TRUE)
    if (length(src_files) > 0) {
        message("Sourcing ", length(src_files), " files from ", opt$src)
        for (f in src_files) source(f)
    } else {
        warning("No R files found in ", opt$src)
    }
} else {
    # Fallback to devtools if myR exists locally
    if (dir.exists("myR")) devtools::load_all("myR")
}

# Parallel Setup
registerDoFuture()
plan(multisession, workers = opt$cores)
message("Parallel backend registered with ", opt$cores, " cores.")

# --- Load Data ---
message("Loading data from: ", opt$input)
if (grepl("\\.qs$", opt$input)) {
    sobj <- qs::qread(opt$input)
} else {
    sobj <- readRDS(opt$input)
}

# --- Parse Arguments ---
covar_vec <- if (nchar(opt$covariates) > 0) unlist(strsplit(opt$covariates, ",")) else NULL
if (!is.null(covar_vec)) covar_vec <- trimws(covar_vec)

methods_vec <- if (opt$methods == "all") {
    c(
        "edgeR-LRT", "edgeR-QLF",
        "DESeq2-Wald", "DESeq2-LRT",
        "limma-voom", "limma-trend",
        "dream", "nebula",
        "muscat-edgeR", "muscat-DESeq2",
        "muscat-limma-voom", "muscat-limma-trend"
    )
} else {
    trimws(unlist(strsplit(opt$methods, ",")))
}

# Create Output Directory
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# --- Metadata Checks ---
meta <- sobj@meta.data
if (!opt$cluster %in% names(meta)) stop(paste("Cluster column", opt$cluster, "not found."))
if (!opt$group %in% names(meta)) stop(paste("Group column", opt$group, "not found."))
if (!is.null(covar_vec)) {
    missing_cov <- covar_vec[!covar_vec %in% names(meta)]
    if (length(missing_cov) > 0) stop(paste("Covariates not found:", paste(missing_cov, collapse = ", ")))
}

# --- Run Analysis ---
message("----------------------------------------------------------------")
message("Starting Analysis")
message("Cluster: ", opt$cluster)
message("Group: ", opt$group)
message("Contrast: ", opt$contrast)
message("Covariates: ", paste(covar_vec, collapse = ", "))
message("Methods: ", paste(methods_vec, collapse = ", "))
message("Output: ", opt$output)
message("----------------------------------------------------------------")

out_file <- file.path(opt$output, "results_consensus.qs")

# 1. Run DEG Pipeline
# Using run_deg_consensus (Low level) and then pipeline visualization
res <- tryCatch(
    {
        run_deg_consensus(
            sobj = sobj,
            contrast = opt$contrast,
            methods = methods_vec,
            cluster_id = opt$cluster,
            sample_id = opt$sampleid, # Hardcoded or add arg? Assuming standard.
            group_id = opt$group,
            batch_id = NULL,
            covar_effects = covar_vec,
            n_cores = opt$cores,
            verbose = TRUE
        )
    },
    error = function(e) {
        message("Error in run_deg_consensus: ", e$message)
        return(NULL)
    }
)

if (is.null(res)) stop("DEG Analysis failed.")

# 2. Standardize
message("Standardizing results...")
std <- lapply(res$methods_run, function(m) {
    standardize_deg_results(res$results[[m]], m)
})
names(std) <- res$methods_run
res$standardized_results <- std

# 3. Build Matrices
mat <- build_deg_matrices(std, fdr_threshold = 0.1, significance_mode = "pvalue", pvalue_threshold = 0.05)
res$deg_matrices <- mat

# 4. Consensus Scores
agree <- compute_agreement_scores(mat$significance)
cons_scores <- compute_consensus_scores(mat, agree)
res$consensus_scores <- cons_scores

# 5. Consensus List
cons_list <- generate_consensus_deg_list(cons_scores, fdr_threshold = 0.05)
res$consensus_deg_list <- cons_list

# Save Result
qs::qsave(res, out_file)
message("Results saved to: ", out_file)

# 6. Plots
plot_prefix <- file.path(opt$output, "AG_Consensus")

# Helper for caption
get_git_hash <- function() {
    tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")
}
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
git_hash <- get_git_hash()

caption_text <- paste0(
    "Cluster: ", opt$cluster, "\n",
    "Formula: ~ ", opt$group, " + ", paste(covar_vec, collapse = " + "), " [Varies]\n",
    "Script: run_deg_consensus_cli.R\n",
    "Time: ", timestamp, "\n",
    "Commit: ", git_hash
)

message("Generating plots...")
try({
    plot_consensus_meta_distributions(cons_scores, file_prefix = plot_prefix, caption = caption_text)
    plot_consensus_volcano(cons_scores, file_prefix = plot_prefix, fdr_threshold = 0.05, caption = caption_text)
    plot_consensus_heatmap(mat, cons_scores, file_prefix = plot_prefix, method_top_n = 20, caption = caption_text)

    pca_res <- perform_deg_pca(mat)
    if (!is.null(pca_res)) plot_method_pca(pca_res, file_prefix = plot_prefix, caption = caption_text)

    plot_gene_umap(mat, cons_scores, file_prefix = plot_prefix, caption = caption_text)
})

# 7. Method Similarity (Statistics)
message("Method Similarity Analysis...")
if (!is.null(mat$statistic)) {
    mat_stat <- mat$statistic
    valid_rows <- rowSums(!is.na(mat_stat)) >= 2
    if (sum(valid_rows) > 10) {
        mat_stat_sub <- mat_stat[valid_rows, ]
        cor_mat <- cor(mat_stat_sub, use = "pairwise.complete.obs", method = "spearman")

        p_file <- file.path(opt$output, "method_similarity_heatmap.png")
        png(p_file, width = 800, height = 800)
        pheatmap::pheatmap(cor_mat,
            main = paste("Method Similarity (Spearman) -", opt$cluster),
            display_numbers = TRUE,
            treeheight_row = 0,
            fontsize = 10
        )
        dev.off()
    }
}

# 8. Summary CSV
export_method_summary_csv(res, file.path(opt$output, "method_summary.csv"))

message("Analysis Complete.")
