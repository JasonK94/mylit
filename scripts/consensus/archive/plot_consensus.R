#!/usr/bin/env Rscript
# plot_consensus.R
#
# Objective: Regenerate consensus plots from saved .qs results.
# Usage: Rscript scripts/consensus/plot_consensus.R [result_file_path]
# If no argument, defaults to processing AG_run1 broad/fine results.

renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

suppressPackageStartupMessages({
    library(qs)
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
    if (dir.exists("myR")) devtools::load_all("myR")
})

# Helper to get git hash
get_git_hash <- function() {
    tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")
}

# Main Plotting Function
generate_plots_for_result <- function(res_path, out_dir = NULL) {
    if (!file.exists(res_path)) {
        warning("File not found: ", res_path)
        return()
    }

    message("Loading result: ", res_path)
    res <- qs::qread(res_path)

    if (is.null(out_dir)) {
        out_dir <- dirname(res_path)
    }
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Extract Data
    cons_scores <- res$consensus_scores
    mat <- res$deg_matrices
    pca_res <- res$pca_result # Might be NULL if not run previously
    if (is.null(pca_res) || is.null(pca_res$coordinates)) {
        # Try re-running PCA if matrix is available
        if (!is.null(mat)) pca_res <- perform_deg_pca(mat)
    }

    # Construct Caption
    # Try to infer metadata
    cluster_info <- "Unknown Cluster"
    if (!is.null(res$raw_results) && !is.null(res$raw_results$run_info)) {
        # Try to find cluster_id in run_info
        # run_info is a list of call args
        c_id <- res$raw_results$run_info$cluster_id
        if (!is.null(c_id)) cluster_info <- as.character(c_id)
    } else {
        # Infer from path
        if (grepl("Broad_anno3big", res_path)) cluster_info <- "Broad_anno3big (anno3big)"
        if (grepl("Fine_anno3", res_path)) cluster_info <- "Fine_anno3 (anno3)"
    }

    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    git_hash <- get_git_hash()

    caption_text <- paste0(
        "Cluster: ", cluster_info, "\n",
        "Formula: ~ g3 + sex + age + (1|hos_no) [Varies by method]\n",
        "Script: scripts/consensus/plot_consensus.R\n",
        "Time: ", timestamp, "\n",
        "Commit: ", git_hash
    )

    plot_prefix <- file.path(out_dir, "AG_run1") # Or standardized name

    message("Generating plots in: ", out_dir)

    # 1. Meta Distribution
    plot_consensus_meta_distributions(cons_scores, file_prefix = plot_prefix, caption = caption_text)

    # 2. Volcano
    plot_consensus_volcano(cons_scores, file_prefix = plot_prefix, fdr_threshold = 0.05, caption = caption_text)

    # 3. Heatmap
    plot_consensus_heatmap(mat, cons_scores, file_prefix = plot_prefix, method_top_n = 10, caption = caption_text)

    # 4. Method PCA
    if (!is.null(pca_res)) {
        plot_method_pca(pca_res, file_prefix = plot_prefix, caption = caption_text)
    }

    # 5. UMAP
    plot_gene_umap(mat, cons_scores, file_prefix = plot_prefix, caption = caption_text)

    # 6. Method Similarity (Statistics)
    if (!is.null(mat$statistic)) {
        mat_stat <- mat$statistic
        valid_rows <- rowSums(!is.na(mat_stat)) >= 2
        if (sum(valid_rows) > 10) {
            mat_stat_sub <- mat_stat[valid_rows, ]
            cor_mat <- cor(mat_stat_sub, use = "pairwise.complete.obs", method = "spearman")

            p_file <- file.path(out_dir, "AG_run1_method_similarity_heatmap.png")
            png(p_file, width = 800, height = 800)
            pheatmap::pheatmap(cor_mat,
                main = paste("Method Similarity (Spearman) -", cluster_info),
                display_numbers = TRUE,
                treeheight_row = 0,
                fontsize = 10
            )
            # Add caption?? pheatmap doesn't easy support caption.
            # We can use grid to add text, but for now standard pheatmap.
            dev.off()
        }
    }

    # 7. Method Summary CSV
    export_method_summary_csv(res, file.path(out_dir, "AG_run1_method_summary.csv"))

    message("Plots generated.")
}

# Manually source functions if package not installed
src_files <- list.files("myR/R/deg_consensus", pattern = "\\.R$", full.names = TRUE)
for (f in src_files) source(f)


# Command Line Args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
    for (f in args) {
        generate_plots_for_result(f)
    }
} else {
    # Default behavior for AG_run1
    base_dir <- "/data/user3/sobj/consensus/AG_run1"

    # Broad
    broad_file <- file.path(base_dir, "Broad_anno3big", "results_consensus.qs")
    if (file.exists(broad_file)) {
        generate_plots_for_result(broad_file)
    } else {
        message("Default Broad file not found: ", broad_file)
    }

    # Fine
    fine_file <- file.path(base_dir, "Fine_anno3", "results_consensus.qs")
    if (file.exists(fine_file)) {
        generate_plots_for_result(fine_file)
    } else {
        message("Default Fine file not found: ", fine_file)
    }
}
