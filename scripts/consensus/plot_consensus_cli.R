#!/usr/bin/env Rscript
# scripts/consensus/plot_consensus_cli.R
#
# CLI script for regenerating plots from DEG Consensus results (.qs)
#
# Usage Example:
# Rscript scripts/consensus/plot_consensus_cli.R \
#   --input /data/user3/sobj/consensus/AG_run1/Broad_anno3big/results_consensus.qs \
#   --output /data/user3/sobj/consensus/AG_run1/Broad_anno3big/plots \
#   --cluster "Broad_Cluster" \
#   --group "Condition" \
#   --covariates "sex,age"

renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

suppressPackageStartupMessages({
    library(optparse)
    library(qs)
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
    library(ggridges)
})

# --- Option Parser ---
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Path to results_consensus.qs file [Required]", metavar = "FILE"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output directory for plots (default: same as input directory)", metavar = "DIR"
    ),
    make_option(c("--cluster"),
        type = "character", default = "Unknown_Cluster",
        help = "Cluster name for caption [default= %default]"
    ),
    make_option(c("--group"),
        type = "character", default = "Unknown_Group",
        help = "Group/Condition name for caption [default= %default]"
    ),
    make_option(c("--covariates"),
        type = "character", default = "Unknown",
        help = "Covariates string for caption [default= %default]"
    ),
    make_option(c("--random_effect"),
        type = "character", default = "hos_no",
        help = "Random effect variable [default= %default]"
    ),
    make_option(c("--src"),
        type = "character", default = "myR/R/deg_consensus",
        help = "Directory containing R source files [default= %default]"
    ),
    make_option(c("--renv"),
        type = "character", default = NULL,
        help = "Path to renv library to use (optional)"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# --- Validation ---
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file argument (--input) is required.")
}
if (!file.exists(opt$input)) {
    stop("Input file not found: ", opt$input)
}

# --- Setup Environment ---
if (!is.null(opt$renv) && dir.exists(opt$renv)) {
    .libPaths(c(opt$renv, .libPaths()))
}

# Source Functions
if (dir.exists(opt$src)) {
    src_files <- list.files(opt$src, pattern = "\\.R$", full.names = TRUE)
    for (f in src_files) source(f)
} else {
    if (dir.exists("myR")) devtools::load_all("myR")
}

# --- Main Execution ---
res_path <- opt$input
out_dir <- if (is.null(opt$output)) dirname(res_path) else opt$output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading result: ", res_path)
res <- qs::qread(res_path)

# Extract Data
cons_scores <- res$consensus_scores
mat <- res$deg_matrices
pca_res <- res$pca_result
if (is.null(pca_res) || is.null(pca_res$coordinates)) {
    if (!is.null(mat)) pca_res <- perform_deg_pca(mat)
}

# Construct Caption
get_git_hash <- function() {
    tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")
}
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
git_hash <- get_git_hash()

# Parse covariates
cov_list <- unlist(strsplit(opt$covariates, ","))
cov_str <- paste(cov_list, collapse = " + ")

caption_text <- paste0(
    "Cluster: ", opt$cluster, "\n",
    "Formula: ~ ", opt$group, " + ", cov_str, " + (1|", opt$random_effect, ") [Varies]\n",
    "Script: plot_consensus_cli.R\n",
    "Time: ", timestamp, "\n",
    "Commit: ", git_hash
)

plot_prefix <- file.path(out_dir, "AG_Consensus_Replot")
message("Generating plots in: ", out_dir)

# 1. Plots
try({
    plot_consensus_meta_distributions(cons_scores, file_prefix = plot_prefix, caption = caption_text)
    plot_consensus_volcano(cons_scores, file_prefix = plot_prefix, fdr_threshold = 0.05, caption = caption_text)
    plot_consensus_heatmap(mat, cons_scores, file_prefix = plot_prefix, method_top_n = 20, caption = caption_text)

    if (!is.null(pca_res)) plot_method_pca(pca_res, file_prefix = plot_prefix, caption = caption_text)
    plot_gene_umap(mat, cons_scores, file_prefix = plot_prefix, caption = caption_text)

    # New: Ridge Plot
    plot_pvalue_ridges(res, file_prefix = plot_prefix, caption = caption_text)
})

# 2. Similarity
if (!is.null(mat$statistic)) {
    mat_stat <- mat$statistic
    valid_rows <- rowSums(!is.na(mat_stat)) >= 2
    if (sum(valid_rows) > 10) {
        mat_stat_sub <- mat_stat[valid_rows, ]

        # Rename columns to include statistic type for legend
        new_cols <- sapply(colnames(mat_stat_sub), function(x) {
            suffix <- "stat"
            if (grepl("LRT", x)) {
                suffix <- "LRT"
            } else if (grepl("Wald", x)) {
                suffix <- "Wald"
            } else if (grepl("QLF", x)) {
                suffix <- "F-test"
            } else if (grepl("voom|trend|dream", x)) {
                suffix <- "t-stat"
            } else if (grepl("nebula", x)) suffix <- "stat"
            paste0(x, " (", suffix, ")")
        })
        colnames(mat_stat_sub) <- new_cols

        cor_mat <- cor(mat_stat_sub, use = "pairwise.complete.obs", method = "spearman")

        p_file <- file.path(out_dir, "method_similarity_heatmap.png")
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

# 3. CSV Export


# 2-B. Similarity: LogFC and P-value based
if (!is.null(mat$beta) && !is.null(mat$pvalue)) {
    # LOGFC Similarity
    mat_beta <- mat$beta
    valid_beta <- rowSums(!is.na(mat_beta)) >= 2
    if (sum(valid_beta) > 10) {
        cor_mat_beta <- cor(mat_beta[valid_beta, ], use = "pairwise.complete.obs", method = "spearman")
        p_file <- file.path(out_dir, "method_similarity_logFC_spearman.png")
        png(p_file, width = 800, height = 800)
        pheatmap::pheatmap(cor_mat_beta,
            main = paste("Method Similarity (LogFC Spearman) -", opt$cluster),
            display_numbers = TRUE, treeheight_row = 0, fontsize = 10
        )
        dev.off()
    }

    # P-value Similarity (using Rank)
    mat_pval <- mat$pvalue
    valid_pval <- rowSums(!is.na(mat_pval)) >= 2
    if (sum(valid_pval) > 10) {
        cor_mat_pval <- cor(mat_pval[valid_pval, ], use = "pairwise.complete.obs", method = "spearman")
        p_file <- file.path(out_dir, "method_similarity_pvalue_spearman.png")
        png(p_file, width = 800, height = 800)
        pheatmap::pheatmap(cor_mat_pval,
            main = paste("Method Similarity (P-value Spearman) -", opt$cluster),
            display_numbers = TRUE, treeheight_row = 0, fontsize = 10
        )
        dev.off()
    }
}

# 3. CSV Export and Image Summary and Text Summary
export_method_summary_csv(res, file.path(out_dir, "method_summary.csv"))
generate_summary_table_image(res, file.path(out_dir, "method_summary_table.png"))
export_method_formulas_text(res, file.path(out_dir, "method_formulas.txt"))

message("Re-plotting Complete.")
