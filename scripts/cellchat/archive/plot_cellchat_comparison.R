#!/usr/bin/env Rscript

# scripts/cellchat/plot_cellchat_comparison.R
# Script to generate comparison plots between two conditions

# --- Pre-load for renv ---
args <- commandArgs(trailingOnly = FALSE)
renv_arg_idx <- grep("--renv=", args)
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu" # Default

if (length(renv_arg_idx) > 0) {
    renv_lib <- sub("--renv=", "", args[renv_arg_idx])
}

if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
} else {
    warning("Library path not found: ", renv_lib)
}

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(CellChat)
    library(ggplot2)
    library(qs)
})

# Source plotting function
plot_func_path <- "/home/user3/data_user3/git_repo/_wt/cellchat/myR/R/cci_cellchat_plotting.R"
if (!file.exists(plot_func_path)) {
    # Fallback absolute path if relative fails
    plot_func_path <- "/home/user3/data_user3/git_repo/_wt/cellchat/myR/R/cci_cellchat_plotting.R"
}
if (!file.exists(plot_func_path)) stop("Cannot find plotting function: ", plot_func_path)
source(plot_func_path)

# Options
option_list <- list(
    make_option(c("--file1"), type = "character", help = "Path to FIRST CellChat object (.qs) [Control/Baseline]"),
    make_option(c("--file2"), type = "character", help = "Path to SECOND CellChat object (.qs) [Condition/Treatment]"),
    make_option(c("--name1"), type = "character", default = "Group1", help = "Name for First Group"),
    make_option(c("--name2"), type = "character", default = "Group2", help = "Name for Second Group"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for comparison plots"),
    make_option(c("--renv"), type = "character", default = NULL, help = "Path to renv library (optional)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file1) || is.null(opt$file2) || is.null(opt$output_dir)) {
    print_help(opt_parser)
    stop("file1, file2, and output_dir are required")
}

message("=== CellChat Comparison Plotting ===")
message("Group 1: ", opt$name1, " (", opt$file1, ")")
message("Group 2: ", opt$name2, " (", opt$file2, ")")

# Load objects
message("Loading objects...")
cc1 <- qs::qread(opt$file1)
cc2 <- qs::qread(opt$file2)

# Helper to validate object
validate_cc <- function(cc, name) {
    if (is.null(cc)) stop("Object ", name, " is NULL")

    # Try update first (handles v1->v2 slot changes)
    tryCatch(
        {
            cc <- CellChat::updateCellChat(cc)
        },
        error = function(e) {
            warning("updateCellChat failed for ", name, ": ", e$message)
        }
    )

    # Logic:
    # 1. If net is a list and has valid 'prob' in one of its elements -> It's valid (merged).
    # 2. If net has 'prob' directly -> It's valid (single).

    has_prob <- !is.null(cc@net$prob)
    if (!has_prob && is.list(cc@net)) {
        # Check if any element has prob, excluding 'pval', 'count', 'weight' if they are at top level
        # Actually in merged object, net$dataset_name$prob exists
        has_prob <- any(sapply(cc@net, function(x) is.list(x) && !is.null(x$prob)))
    }

    if (!has_prob) {
        warning("Object ", name, " has no interaction probability data (@net$prob is NULL or empty).")
        return(NULL)
    }

    # If count is missing, try aggregateNet (only for single objects or specific structure)
    if (is.null(cc@net$count) && !is.list(cc@net[[1]])) {
        tryCatch(
            {
                cc <- CellChat::aggregateNet(cc)
            },
            error = function(e) NULL
        )
    }

    return(cc)
}

cc1 <- validate_cc(cc1, opt$name1)
cc2 <- validate_cc(cc2, opt$name2)

if (is.null(cc1) || is.null(cc2)) {
    stop(
        "One or both CellChat objects contain no significant interactions.\n",
        "The input files likely come from a run where no interactions were detected (signal too weak).\n",
        "Suggestion: Re-run the analysis with more permissive parameters as per commands_KR.md:\n",
        "  --thresh_p_expression 0.1 --min_cells 5"
    )
}

# 2. Merge for comparison
message("Merging objects for comparison...")
object.list <- list(cc1, cc2)
names(object.list) <- c(opt$name1, opt$name2)

# Note: mergeCellChat for comparison takes a LIST of objects
cc_merged <- CellChat::mergeCellChat(object.list, add.names = names(object.list))

# 3. Generate Comparison Plots
message("Generating plots...")
plot_cellchat_comparison(cc_merged, opt$output_dir)

message("Done.")
