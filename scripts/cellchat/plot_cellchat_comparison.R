#!/usr/bin/env Rscript

# scripts/cellchat/plot_cellchat_comparison.R
# Script to generate comparison plots between two conditions

# Set library path
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
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
if (!file.exists(plot_func_path)) stop("Cannot find plotting function: ", plot_func_path)
source(plot_func_path)

# Options
option_list <- list(
    make_option(c("--file1"), type = "character", help = "Path to FIRST CellChat object (.qs) [Control/Baseline]"),
    make_option(c("--file2"), type = "character", help = "Path to SECOND CellChat object (.qs) [Condition/Treatment]"),
    make_option(c("--name1"), type = "character", default = "Group1", help = "Name for First Group"),
    make_option(c("--name2"), type = "character", default = "Group2", help = "Name for Second Group"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for comparison plots")
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

# Check and fix aggregateNet if missing
if (is.null(cc1@net$count)) {
    message("Running aggregateNet for Group 1...")
    cc1 <- CellChat::aggregateNet(cc1)
}
if (is.null(cc2@net$count)) {
    message("Running aggregateNet for Group 2...")
    cc2 <- CellChat::aggregateNet(cc2)
}

# 1. Update/Lift objects if needed (CellChat::updateCellChat)
cc1 <- CellChat::updateCellChat(cc1)
cc2 <- CellChat::updateCellChat(cc2)

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
