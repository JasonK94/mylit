#!/usr/bin/env Rscript

# CellChat CLI with Proper Aggregation Method

# Set library path
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
}

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(CellChat)
    library(future)
    library(dplyr)
    if (requireNamespace("qs", quietly = TRUE)) library(qs)
})

options(future.globals.maxSize = 200 * 1024^3)

# Options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Input Seurat object (.qs/.rds) [Required]"
    ),
    make_option(c("-g", "--group_by"),
        type = "character", default = NULL,
        help = "Cell type grouping column [Required]"
    ),
    make_option(c("-s", "--split_by"),
        type = "character", default = NULL,
        help = "Sample-level split (e.g., hos_no for patients). Each gets independent analysis. [Recommended]"
    ),
    make_option(c("-a", "--aggregate_by"),
        type = "character", default = NULL,
        help = "Condition for aggregation (e.g., g3). Merges split samples by this. [Optional]"
    ),
    make_option(c("--subset_aggregate"),
        type = "character", default = NULL,
        help = "Comma-separated conditions to analyze (e.g., '2,1') [Optional]"
    ),
    make_option(c("-d", "--db_use"),
        type = "character", default = NULL,
        help = "DB types (comma-separated): 'Secreted Signaling', 'ECM-Receptor', 'Cell-Cell Contact'. Runs each separately. [Optional, uses all if not specified]"
    ),
    make_option(c("-p", "--prob_threshold"),
        type = "numeric", default = 0.05,
        help = "Probability threshold [default: %default]"
    ),
    make_option(c("-m", "--min_cells"),
        type = "integer", default = 10,
        help = "Minimum cells per cell type [default: %default]"
    ),
    make_option(c("--species"),
        type = "character", default = "human",
        help = "Species [default: %default]"
    ),
    make_option(c("--assay"),
        type = "character", default = "RNA",
        help = "Assay [default: %default]"
    ),
    make_option(c("-c", "--cores"),
        type = "integer", default = 8,
        help = "Cores [default: %default]"
    ),
    make_option(c("-o", "--output_dir"),
        type = "character", default = NULL,
        help = "Output directory [Optional]"
    ),
    make_option(c("-n", "--name"),
        type = "character", default = NULL,
        help = "Run name [Optional, auto-generated]"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$group_by)) {
    print_help(opt_parser)
    stop("--input and --group_by are required")
}

# Parse subset_aggregate
subset_agg <- NULL
if (!is.null(opt$subset_aggregate)) {
    subset_agg <- strsplit(opt$subset_aggregate, ",")[[1]]
    subset_agg <- trimws(subset_agg)
    subset_agg_num <- suppressWarnings(as.numeric(subset_agg))
    if (!any(is.na(subset_agg_num))) subset_agg <- subset_agg_num
}

# Parse db_use
db_use_val <- NULL
if (!is.null(opt$db_use)) {
    db_use_val <- strsplit(opt$db_use, ",")[[1]]
    db_use_val <- trimws(db_use_val)
}

# Source wrapper
wrapper_path <- "myR/R/cci_cellchat_wrapper.R"
if (!file.exists(wrapper_path)) {
    wrapper_path <- "/home/user3/data_user3/git_repo/_wt/cellchat/myR/R/cci_cellchat_wrapper.R"
}
if (file.exists(wrapper_path)) {
    source(wrapper_path)
} else {
    stop("Cannot find wrapper")
}

# Run
message("=== CellChat CLI (Proper Aggregation Method) ===")
message("Input: ", opt$input)
message("Group by: ", opt$group_by)
if (!is.null(opt$split_by)) {
    message("Split by (sample-level): ", opt$split_by)
}
if (!is.null(opt$aggregate_by)) {
    message("Aggregate by (condition): ", opt$aggregate_by)
    if (!is.null(subset_agg)) {
        message("  Subset: ", paste(subset_agg, collapse = ", "))
    }
}
if (!is.null(db_use_val)) {
    message("DB types (separate analyses): ", paste(db_use_val, collapse = ", "))
}
message("Probability threshold: ", opt$prob_threshold)

run_cellchat_analysis(
    input_data = opt$input,
    group.by = opt$group_by,
    split.by = opt$split_by,
    aggregate.by = opt$aggregate_by,
    subset.aggregate = subset_agg,
    db.use = db_use_val,
    species = opt$species,
    assay_name = opt$assay,
    min.cells.group = opt$min_cells,
    prob.threshold = opt$prob_threshold,
    output_dir = opt$output_dir,
    n_cores = opt$cores,
    run_name = opt$name,
    verbose = TRUE
)

message("Done.")
