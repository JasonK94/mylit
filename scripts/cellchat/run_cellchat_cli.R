#!/usr/bin/env Rscript

# scripts/cellchat/run_cellchat_cli.R
# Command line interface for CellChat analysis

# Set library path to shared renv
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
} else {
    warning("Shared renv library not found at: ", renv_lib)
}

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(CellChat)
    library(future)
    if (requireNamespace("qs", quietly = TRUE)) {
        library(qs)
    }
})

# Define options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Path to Seurat object (.qs or .rds) [Required]", metavar = "file"
    ),
    make_option(c("-g", "--group_by"),
        type = "character", default = NULL,
        help = "Metadata column for cell grouping [Required]", metavar = "column"
    ),
    make_option(c("-o", "--output_dir"),
        type = "character", default = NULL,
        help = "Output directory [Optional, defaults to docs/cellchat/runX]", metavar = "dir"
    ),
    make_option(c("-s", "--species"),
        type = "character", default = "human",
        help = "Species (human or mouse) [default: %default]", metavar = "species"
    ),
    make_option(c("-a", "--assay"),
        type = "character", default = "RNA",
        help = "Assay name [default: %default]", metavar = "assay"
    ),
    make_option(c("-c", "--cores"),
        type = "integer", default = 16,
        help = "Number of cores [default: %default]", metavar = "int"
    ),
    make_option(c("-n", "--name"),
        type = "character", default = NULL,
        help = "Run name (e.g. run1) [Optional, auto-generated]", metavar = "name"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$group_by)) {
    print_help(opt_parser)
    stop("Input file and group_by column are required.", call. = FALSE)
}

# Source the wrapper function
# Assuming the script is run from the project root or we can find the file
wrapper_path <- "myR/R/cci_cellchat_wrapper.R"
if (!file.exists(wrapper_path)) {
    # Try absolute path based on known structure
    wrapper_path <- "/home/user3/data_user3/git_repo/_wt/cellchat/myR/R/cci_cellchat_wrapper.R"
}

if (file.exists(wrapper_path)) {
    source(wrapper_path)
} else {
    stop("Could not find cci_cellchat_wrapper.R")
}

# Run analysis
message("Starting CellChat CLI...")
message("Input: ", opt$input)
message("Group By: ", opt$group_by)

# The wrapper handles loading the input file
run_cellchat_analysis(
    input_data = opt$input,
    group.by = opt$group_by,
    species = opt$species,
    assay_name = opt$assay,
    output_dir = opt$output_dir,
    n_cores = opt$cores,
    run_name = opt$name,
    verbose = TRUE
)

message("Done.")
