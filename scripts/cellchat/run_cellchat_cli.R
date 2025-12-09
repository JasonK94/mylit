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

options(future.globals.maxSize = 200 * 1024^3)

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
    make_option(c("-b", "--split_by"),
        type = "character", default = NULL,
        help = "Metadata column to split analysis by (e.g. patient, condition) [Optional]", metavar = "column"
    ),
    make_option(c("--subset_split"),
        type = "character", default = NULL,
        help = "Comma-separated list of split groups to analyze (e.g. '2,1') [Optional, analyzes all if not specified]", metavar = "groups"
    ),
    make_option(c("-d", "--db_use"),
        type = "character", default = NULL,
        help = "Interaction DB subset: 'Secreted Signaling', 'ECM-Receptor', or 'Cell-Cell Contact'. Can specify multiple comma-separated. Uses all if not specified. [Optional]", metavar = "db_type"
    ),
    make_option(c("-m", "--min_cells"),
        type = "integer", default = 10,
        help = "Minimum cells per group [default: %default]", metavar = "int"
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

# Parse subset_split if provided
subset_split_val <- NULL
if (!is.null(opt$subset_split)) {
    subset_split_val <- strsplit(opt$subset_split, ",")[[1]]
    # Trim whitespace
    subset_split_val <- trimws(subset_split_val)
    # Try to convert to numeric if possible
    subset_split_numeric <- suppressWarnings(as.numeric(subset_split_val))
    if (!any(is.na(subset_split_numeric))) {
        subset_split_val <- subset_split_numeric
    }
}

# Parse db_use if provided
db_use_val <- NULL
if (!is.null(opt$db_use)) {
    db_use_val <- strsplit(opt$db_use, ",")[[1]]
    db_use_val <- trimws(db_use_val)
}

# Source the wrapper function
wrapper_path <- "myR/R/cci_cellchat_wrapper.R"
if (!file.exists(wrapper_path)) {
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
if (!is.null(opt$split_by)) {
    message("Split By: ", opt$split_by)
    if (!is.null(subset_split_val)) {
        message("Subset Split: ", paste(subset_split_val, collapse = ", "))
    }
}
if (!is.null(db_use_val)) {
    message("DB Use: ", paste(db_use_val, collapse = ", "))
}

run_cellchat_analysis(
    input_data = opt$input,
    group.by = opt$group_by,
    split.by = opt$split_by,
    subset.split = subset_split_val,
    db.use = db_use_val,
    min.cells.group = opt$min_cells,
    species = opt$species,
    assay_name = opt$assay,
    output_dir = opt$output_dir,
    n_cores = opt$cores,
    run_name = opt$name,
    verbose = TRUE
)

message("Done.")
