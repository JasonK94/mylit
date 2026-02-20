#!/usr/bin/env Rscript
# Generic NicheNet Analysis Script
# Usage: Rscript run_nichenet.R -i <sobj_path> ...

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(dplyr)
    library(qs)
    library(nichenetr)
})

# Define command line options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Input Seurat object path (.qs or .rds)"
    ),
    make_option(c("-s", "--senders"),
        type = "character", default = NULL,
        help = "Comma-separated sender cell types"
    ),
    make_option(c("-r", "--receiver"),
        type = "character", default = NULL,
        help = "Receiver cell type"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "nichenet_results",
        help = "Output directory"
    ),
    make_option(c("--condition_col"),
        type = "character", default = NULL,
        help = "Metadata column for condition (e.g. 'group')"
    ),
    make_option(c("--condition_oi"),
        type = "character", default = NULL,
        help = "Condition of interest (e.g. 'Tumor')"
    ),
    make_option(c("--condition_ref"),
        type = "character", default = NULL,
        help = "Reference condition (e.g. 'Normal')"
    ),
    make_option(c("--cluster_col"),
        type = "character", default = "cell_type",
        help = "Metadata column for cell types [default: %default]"
    ),
    make_option(c("--assay"),
        type = "character", default = "SCT",
        help = "Assay to use [default: %default]"
    ),
    make_option(c("--top_n_ligands"),
        type = "integer", default = 20,
        help = "Number of top ligands [default: %default]"
    ),
    make_option(c("--top_n_targets"),
        type = "integer", default = 200,
        help = "Number of top targets per ligand [default: %default]"
    ),
    make_option(c("--species"),
        type = "character", default = "human",
        help = "Species (human/mouse) [default: %default]"
    ),
    make_option(c("--verbose"),
        type = "logical", default = TRUE,
        help = "Verbose output [default: %default]"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$senders) || is.null(opt$receiver)) {
    print_help(opt_parser)
    stop("Input arguments missing.")
}

# Source Functions
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/core_nichenet.R")

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║            NicheNet Analysis - Run Script                  ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# Load Data
cat("Loading Seurat object...\n")
if (grepl("\\.qs$", opt$input)) {
    sobj <- qs::qread(opt$input)
} else {
    sobj <- readRDS(opt$input)
}

senders <- strsplit(opt$senders, ",")[[1]]
receiver <- opt$receiver

# Save parameters log
params <- as.list(opt)
save_parameters_log(params, opt$output, prefix = "nichenet_params")

cat("Running NicheNet Analysis...\n")
cat("  Senders: ", paste(senders, collapse = ", "), "\n")
cat("  Receiver: ", receiver, "\n")
cat("  Condition: ", opt$condition_oi, " vs ", opt$condition_ref, "\n")

results <- run_nichenet_analysis(
    seurat_obj = sobj,
    species = opt$species,
    sender_celltypes = senders,
    receiver_celltype = receiver,
    cluster_col = opt$cluster_col,
    assay_name = opt$assay,
    receiver_DE_group_by = opt$condition_col,
    receiver_DE_ident1 = opt$condition_oi,
    receiver_DE_ident2 = opt$condition_ref,
    top_n_ligands = opt$top_n_ligands,
    top_n_targets_per_ligand = opt$top_n_targets,
    output_dir = opt$output,
    run_circos = TRUE,
    verbose = opt$verbose
)

# Save Results
if (!dir.exists(opt$output)) dir.create(opt$output, recursive = TRUE)
qs_path <- get_safe_filename(file.path(opt$output, "nichenet_results.qs"))
qs::qsave(results, qs_path)

cat("\n✅ NicheNet Analysis Complete!\n")
cat("Results saved to: ", qs_path, "\n")
