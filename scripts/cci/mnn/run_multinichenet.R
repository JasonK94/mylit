#!/usr/bin/env Rscript
# Generic MultiNicheNet Analysis Script
# Usage: Rscript run_multinichenet.R -i <data_path> -g <group_id> -s <sample_id> -c <celltype_id> -f <contrasts>

# Set library path
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
} else {
    warning("Shared renv library not found at: ", renv_lib)
}

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(dplyr)
    library(qs)
})

# Define command line options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = "/data/user3/sobj/is2_IS_3_1_plots.qs",
        help = "Input Seurat object path [default: %default]"
    ),
    make_option(c("-s", "--sample_id"),
        type = "character", default = "hos_no",
        help = "Sample ID column name [default: %default]"
    ),
    make_option(c("-g", "--group_id"),
        type = "character", default = "g3",
        help = "Group/condition column name [default: %default]"
    ),
    make_option(c("-c", "--celltype_id"),
        type = "character", default = "anno3.scvi",
        help = "Cell type column name [default: %default]"
    ),
    make_option(c("-f", "--contrasts"),
        type = "character", default = "X2-X1",
        help = "Contrasts formula (after make.names) [default: %default]"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "multinichenet_results",
        help = "Output directory name [default: %default]"
    ),
    make_option("--min_cells",
        type = "integer", default = 10,
        help = "Minimum cells per cell type per sample [default: %default]"
    ),
    make_option("--logfc_threshold",
        type = "double", default = 0.10,
        help = "Log fold-change threshold for DE [default: %default]"
    ),
    make_option("--fraction_cutoff",
        type = "double", default = 0.05,
        help = "Fraction cutoff for gene expression (min.pct) [default: %default]"
    ),
    make_option("--p_val_threshold",
        type = "double", default = 0.05,
        help = "P-value threshold for DE [default: %default]"
    ),
    make_option("--p_val_adj",
        type = "logical", default = FALSE,
        help = "Use adjusted p-value [default: %default]"
    ),
    make_option("--cores",
        type = "integer", default = 8,
        help = "Number of cores to use [default: %default]"
    ),
    make_option("--senders",
        type = "character", default = NULL,
        help = "Comma-separated sender cell types (NULL = all) [default: %default]"
    ),
    make_option("--receivers",
        type = "character", default = NULL,
        help = "Comma-separated receiver cell types (NULL = all) [default: %default]"
    ),
    make_option("--verbose",
        type = "logical", default = TRUE,
        help = "Verbose output [default: %default]"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║         MultiNicheNet Analysis - Generic Script           ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("Parameters:\n")
cat("  Input data:        ", opt$input, "\n")
cat("  Sample ID:         ", opt$sample_id, "\n")
cat("  Group ID:          ", opt$group_id, "\n")
cat("  Cell type ID:      ", opt$celltype_id, "\n")
cat("  Contrasts:         ", opt$contrasts, "\n")
cat("  Output directory:  ", opt$output, "\n")
cat("  Min cells:         ", opt$min_cells, "\n")
cat("  LogFC threshold:   ", opt$logfc_threshold, "\n")
cat("  Fraction cutoff:   ", opt$fraction_cutoff, "\n")
cat("  P-value threshold: ", opt$p_val_threshold, "\n")
cat("  P-value adjusted:  ", opt$p_val_adj, "\n")
cat("  Cores:             ", opt$cores, "\n")
cat("\n")

# Source wrapper function
wrapper_path <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci_multinichenet_wrapper.R"
if (!file.exists(wrapper_path)) {
    stop("Wrapper function not found: ", wrapper_path)
}
source(wrapper_path)

# Load data
cat("Loading data...\n")
if (!file.exists(opt$input)) {
    stop("Input file not found: ", opt$input)
}
sobj <- qs::qread(opt$input)
cat("  Total cells: ", ncol(sobj), "\n")
cat("  Total features: ", nrow(sobj), "\n\n")

# Parse sender/receiver lists
senders_oi <- if (!is.null(opt$senders)) {
    strsplit(opt$senders, ",")[[1]]
} else {
    NULL
}

receivers_oi <- if (!is.null(opt$receivers)) {
    strsplit(opt$receivers, ",")[[1]]
} else {
    NULL
}

# Run analysis
cat("Starting MultiNicheNet analysis...\n\n")
start_time <- Sys.time()

results <- tryCatch(
    {
        run_multinichenet_analysis(
            sobj = sobj,
            sample_id = opt$sample_id,
            group_id = opt$group_id,
            celltype_id = opt$celltype_id,
            contrasts_oi = opt$contrasts,
            senders_oi = senders_oi,
            receivers_oi = receivers_oi,
            min_cells = opt$min_cells,
            species = "human",
            output_dir = opt$output,
            verbose = opt$verbose,
            cores = opt$cores,
            # Threshold parameters
            min_pct = opt$fraction_cutoff,
            logfc_thresh = opt$logfc_threshold,
            p_val_thresh = opt$p_val_threshold,
            p_val_adj = opt$p_val_adj
        )
    },
    error = function(e) {
        cat("\n❌ Analysis failed:\n")
        cat(e$message, "\n")
        quit(status = 1)
    }
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║                    Analysis Complete                       ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("Elapsed time: ", round(elapsed, 2), " minutes\n")
cat("Results saved to: ", opt$output, "/multinichenet_results.rds\n")

# Summary
if (!is.null(results$prioritization_tables)) {
    tbl <- results$prioritization_tables$group_prioritization_tbl
    cat("\nPrioritization table: ", nrow(tbl), " interactions\n")

    if (nrow(tbl) > 0) {
        cat("  Unique senders: ", length(unique(tbl$sender)), "\n")
        cat("  Unique receivers: ", length(unique(tbl$receiver)), "\n")
        cat("  Unique ligands: ", length(unique(tbl$ligand)), "\n")
        cat("  Unique receptors: ", length(unique(tbl$receptor)), "\n")
    } else {
        cat("\n⚠️  Warning: Prioritization table is empty.\n")
        cat("   This may indicate:\n")
        cat("   - No significant interactions found\n")
        cat("   - Thresholds too stringent\n")
        cat("   - Insufficient statistical power\n")
    }
}

cat("\nDE genes: ", nrow(results$celltype_de), "\n")
cat("Ligand activities: ", nrow(results$ligand_activities_targets_DEgenes$ligand_activities), "\n")

cat("\n✅ Done.\n")
