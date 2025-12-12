#!/usr/bin/env Rscript
# Plot NicheNet Results
# Usage: Rscript plot_nichenet.R -f <results.qs> -o <output_dir>

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(qs)
    library(cowplot)
    library(gridExtra)
})

# Source Utilities
if (file.exists("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_cci.R")) {
    source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_cci.R")
}

option_list <- list(
    make_option(c("-f", "--file"),
        type = "character", default = NULL,
        help = "Path to NicheNet results file (.qs)"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = "nichenet_plots",
        help = "Output directory"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)) {
    print_help(opt_parser)
    stop("Input file is required.")
}

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("Loading results from:", opt$file, "\n")
results <- qs::qread(opt$file)

cat("Generating and saving plots...\n")

# 1. Ligand Activity Histogram
if (!is.null(results$plot_ligand_activity_hist)) {
    save_plot_formats("ligand_activity_histogram", opt$outdir, 8, 6, plot_obj = results$plot_ligand_activity_hist)
}

# 2. Ligand AUPR Heatmap
if (!is.null(results$plot_ligand_aupr_heatmap)) {
    save_plot_formats("ligand_aupr_heatmap", opt$outdir, 6, 8, plot_obj = results$plot_ligand_aupr_heatmap)
}

# 3. Ligand-Target Network Heatmap
if (!is.null(results$plot_ligand_target_network)) {
    save_plot_formats("ligand_target_heatmap", opt$outdir, 12, 10, plot_obj = results$plot_ligand_target_network)
}

# 4. Ligand-Receptor Network Heatmap
if (!is.null(results$plot_ligand_receptor_network)) {
    save_plot_formats("ligand_receptor_heatmap", opt$outdir, 10, 8, plot_obj = results$plot_ligand_receptor_network)
}

# 5. Circos Plot (Recorded Plot)
if (!is.null(results$plot_circos)) {
    cat("Saving Circos plot...\n")
    # For recorded plots, we use replayPlot
    save_plot_formats("circos_plot", opt$outdir, 10, 10, plot_obj = results$plot_circos)

    # Try to save a version with legend if possible (often needs complex handling)
}

# 6. Dotplot of Top Ligands in Senders (if Seurat object is available, but usually it's not in the result object)
# Note: The result object doesn't contain the Seurat object to save space.
# We can only plot what is available.

# 7. Dotplot of Top Targets (using expression data if available in results?)
# results doesn't seem to store expression data directly, only inference results.

cat("✅ Plots saved to:", opt$outdir, "\n")
