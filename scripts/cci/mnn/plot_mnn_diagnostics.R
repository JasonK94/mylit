#!/usr/bin/env Rscript
# Plot MultiNicheNet Diagnostics (Bubble Plots & Heatmaps)
# Usage: Rscript plot_mnn_diagnostics.R -f <results.qs> -o <outdir>

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
    library(multinichenetr)
    library(qs)
})

# Source Utilities
if (file.exists("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_cci.R")) {
    source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_cci.R")
}

option_list <- list(
    make_option(c("-f", "--file"),
        type = "character", default = NULL,
        help = "Path to MultiNicheNet results file (.qs)"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = "mnn_plots",
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
mnn_results <- qs::qread(opt$file)

cat("Generating Diagnostic Plots...\n")

# Check if prioritization tables exist
if (!is.null(mnn_results$prioritization_tables)) {
    # 1. Bubble Plot: Top Interactions
    # Typically this involves make_sample_lr_prod_activity_plots or make_bubble_plot
    # We will try to make a bubble plot of the top interactions grouped by contrast/group

    # Extract top interactions
    prior_tbl <- mnn_results$prioritization_tables$group_prioritization_tbl

    if (!is.null(prior_tbl)) {
        # Filter top 50
        top_tbl <- prior_tbl %>%
            arrange(desc(score)) %>%
            head(50)

        # Plotting L-R expression (Bubble Plot)
        # Note: This usually requires the original SCE/Seurat for expression values
        # effectively, mnn_results might store some expression averages if configured.

        # If not, we can only plot the scores

        p_scores <- ggplot(top_tbl, aes(x = id, y = score, fill = group)) +
            geom_col() +
            coord_flip() +
            theme_classic() +
            labs(title = "Top 50 Interactions (MNN Score)", x = "Interaction", y = "Score")

        save_plot_formats("mnn_top_scores", opt$outdir, 8, 12, plot_obj = p_scores)
    }

    # 2. Ligand Activities Heatmap
    if (!is.null(mnn_results$ligand_activities_target_de_tbl)) {
        # Aggregate Ligand Activities
        la_tbl <- mnn_results$ligand_activities_target_de_tbl %>%
            distinct(ligand, activity) %>%
            arrange(desc(activity)) %>%
            head(30)

        p_la <- ggplot(la_tbl, aes(x = reorder(ligand, activity), y = activity)) +
            geom_col(fill = "steelblue") +
            coord_flip() +
            theme_classic() +
            labs(title = "Top Ligand Activities", x = "Ligand", y = "Activity (Pearson)")

        save_plot_formats("mnn_ligand_activities", opt$outdir, 6, 8, plot_obj = p_la)
    }
} else {
    cat("Warning: Prioritization tables not found in results.\n")
}

# 3. Combined PDF
cat("Saving combined PDF...\n")
combined_pdf <- file.path(opt$outdir, "mnn_diagnostics_combined.pdf")
pdf(combined_pdf, width = 10, height = 12)

if (exists("p_scores")) print(p_scores)
if (exists("p_la")) print(p_la)

dev.off()
cat("Combined PDF saved to:", combined_pdf, "\n")

cat("✅ Diagnostic Plots Generated!\n")
