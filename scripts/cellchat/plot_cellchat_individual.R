#!/usr/bin/env Rscript

# scripts/cellchat/run_cellchat_plots.R
# Script to generate comparison plots from merged CellChat objects

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
plot_func_path <- "myR/R/cci_cellchat_plotting.R"
if (!file.exists(plot_func_path)) {
    plot_func_path <- "/home/user3/data_user3/git_repo/_wt/cellchat/myR/R/cci_cellchat_plotting.R"
}
if (file.exists(plot_func_path)) {
    source(plot_func_path)
} else {
    stop("Cannot find plotting function: ", plot_func_path)
}

# Options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Path to MERGED CellChat object (.qs) or directory containing them [Required]"
    ),
    make_option(c("-o", "--output_dir"),
        type = "character", default = NULL,
        help = "Output directory for plots [Required]"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output_dir)) {
    print_help(opt_parser)
    stop("--input and --output_dir are required")
}

# Check input
if (dir.exists(opt$input)) {
    # If input is directory, look for cellchat_merged.qs
    qs_files <- list.files(opt$input, pattern = "cellchat_merged.qs", recursive = TRUE, full.names = TRUE)
    if (length(qs_files) == 0) {
        stop("No cellchat_merged.qs files found in directory: ", opt$input)
    }
} else if (file.exists(opt$input)) {
    qs_files <- c(opt$input)
} else {
    stop("Input not found: ", opt$input)
}

message("Found ", length(qs_files), " merged object(s) to process.")

# Process each merged object
for (qs_file in qs_files) {
    message("\nProcessing: ", qs_file)

    # Load object
    cc_merged <- qs::qread(qs_file)

    # Determine logical output subdirectory
    # If recursing, keep structure
    stats <- file.info(qs_file)
    message("  Size: ", round(stats$size / 1024 / 1024, 2), " MB")

    # Check if it's actually a merged object (list of cellchats inside or merged struct)
    # CellChat mergeCellChat returns a CellChat object, but with
    # @meta$datasets indicating original datasets if merged properly.
    # WAIT: aggregate.by implementation produces ONE merged object per condition (X1, X2).
    # BUT compareInteractions needs a LIST of CellChat objects (one per condition).

    # The user request said: "X2, X1 각각의 merged.qs에 대하여"
    # But comparison plots usually compare X1 vs X2.
    # IF we just want to plot X1 alone and X2 alone, that's fine.
    # Comparison functions (netVisual_diffInteraction) require a MERGED object containing BOTH conditions.
    # OR a list of two objects.

    # Let's check what we have.
    # X1/cellchat_merged.qs = Merged 8 control patients (Condition 1 Aggregated)
    # X2/cellchat_merged.qs = Merged 15 stroke patients (Condition 2 Aggregated)

    # The comparison needs to JOIN these two.

    # For now, let's implement PLOTTING FOR SINGLE CONDITION (since input is single file)
    # Standard plots for aggregated object: Circle, Bubble, NetVisual

    # Process merged object
    # If net$count is missing (common after merge), run aggregateNet to consolidate
    if (is.null(cc_merged@net$count)) {
        message("  Running aggregateNet to consolidate merged interactions...")
        cc_merged <- tryCatch(
            {
                CellChat::aggregateNet(cc_merged)
            },
            error = function(e) {
                message("  WARNING: aggregateNet failed (likely no interactions detected): ", e$message)
                return(NULL)
            }
        )
    }

    if (is.null(cc_merged) || is.null(cc_merged@net$count)) {
        message("  Skipping plotting for ", basename(qs_file), " - No aggregated interactions found.")
        next
    }

    # Define sub-output dir based on file name/parent dir to avoid overwrite
    parent_name <- basename(dirname(qs_file))
    plot_out_dir <- file.path(opt$output_dir, parent_name)
    if (!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE)

    message("  Plotting individual plots to: ", plot_out_dir)

    # Handle idents (merged object has list)
    idents_to_use <- cc_merged@idents
    if (is.list(idents_to_use) && "joint" %in% names(idents_to_use)) {
        idents_to_use <- idents_to_use$joint
    }

    # Circle Plot
    pdf(file.path(plot_out_dir, "net_circle_aggregated.pdf"), width = 10, height = 10)
    groupSize <- as.numeric(table(idents_to_use))
    CellChat::netVisual_circle(cc_merged@net$count,
        vertex.weight = groupSize,
        weight.scale = T, label.edge = F,
        title.name = paste0("Interactions (", parent_name, ")")
    )
    dev.off()

    # Bubble Plot
    if (length(cc_merged@netP$pathways) > 0) {
        pdf(file.path(plot_out_dir, "net_bubble_aggregated.pdf"), width = 12, height = 8)
        tryCatch(
            {
                print(CellChat::netVisual_bubble(cc_merged, remove.isolate = FALSE))
            },
            error = function(e) message("  Bubble error: ", e$message)
        )
        dev.off()
    }
}

message("\nDone.")
