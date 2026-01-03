#!/usr/bin/env Rscript

# scripts/cellchat/run_cellchat_comparison_analysis.R

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
    library(patchwork)
    library(ComplexHeatmap)
    library(qs)
})

option_list <- list(
    make_option(c("--file1"), type = "character", help = "Path to first CellChat object [Required]"),
    make_option(c("--file2"), type = "character", help = "Path to second CellChat object [Required]"),
    make_option(c("--name1"), type = "character", default = "Control", help = "Name for first group"),
    make_option(c("--name2"), type = "character", default = "Stroke", help = "Name for second group"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output directory [Required]"),
    make_option(c("--collapse_replicates"), action = "store_true", default = TRUE, help = "Collapse replicates into single group for 1-vs-1 comparison plots [Default: TRUE]"),
    make_option(c("--renv"), type = "character", default = NULL, help = "Path to renv library (optional)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$file1) || is.null(opt$file2) || is.null(opt$output_dir)) {
    stop("Missing arguments.")
}

message("Loading objects...")
cc1 <- qs::qread(opt$file1)
cc2 <- qs::qread(opt$file2)

# Function to collapse a merged object (list of replicates) into a single representative object
collapse_cellchat <- function(cc, name) {
    # If net is already a single matrix (not list of matrices), it's already collapsed/single
    if (!is.null(cc@net$count) && is.matrix(cc@net$count)) {
        return(cc)
    }

    # If net is a list (replicates), we need to aggregate
    if (is.list(cc@net) && !is.null(names(cc@net))) {
        message(paste0("  Collapsing ", length(names(cc@net)), " replicates in ", name, " into one..."))

        # 1. Identify union of all cell types
        all_cells <- unique(unlist(lapply(cc@net, function(x) rownames(x$count))))
        all_cells <- sort(all_cells)
        n_cells <- length(all_cells)

        message(paste0("    Union of cell types: ", n_cells))

        # 2. Lift and Sum
        net_count_sum <- matrix(0, nrow = n_cells, ncol = n_cells, dimnames = list(all_cells, all_cells))
        weight_sum <- matrix(0, nrow = n_cells, ncol = n_cells, dimnames = list(all_cells, all_cells))

        # Helper to lift a matrix
        lift_matrix <- function(mat, cells) {
            new_mat <- matrix(0, nrow = length(cells), ncol = length(cells), dimnames = list(cells, cells))
            common <- intersect(rownames(mat), cells)
            if (length(common) > 0) {
                new_mat[common, common] <- mat[common, common]
            }
            return(new_mat)
        }

        for (n in names(cc@net)) {
            # Handle NAs
            cnt <- cc@net[[n]]$count
            cnt[is.na(cnt)] <- 0
            wt <- cc@net[[n]]$weight
            wt[is.na(wt)] <- 0

            # Lift
            cnt_lifted <- lift_matrix(cnt, all_cells)
            wt_lifted <- lift_matrix(wt, all_cells)

            net_count_sum <- net_count_sum + cnt_lifted
            weight_sum <- weight_sum + wt_lifted
        }

        # Mean for weight/prob
        weight_mean <- weight_sum / length(names(cc@net))

        # Create a new object structure (shallow copy)
        new_cc <- cc
        new_cc@net <- list(
            count = net_count_sum,
            weight = weight_mean,
            prob = weight_mean, # Approximate prob as weight for visualization
            pval = weight_mean # Dummy
        )

        # Note: We do NOT subset meta or data.signaling here because
        # all_cells contains Cluster Names (e.g. 'Mono'), while meta contains Barcodes.
        # The original object should already have valid meta/data for all cells.
        # mergeCellChat needs the full meta/data to verify barcodes if it checks them.

        return(new_cc)
    }
    return(cc)
}

# Fix metadata barcode mismatch issue manually if present
fix_meta <- function(cc) {
    # Ensure data matrix colnames match meta rownames
    if (!identical(colnames(cc@data.signaling), rownames(cc@meta))) {
        # Check if they are just overlapping
        common <- intersect(colnames(cc@data.signaling), rownames(cc@meta))
        if (length(common) > 0) {
            cc@data.signaling <- cc@data.signaling[, common, drop = FALSE]
            cc@meta <- cc@meta[common, , drop = FALSE]
        }
    }
    return(cc)
}

cc1 <- fix_meta(cc1)
cc2 <- fix_meta(cc2)

# Collapse if requested
if (opt$collapse_replicates) {
    cc1_use <- collapse_cellchat(cc1, opt$name1)
    cc2_use <- collapse_cellchat(cc2, opt$name2)
} else {
    cc1_use <- cc1
    cc2_use <- cc2
}

if (is.null(cc1_use) || is.null(cc2_use)) {
    stop("Collapsing failed. check inputs.")
}

# Merge
message("Merging for comparison...")
object.list <- list()
object.list[[opt$name1]] <- cc1_use
object.list[[opt$name2]] <- cc2_use

# Use mergeCellChat
# Catch error if barcodes still mismatch
cellchat <- tryCatch(
    {
        mergeCellChat(object.list, add.names = names(object.list))
    },
    error = function(e) {
        message("Error in mergeCellChat: ", e$message)
        message("Attempting force merge using cell.prefix=TRUE...")
        # Force merge using cell.prefix = TRUE to avoid barcode collision errors
        mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
    }
)

# Create output dir
if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

message("Performing Downstream Analysis...")

# --- A. General Stats Comparison ---
pdf(file.path(opt$output_dir, "comparison_interaction_stats.pdf"), width = 8, height = 6)
try({
    p1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
    p2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
    print(p1 + p2)
})
dev.off()

# --- B. Differential Interaction Network (Circle Plots) ---
pdf(file.path(opt$output_dir, "comparison_diffInteraction_circle.pdf"), width = 10, height = 10)
par(mfrow = c(1, 2))
try({
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
})
dev.off()

# Heatmap
pdf(file.path(opt$output_dir, "comparison_diffInteraction_heatmap.pdf"), width = 10, height = 10)
par(mfrow = c(1, 2))
try({
    netVisual_heatmap(cellchat)
    netVisual_heatmap(cellchat, measure = "weight")
})
dev.off()

# --- C. Information Flow ---
pdf(file.path(opt$output_dir, "comparison_information_flow.pdf"), width = 10, height = 12)
try({
    rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
})
dev.off()

# --- D. Bubble Plots ---
pdf(file.path(opt$output_dir, "comparison_bubble_upregulated.pdf"), width = 12, height = 10)
try({
    # Simple bubble plot of all pathways (if not too huge)
    print(netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE))
})
dev.off()

message("Done. Results saved to: ", opt$output_dir)
