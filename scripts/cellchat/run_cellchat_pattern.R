#!/usr/bin/env Rscript

# scripts/cellchat/run_cellchat_pattern.R
# Communication Pattern Analysis using NMF (Non-negative Matrix Factorization)

# --- 0. Parse Arguments Pre-load for renv ---
args <- commandArgs(trailingOnly = FALSE)
renv_arg_idx <- grep("--renv=", args)
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu" # Default

if (length(renv_arg_idx) > 0) {
    renv_lib <- sub("--renv=", "", args[renv_arg_idx])
}

if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
    message("Using library path: ", renv_lib)
} else {
    warning("Library path not found: ", renv_lib)
}

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(CellChat)
    library(ggplot2)
    library(qs)
    library(NMF)
    library(ggalluvial)
})

# Options
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to CellChat object (.qs)"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for plots"),
    make_option(c("-k", "--n_patterns"), type = "integer", default = 5, help = "Number of patterns (latent topics) [default: 5]"),
    make_option(c("--estimate_k"), action = "store_true", default = FALSE, help = "Run selectK to estimate optimal k"),
    make_option(c("--renv"), type = "character", default = NULL, help = "Path to renv library (optional)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output_dir)) {
    print_help(opt_parser)
    stop("Input and output directory required.")
}

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

message("=== CellChat Pattern Analysis (NMF) ===")
message("Input: ", opt$input)

# Load Data
cc <- qs::qread(opt$input)

# --- Validation ---
if (is.null(cc@net$prob)) {
    stop(
        "Error: The CellChat object has no interaction probability matrix (@net$prob is NULL).\n",
        "  Likely causes:\n",
        "  1. No interactions were detected in the upstream analysis.\n",
        "  2. 'thresh.p.expression' was too strict (try 0.1).\n",
        "  3. 'min.cells' was too high."
    )
}

# Helper to save plots
save_plot <- function(filename, plot_func, w = 10, h = 8) {
    pdf_path <- file.path(opt$output_dir, paste0(filename, ".pdf"))
    png_path <- file.path(opt$output_dir, paste0(filename, ".png"))

    message("  Saving: ", pdf_path)
    # Use Cairo if available or default
    tryCatch(
        {
            pdf(pdf_path, width = w, height = h)
            plot_func()
            dev.off()
        },
        error = function(e) warning("Error saving PDF ", filename, ": ", e$message)
    )

    tryCatch(
        {
            png(png_path, width = w * 100, height = h * 100, res = 100)
            plot_func()
            dev.off()
        },
        error = function(e) warning("Error saving PNG ", filename, ": ", e$message)
    )
}

# 1. Estimate K if requested
if (opt$estimate_k) {
    message("Estimating optimal number of patterns (selectK)... this may take a while.")
    save_plot("selectK_outgoing", function() {
        CellChat::selectK(cc, pattern = "outgoing")
    })
    save_plot("selectK_incoming", function() {
        CellChat::selectK(cc, pattern = "incoming")
    })
}

# 2. Analyze Outgoing Patterns (Sender)
message("Analyzing Outgoing (Sender) Patterns...")
tryCatch(
    {
        # Check if k is valid for matrix size
        if (opt$n_patterns >= min(dim(cc@net$prob)[1:2])) {
            new_k <- max(2, min(dim(cc@net$prob)[1:2]) - 1)
            warning("Adjusting k from ", opt$n_patterns, " to ", new_k, " due to small matrix dimensionality.")
            opt$n_patterns <- new_k
        }

        cc <- CellChat::identifyCommunicationPatterns(cc, pattern = "outgoing", k = opt$n_patterns)

        save_plot("patterns_outgoing_river", function() {
            CellChat::netAnalysis_river(cc, pattern = "outgoing")
        }, w = 12, h = 8)
        save_plot("patterns_outgoing_dot", function() {
            CellChat::netAnalysis_dot(cc, pattern = "outgoing")
        }, w = 10, h = 8)
    },
    error = function(e) {
        message("Error in Outgoing Pattern Analysis: ", e$message)
        message("  (If 'dim(X) must have a positive length', it means zero interactions exist in the matrix.)")
    }
)

# 3. Analyze Incoming Patterns (Receiver)
message("Analyzing Incoming (Receiver) Patterns...")
tryCatch(
    {
        cc <- CellChat::identifyCommunicationPatterns(cc, pattern = "incoming", k = opt$n_patterns)

        save_plot("patterns_incoming_river", function() {
            CellChat::netAnalysis_river(cc, pattern = "incoming")
        }, w = 12, h = 8)
        save_plot("patterns_incoming_dot", function() {
            CellChat::netAnalysis_dot(cc, pattern = "incoming")
        }, w = 10, h = 8)
    },
    error = function(e) {
        message("Error in Incoming Pattern Analysis: ", e$message)
    }
)

message("Done.")
