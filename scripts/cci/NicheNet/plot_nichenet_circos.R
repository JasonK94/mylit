#!/usr/bin/env Rscript
# Plot NicheNet Circos
# Usage: Rscript plot_nichenet_circos.R -f <results.qs>

suppressPackageStartupMessages({
    library(optparse)
    library(circlize)
    library(qs)
})

option_list <- list(
    make_option(c("-f", "--file"),
        type = "character", default = NULL,
        help = "Path to NicheNet results file (.qs)"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "nichenet_circos.pdf",
        help = "Output PDF file path"
    ),
    make_option(c("--format"),
        type = "character", default = "pdf",
        help = "Output format (pdf/png)"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)) {
    print_help(opt_parser)
    stop("Input file is required.")
}

source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_cci.R")

cat("Generating NicheNet Circos Plot...\n")
cat("  File: ", opt$file, "\n")
cat("  Output: ", opt$output, "\n")

replot_nichenet_circos(
    results_file = opt$file,
    output_file = opt$output,
    format = opt$format,
    verbose = TRUE
)

cat("\n✅ Plot saved!\n")
