#!/usr/bin/env Rscript

# scripts/cellchat/inspect_cellchat.R
# Tool to inspect internal CellChat data structures

suppressPackageStartupMessages({
    library(optparse)
    library(CellChat)
    library(qs)
})

option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to CellChat object (.qs)"),
    make_option(c("-n", "--top_n"), type = "integer", default = 20, help = "Number of top items to show")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("Input required (-i)")

message("Loading: ", opt$input)
cc <- qs::qread(opt$input)

message("\n=== 1. Basic Info ===")
print(cc)

if ("ids" %in% slotNames(cc)) {
    message("  Samples/Datasets: ", length(names(cc@ids)))
    if (length(names(cc@ids)) > 0) message("  Names: ", paste(head(names(cc@ids)), collapse = ", "), "...")
} else {
    message("  (Single dataset object or 'ids' slot not present)")
}

if ("idents" %in% slotNames(cc)) {
    message("  Cell Groups: ", length(levels(cc@idents)))
    message("  Levels: ", paste(head(levels(cc@idents)), collapse = ", "), "...")
}

message("\n=== 2. Interaction Data (@net) ===")
if (!is.null(cc@net$prob)) {
    prob <- cc@net$prob
    message("  Prob Array Dim: ", paste(dim(prob), collapse = " x "))

    # Flatten non-zero probs
    vals <- prob[prob > 0]
    message("  Non-zero Interactions: ", length(vals))
    message("  Probability Stats:")
    print(summary(vals))

    message("  Top ", opt$top_n, " Strongest Interactions:")
    # Find indices of top values
    idx <- order(vals, decreasing = TRUE)[1:min(length(vals), opt$top_n)]
    top_vals <- vals[idx]
    # We need to map back to names (complex for 3D array) -> Just show values
    print(head(top_vals, opt$top_n))
} else {
    message("  @net$prob is NULL (Aggregation might be needed)")
}

if (!is.null(cc@net$count)) {
    message("\n  Aggregated Count Matrix Dim: ", paste(dim(cc@net$count), collapse = " x "))
    message("  Total Count Sum: ", sum(cc@net$count))
}

message("\n=== 3. Ligand-Receptor Database (@LR) ===")
if (!is.null(cc@LR$LRsig)) {
    lrsig <- cc@LR$LRsig
    message("  Significant LR pairs used: ", nrow(lrsig))
    message("  Columns: ", paste(colnames(lrsig), collapse = ", "))
    print(head(lrsig, opt$top_n))
}

message("\n=== 4. Pathways (@netP) ===")
if (!is.null(cc@netP$pathways)) {
    pathways <- cc@netP$pathways
    message("  Significant Pathways: ", length(pathways))
    print(head(pathways, opt$top_n))

    if (!is.null(cc@netP$prob)) {
        message("  Pathway Prob Stats (Non-zero):")
        print(summary(cc@netP$prob[cc@netP$prob > 0]))
    }
}

message("\n=== Done ===")
