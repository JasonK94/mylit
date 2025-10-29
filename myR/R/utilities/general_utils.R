#' General Utility Functions
#'
#' This module provides general-purpose utility functions used across the package.
#'
#' @name general_utils
NULL

#' Downsample a Seurat object by a specified ratio
#'
#' This function randomly samples cells from a Seurat object to create a smaller subset.
#' The sampling is done without replacement and maintains the original data structure.
#'
#' @param sobj A Seurat object to be downsampled
#' @param ratio Integer indicating the downsampling ratio (1:ratio). Default is 10, meaning
#'              the output will contain 1/10th of the original cells
#' @param seed Integer specifying the random seed for reproducibility. Default is 1234
#'
#' @return A downsampled Seurat object containing a subset of cells from the input object
#'
#' @examples
#' \dontrun{
#' # Downsample a Seurat object to 1/10th of its original size
#' sobj_down <- downsample_sobj(sobj, ratio = 10)
#' 
#' # Downsample with a different ratio and seed
#' sobj_down <- downsample_sobj(sobj, ratio = 5, seed = 42)
#' }
#' @export
downsample_sobj <- function(sobj, ratio = 10, seed = 1234) {
  set.seed(seed)
  barcodes <- colnames(sobj)
  barcodes_down <- sample(barcodes, length(barcodes) %/% ratio)
  sobj_down <- sobj[, barcodes_down]
  return(sobj_down)
}

#' Print Marker Genes in Comma-Separated Format
#'
#' Convenience function to print marker genes in a format ready for copy-paste.
#'
#' @param markers Data frame with marker gene results
#' @param sign Sign of log2FC to print ("+" for positive, "-" for negative)
#' @param num Number of genes to print (default: 100)
#' @param pseudobulk Whether markers are from pseudobulk analysis (uses "logFC" instead of "avg_log2FC")
#'
#' @return Prints genes to console
#' @export
printmy <- function(markers, sign = "+", num = 100, pseudobulk = FALSE) {
  if (!"gene" %in% names(markers)) {
    markers$gene <- rownames(markers)
  }
  
  if (!pseudobulk) {
    if (sign == "-") {
      print(paste(markers[markers$avg_log2FC < 0, ][1:num, ]$gene, collapse = ", "))
    } else {
      print(paste(markers[markers$avg_log2FC > 0, ][1:num, ]$gene, collapse = ", "))
    } 
  } else {
    if (sign == "-") {
      print(paste(markers[markers$logFC < 0, ][1:num, ]$gene, collapse = ", "))
    } else {
      print(paste(markers[markers$logFC > 0, ][1:num, ]$gene, collapse = ", "))
    } 
  }
}

#' Print Multiple Marker Lists
#'
#' Print marker genes from a named list of marker data frames.
#'
#' @param markers_list Named list of marker data frames
#' @param ... Additional arguments passed to printmy()
#'
#' @return Prints genes to console with list names as headers
#' @export
printMy <- function(markers_list, ...) {
  args <- list(...)
  for (name in names(markers_list)) {
    print(name)
    printmy(markers_list[[name]], pseudobulk = args[["pseudobulk"]])
  }
}

#' Print Top Genes per Set Combination
#'
#' For a list of gene sets, print genes that are unique to each combination
#' (e.g., genes in A only, genes in A & B only, etc.)
#'
#' @param gene_list Named list of character vectors. Each element is a gene vector.
#' @param num_print Integer. Maximum number of genes to print per combination (default: 100).
#'
#' @examples
#' \dontrun{
#' L <- list(
#'   A = c("TP53","EGFR","MYC","BRCA1"),
#'   B = c("EGFR","MYC","PTEN","BRCA2"),
#'   C = c("MYC","PTEN","ALK")
#' )
#' print_gene_combinations(L, num_print = 2)
#' }
#' @export
print_gene_combinations <- function(gene_list, num_print = 100) {
  if (!is.list(gene_list) || is.null(names(gene_list))) {
    stop("gene_list must be a named list")
  }
  set_names <- names(gene_list)
  n_sets <- length(set_names)
  
  cat(sprintf("Printing top %d genes per combination...\n\n", num_print))
  
  # For each k-combination
  for (k in seq_len(n_sets)) {
    combos <- combn(set_names, k, simplify = FALSE)
    for (comb in combos) {
      # Combination name
      comb_name <- paste(comb, collapse = "&")
      
      # Genes common to all sets in this combination
      genes_in <- Reduce(intersect, gene_list[comb])
      
      # Exclude genes found in other sets
      other_sets <- setdiff(set_names, comb)
      if (length(other_sets) > 0) {
        genes_out <- unique(unlist(gene_list[other_sets], use.names = FALSE))
        genes_in <- setdiff(genes_in, genes_out)
      }
      
      # Sort and sample
      genes_in <- sort(genes_in)
      to_print <- head(genes_in, num_print)
      
      # Print
      cat(comb_name, ":\n")
      if (length(to_print) == 0) {
        cat("  (none)\n\n")
      } else {
        cat(" ", paste(to_print, collapse = ", "), "\n\n")
      }
    }
  }
}

#' NULL Coalescing Operator
#'
#' Returns the first non-NULL value.
#'
#' @param a First value
#' @param b Second value (returned if a is NULL)
#'
#' @return a if not NULL, otherwise b
#' @keywords internal
#' @export
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


