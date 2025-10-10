#' Gene Signature Scoring Functions
#'
#' This module provides functions for calculating and visualizing gene signature
#' scores using various methods including AddModuleScore and enrichIt.
#'
#' @name signature_scoring
NULL

#' Add Module Scores for Multiple Gene Sets
#'
#' Calculates module scores for multiple gene sets using Seurat's AddModuleScore.
#' Automatically handles naming and removes the default "1" suffix.
#'
#' @param seurat_object A Seurat object
#' @param feature_sets Named or unnamed list of character vectors (gene sets)
#' @param assay Name of the assay to use (default: DefaultAssay)
#' @param slot Slot to pull expression data from (default: "data")
#' @param nbin Number of bins for AddModuleScore (default: 24)
#' @param ctrl Number of control features (default: 100)
#' @param seed Random seed (default: 1)
#' @param search Passed to Seurat::[.Assay (default: FALSE)
#' @param ... Additional arguments passed to Seurat::AddModuleScore
#'
#' @return Seurat object with added module scores in metadata
#'
#' @examples
#' \dontrun{
#' gene_sets <- list(
#'   Tcell_activation = c("CD69", "IFNG", "TNF"),
#'   Monocyte_markers = c("CD14", "LYZ")
#' )
#' pbmc <- AddMultipleModuleScores(pbmc, gene_sets)
#' }
#'
#' @export
AddMultipleModuleScores <- function(seurat_object,
                                    feature_sets,
                                    assay = NULL,
                                    slot = "data",
                                    nbin = 24,
                                    ctrl = 100,
                                    seed = 1,
                                    search = FALSE,
                                    ...) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  
  if (!is.list(feature_sets)) {
    stop("'feature_sets' must be a list of character vectors.")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_object)
  }
  
  all_genes <- rownames(Seurat::GetAssayData(seurat_object, assay = assay, slot = slot))
  added_cols <- character(0)
  
  for (i in seq_along(feature_sets)) {
    genes_raw <- feature_sets[[i]]
    genes_use <- intersect(genes_raw, all_genes)
    
    if (length(genes_use) == 0) {
      warning(sprintf("feature set %d: no genes found – skipped", i))
      next
    }
    
    set_name <- names(feature_sets)[i]
    if (is.null(set_name) || set_name == "") {
      set_name <- paste(gsub("-", "_", genes_use), collapse = "+")
    }
    
    # Run AddModuleScore
    before_cols <- colnames(seurat_object[[]])
    seurat_object <- Seurat::AddModuleScore(
      object = seurat_object,
      features = list(genes_use),
      name = set_name,
      assay = assay,
      slot = slot,
      nbin = nbin,
      ctrl = ctrl,
      seed = seed,
      search = search,
      ...
    )
    after_cols <- colnames(seurat_object[[]])
    new_col <- setdiff(after_cols, before_cols)
    
    # Remove trailing "1" and rename
    tidy_col <- sub("1$", "", new_col)
    if (tidy_col != new_col) {
      colnames(seurat_object[[]])[match(new_col, after_cols)] <- tidy_col
      message(sprintf("renamed '%s' -> '%s'", new_col, tidy_col))
    }
    
    added_cols <- c(added_cols, tidy_col)
  }
  
  # Print usage example
  if (length(added_cols)) {
    msg1 <- "# Copy-&-paste for FeaturePlot:"
    msg2 <- sprintf("FeaturePlot(obj, features = c(%s))",
                    paste(sprintf("'%s'", added_cols), collapse = ", "))
    message("\n", msg1, "\n", msg2, "\n")
    flush.console()
  } else {
    warning("No module scores were added.")
  }
  
  invisible(seurat_object)
}

#' Add Gene Signature Score using enrichIt
#'
#' Calculates a gene signature score using escape::enrichIt and adds it to 
#' the Seurat object's metadata. Flexibly accepts different gene ID types.
#'
#' @param seurat_obj A Seurat object
#' @param gene_source File path, data.frame, or vector containing gene IDs
#' @param signature_name Name for the new metadata column
#' @param input_keytype Type of input gene IDs (default: "ENSEMBL")
#' @param gene_col Column index/name if gene_source is file or data.frame (default: 1)
#' @param sheet_name Sheet name/index if xlsx file (default: 1)
#' @param assay Assay to use (default: "RNA")
#' @param layer Layer (slot) to use (default: "data")
#' @param ... Additional arguments passed to escape::enrichIt
#'
#' @return Seurat object with signature score added to metadata
#'
#' @export
add_signature_enrichit <- function(seurat_obj,
                                   gene_source,
                                   signature_name,
                                   input_keytype = "ENSEMBL",
                                   gene_col = 1,
                                   sheet_name = 1,
                                   assay = "RNA",
                                   layer = "data",
                                   ...) {
  
  if (!requireNamespace("escape", quietly = TRUE)) {
    stop("escape package required. Install with: BiocManager::install('escape')")
  }
  
  if (!requireNamespace("GSEABase", quietly = TRUE)) {
    stop("GSEABase package required. Install with: BiocManager::install('GSEABase')")
  }
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package required. Install with: BiocManager::install('org.Hs.eg.db')")
  }
  
  # Load gene list
  if (is.character(gene_source) && file.exists(gene_source)) {
    ext <- tools::file_ext(gene_source)
    gene_list_raw <- switch(
      ext,
      xlsx = readxl::read_xlsx(gene_source, sheet = sheet_name)[[gene_col]],
      csv = read.csv(gene_source, stringsAsFactors = FALSE)[[gene_col]],
      txt = read.table(gene_source, stringsAsFactors = FALSE)[[gene_col]],
      stop("Unsupported file type.")
    )
  } else if (is.data.frame(gene_source)) {
    gene_list_raw <- gene_source[[gene_col]]
  } else if (is.vector(gene_source)) {
    gene_list_raw <- gene_source
  } else {
    stop("`gene_source` must be a valid file path, data.frame, or vector.")
  }
  
  # Convert to SYMBOL
  message(paste("Input keytype is", input_keytype, ". Converting to Gene Symbols..."))
  gene_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(na.omit(as.character(gene_list_raw))),
    keytype = input_keytype,
    column = "SYMBOL",
    multiVals = "first"
  )
  gene_symbols <- na.omit(gene_symbols)
  
  if (length(gene_symbols) == 0) {
    stop(paste("No valid gene symbols could be mapped using keytype:", input_keytype))
  }
  message(paste(length(gene_symbols), "gene symbols were successfully mapped."))
  
  # Create GeneSet
  gene_set <- GSEABase::GeneSet(gene_symbols, setName = signature_name)
  gene_sets_collection <- GSEABase::GeneSetCollection(gene_set)
  
  # Run enrichIt
  message("Running enrichIt...")
  expr_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  enrichment_scores <- escape::enrichIt(
    obj = expr_matrix,
    gene.sets = gene_sets_collection,
    ...
  )
  
  # Add to metadata
  message("Adding scores to Seurat metadata...")
  seurat_obj <- Seurat::AddMetaData(
    object = seurat_obj,
    metadata = as.data.frame(enrichment_scores)
  )
  
  return(seurat_obj)
}

#' Add Pathway Activity Scores using progeny
#'
#' Infers pathway activities using progeny and adds them to the Seurat 
#' object's metadata.
#'
#' @param seurat_obj A Seurat object
#' @param organism Organism: "Human" or "Mouse" (default: "Human")
#' @param topn Number of top genes per pathway (default: 100)
#' @param ... Additional arguments passed to progeny::progeny
#'
#' @return Seurat object with pathway activity scores in metadata
#'
#' @export
add_progeny_scores <- function(seurat_obj, organism = "Human", topn = 100, ...) {
  
  if (!requireNamespace("progeny", quietly = TRUE)) {
    stop("progeny package required. Install with: BiocManager::install('progeny')")
  }
  
  # Run progeny
  message("Running progeny...")
  seurat_obj <- progeny::progeny(
    seurat_obj,
    scale = FALSE,
    organism = organism,
    topn = topn,
    return_assay = TRUE,
    ...
  )
  
  # Scale progeny assay
  message("Scaling progeny assay...")
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = "progeny")
  
  # Add to metadata
  message("Adding scores to Seurat metadata...")
  progeny_scores <- as.data.frame(
    t(Seurat::GetAssayData(seurat_obj, assay = "progeny", slot = "scale.data"))
  )
  
  seurat_obj <- Seurat::AddMetaData(
    object = seurat_obj,
    metadata = progeny_scores
  )
  
  return(seurat_obj)
}

#' Score New Data with Existing Signature
#'
#' Applies a gene signature (from find_gene_signature) to score new expression data.
#'
#' @param expr_data Seurat object or expression matrix
#' @param signature gene_signature object from find_gene_signature
#' @param normalize Whether to z-score normalize scores (default: TRUE)
#'
#' @return Named numeric vector of signature scores
#'
#' @export
score_signature <- function(expr_data, signature, normalize = TRUE) {
  genes <- signature$genes
  weights <- signature$weights
  
  # Extract expression matrix
  if (inherits(expr_data, "Seurat")) {
    expr_mat <- as.matrix(Seurat::GetAssayData(expr_data, slot = "data"))
  } else {
    expr_mat <- as.matrix(expr_data)
  }
  
  # Check gene availability
  available_genes <- intersect(genes, rownames(expr_mat))
  if (length(available_genes) == 0) {
    stop("None of the signature genes found in data")
  }
  if (length(available_genes) < length(genes)) {
    warning(sprintf("%d/%d signature genes not found in data",
                    length(genes) - length(available_genes), length(genes)))
  }
  
  # Calculate scores
  weights <- weights[available_genes]
  scores <- colSums(expr_mat[available_genes, , drop = FALSE] * weights)
  
  if (normalize) {
    scores <- scale(scores)[, 1]
  }
  
  return(scores)
}

