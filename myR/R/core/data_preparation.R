#' Core Data Preparation Functions
#'
#' This module provides core functions for preparing data from Seurat objects
#' for various downstream analyses.
#'
#' @name data_preparation
NULL

#' Extract Feature from Seurat Object
#'
#' Internal helper to extract a feature (gene or metadata column) from a Seurat object.
#'
#' @param object Seurat object
#' @param feature Feature name (gene or metadata column)
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use (default: "data")
#' @param cells Subset of cells to extract (default: all cells)
#'
#' @return Numeric vector of feature values
#'
#' @keywords internal
.get_feature_vector <- function(object, 
                                feature, 
                                assay = NULL, 
                                slot = "data",
                                cells = NULL) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  if (is.null(cells)) {
    cells <- colnames(object)
  }
  
  # Check if feature is in metadata
  if (feature %in% colnames(object@meta.data)) {
    values <- object@meta.data[cells, feature]
    return(as.numeric(values))
  }
  
  # Check if feature is a gene
  if (feature %in% rownames(object[[assay]])) {
    values <- Seurat::GetAssayData(object, assay = assay, slot = slot)[feature, cells]
    return(as.numeric(values))
  }
  
  stop("Feature '", feature, "' not found in metadata or assay '", assay, "'")
}

#' Extract Feature Vector (User-Facing)
#'
#' Extracts a feature (gene or metadata column) from a Seurat object.
#'
#' @param object Seurat object
#' @param feature Feature name (gene or metadata column)
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use (default: "data")
#' @param cells Subset of cells to extract (default: all cells)
#'
#' @return Numeric vector of feature values
#'
#' @examples
#' \dontrun{
#' get_feature_vec(seurat_obj, "CD3D")
#' get_feature_vec(seurat_obj, "nCount_RNA")
#' }
#'
#' @export
get_feature_vec <- function(object, feature, assay = NULL, slot = "data", cells = NULL) {
  .get_feature_vector(object, feature, assay, slot, cells)
}

#' Extract Multiple Features
#'
#' Extracts multiple features from a Seurat object and returns as a matrix.
#'
#' @param object Seurat object
#' @param features Character vector of feature names
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use (default: "data")
#' @param cells Subset of cells to extract (default: all cells)
#'
#' @return Matrix with features as rows and cells as columns
#'
#' @export
get_feature_matrix <- function(object, features, assay = NULL, slot = "data", cells = NULL) {
  
  if (is.null(cells)) {
    cells <- colnames(object)
  }
  
  feature_data <- sapply(features, function(feat) {
    tryCatch(
      .get_feature_vector(object, feat, assay, slot, cells),
      error = function(e) {
        warning("Failed to extract feature '", feat, "': ", e$message)
        rep(NA_real_, length(cells))
      }
    )
  })
  
  # Ensure matrix format
  if (!is.matrix(feature_data)) {
    feature_data <- as.matrix(feature_data)
  }
  
  rownames(feature_data) <- cells
  colnames(feature_data) <- features
  
  return(t(feature_data))
}

#' Prepare Metadata Table
#'
#' Extracts and optionally filters metadata from a Seurat object.
#'
#' @param object Seurat object
#' @param columns Specific columns to extract (NULL = all columns)
#' @param cells Subset of cells (NULL = all cells)
#' @param drop_na Whether to drop rows with NA values (default: FALSE)
#'
#' @return Data frame of metadata
#'
#' @export
prepare_metadata_table <- function(object, columns = NULL, cells = NULL, drop_na = FALSE) {
  
  if (is.null(cells)) {
    cells <- colnames(object)
  }
  
  metadata <- object@meta.data[cells, , drop = FALSE]
  
  if (!is.null(columns)) {
    missing_cols <- setdiff(columns, colnames(metadata))
    if (length(missing_cols) > 0) {
      stop("Columns not found in metadata: ", paste(missing_cols, collapse = ", "))
    }
    metadata <- metadata[, columns, drop = FALSE]
  }
  
  if (drop_na) {
    metadata <- metadata[complete.cases(metadata), , drop = FALSE]
  }
  
  return(metadata)
}

#' Aggregate Expression by Group
#'
#' Aggregates gene expression across cells within groups (e.g., clusters, samples).
#'
#' @param object Seurat object
#' @param features Features (genes) to aggregate
#' @param group_by Metadata column to group by
#' @param method Aggregation method: "mean", "median", "sum" (default: "mean")
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use (default: "data")
#'
#' @return Matrix with features as rows and groups as columns
#'
#' @export
aggregate_expression_by_group <- function(object, 
                                          features, 
                                          group_by, 
                                          method = "mean",
                                          assay = NULL,
                                          slot = "data") {
  
  # Validate inputs
  if (!group_by %in% colnames(object@meta.data)) {
    stop("group_by '", group_by, "' not found in metadata")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Get grouping variable
  groups <- object@meta.data[[group_by]]
  unique_groups <- unique(groups)
  
  # Initialize result matrix
  result <- matrix(
    NA_real_, 
    nrow = length(features), 
    ncol = length(unique_groups),
    dimnames = list(features, as.character(unique_groups))
  )
  
  # Aggregate for each group
  for (grp in unique_groups) {
    cells_in_group <- colnames(object)[groups == grp]
    
    if (length(cells_in_group) == 0) next
    
    for (feat in features) {
      tryCatch({
        values <- .get_feature_vector(object, feat, assay, slot, cells_in_group)
        
        result[feat, as.character(grp)] <- switch(
          method,
          "mean" = mean(values, na.rm = TRUE),
          "median" = median(values, na.rm = TRUE),
          "sum" = sum(values, na.rm = TRUE),
          stop("Unknown method: ", method)
        )
      }, error = function(e) {
        warning("Failed to aggregate feature '", feat, "' for group '", grp, "'")
      })
    }
  }
  
  return(result)
}

#' Prepare Count Matrix
#'
#' Extracts a count matrix from a Seurat object with optional filtering.
#'
#' @param object Seurat object
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#' @param features Features to include (NULL = all)
#' @param cells Cells to include (NULL = all)
#' @param min_cells Minimum number of cells expressing a feature (default: 0)
#' @param min_features Minimum number of features expressed in a cell (default: 0)
#'
#' @return Sparse or dense matrix of counts
#'
#' @export
prepare_count_matrix <- function(object, 
                                 assay = "RNA",
                                 slot = "counts",
                                 features = NULL,
                                 cells = NULL,
                                 min_cells = 0,
                                 min_features = 0) {
  
  # Extract counts
  counts <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  
  # Subset to specified features/cells
  if (!is.null(features)) {
    features <- intersect(features, rownames(counts))
    counts <- counts[features, , drop = FALSE]
  }
  
  if (!is.null(cells)) {
    cells <- intersect(cells, colnames(counts))
    counts <- counts[, cells, drop = FALSE]
  }
  
  # Filter by expression
  if (min_cells > 0) {
    n_cells_per_feature <- Matrix::rowSums(counts > 0)
    keep_features <- n_cells_per_feature >= min_cells
    counts <- counts[keep_features, , drop = FALSE]
  }
  
  if (min_features > 0) {
    n_features_per_cell <- Matrix::colSums(counts > 0)
    keep_cells <- n_features_per_cell >= min_features
    counts <- counts[, keep_cells, drop = FALSE]
  }
  
  return(counts)
}

#' Convert Expression to Long Format
#'
#' Converts expression data to long (tidy) format for ggplot2.
#'
#' @param object Seurat object
#' @param features Features to include
#' @param metadata_cols Metadata columns to include
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use (default: "data")
#' @param cells Cells to include (NULL = all)
#'
#' @return Data frame in long format with columns: cell, feature, expression, and metadata
#'
#' @export
convert_to_long_format <- function(object, 
                                   features, 
                                   metadata_cols = NULL,
                                   assay = NULL,
                                   slot = "data",
                                   cells = NULL) {
  
  if (is.null(cells)) {
    cells <- colnames(object)
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Get expression matrix
  expr_mat <- get_feature_matrix(object, features, assay, slot, cells)
  
  # Convert to long format
  expr_long <- as.data.frame(as.table(expr_mat))
  colnames(expr_long) <- c("feature", "cell", "expression")
  
  # Add metadata
  if (!is.null(metadata_cols)) {
    metadata <- prepare_metadata_table(object, metadata_cols, cells)
    expr_long <- merge(expr_long, metadata, by.x = "cell", by.y = "row.names")
  }
  
  return(expr_long)
}

#' Check Data Quality
#'
#' Performs basic quality checks on expression data.
#'
#' @param data Matrix or data frame of expression data
#' @param check_finite Check for infinite values (default: TRUE)
#' @param check_na Check for NA values (default: TRUE)
#' @param check_negative Check for negative values (default: TRUE)
#'
#' @return List with logical flags and counts of issues found
#'
#' @export
check_data_quality <- function(data, 
                               check_finite = TRUE,
                               check_na = TRUE,
                               check_negative = TRUE) {
  
  issues <- list()
  
  if (check_finite) {
    n_infinite <- sum(!is.finite(as.matrix(data)))
    issues$has_infinite <- n_infinite > 0
    issues$n_infinite <- n_infinite
  }
  
  if (check_na) {
    n_na <- sum(is.na(data))
    issues$has_na <- n_na > 0
    issues$n_na <- n_na
  }
  
  if (check_negative) {
    n_negative <- sum(as.matrix(data) < 0, na.rm = TRUE)
    issues$has_negative <- n_negative > 0
    issues$n_negative <- n_negative
  }
  
  issues$n_total <- length(as.matrix(data))
  issues$is_clean <- !any(unlist(issues[grepl("^has_", names(issues))]))
  
  return(issues)
}

