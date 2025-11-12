#' Data Conversion and Validation Utilities
#'
#' Standardized functions for converting Seurat objects to dataframes
#' and validating data integrity.
#'
#' @name utils_data
NULL

#' Convert Seurat Object to Dataframe
#'
#' Converts a Seurat object to a dataframe by combining expression data
#' and metadata. Supports flexible feature selection from both assays and metadata.
#'
#' @param sobj Seurat object
#' @param features Character vector of features to extract. Can include:
#'   - Gene names (from assay)
#'   - Metadata column names
#'   - If NULL, returns only metadata
#' @param assay Character. Assay to use (default: DefaultAssay(sobj))
#' @param layer Character. Layer to use ("counts", "data", "scale.data")
#' @param metadata_only Logical. If TRUE, only return metadata without expression
#' @param add_cell_id Logical. If TRUE, add cell IDs as a column (default: TRUE)
#' @param metadata_df Optional data.frame. Additional metadata to join.
#'   Must have rownames or a column matching cell IDs.
#' @param match_id Character. Column name in metadata_df to match with cell IDs.
#'   If NULL, uses rownames of metadata_df.
#'
#' @return data.frame with cells as rows and features + metadata as columns
#'
#' @examples
#' \dontrun{
#' # Basic conversion
#' df <- .sobj_to_df(sobj, features = c("CD3D", "CD4", "nFeature_RNA"))
#'
#' # With additional metadata
#' extra_meta <- data.frame(cell_id = colnames(sobj), new_var = rnorm(ncol(sobj)))
#' df <- .sobj_to_df(sobj, features = "CD3D", metadata_df = extra_meta, match_id = "cell_id")
#' }
#'
#' @keywords internal
.sobj_to_df <- function(sobj,
                        features = NULL,
                        assay = NULL,
                        layer = "data",
                        metadata_only = FALSE,
                        add_cell_id = TRUE,
                        metadata_df = NULL,
                        match_id = NULL) {
  
  if (!inherits(sobj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  `%||%` <- rlang::`%||%`
  
  # Get assay
  assay <- assay %||% Seurat::DefaultAssay(sobj)
  if (!assay %in% names(sobj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object")
  }
  
  # Start with metadata
  result_df <- as.data.frame(sobj@meta.data)
  
  # Add cell IDs as column if requested
  if (add_cell_id) {
    result_df$.cell_id <- rownames(result_df)
  }
  
  # Extract expression features if requested
  if (!metadata_only && !is.null(features)) {
    meta_cols <- names(sobj@meta.data)
    meta_features <- intersect(features, meta_cols)
    assay_features <- setdiff(features, meta_cols)
    
    # Get expression data for assay features
    if (length(assay_features) > 0) {
      # Check if features exist in assay
      available_features <- intersect(assay_features, rownames(sobj[[assay]]))
      missing_features <- setdiff(assay_features, available_features)
      
      if (length(missing_features) > 0) {
        warning("Features not found in assay: ", paste(missing_features, collapse = ", "))
      }
      
      if (length(available_features) > 0) {
        expr_data <- Seurat::GetAssayData(sobj, assay = assay, layer = layer)
        expr_df <- as.data.frame(t(as.matrix(expr_data[available_features, , drop = FALSE])))
        
        # Merge with result
        result_df <- cbind(result_df, expr_df)
      }
    }
    
    # Meta features are already in result_df
    if (length(meta_features) > 0) {
      # Verify they exist (already checked above)
      missing_meta <- setdiff(meta_features, names(result_df))
      if (length(missing_meta) > 0) {
        warning("Metadata features not found: ", paste(missing_meta, collapse = ", "))
      }
    }
  }
  
  # Join additional metadata if provided
  if (!is.null(metadata_df)) {
    if (!is.data.frame(metadata_df)) {
      stop("metadata_df must be a data.frame")
    }
    
    # Determine join key
    if (is.null(match_id)) {
      # Use rownames
      if (is.null(rownames(metadata_df))) {
        stop("metadata_df must have rownames or match_id must be specified")
      }
      join_key <- rownames(result_df)
      metadata_df$.join_key <- rownames(metadata_df)
      match_id <- ".join_key"
    } else {
      if (!match_id %in% names(metadata_df)) {
        stop("match_id '", match_id, "' not found in metadata_df")
      }
      join_key <- result_df$.cell_id %||% rownames(result_df)
    }
    
    # Perform join
    result_df$.join_key <- join_key
    result_df <- dplyr::left_join(
      result_df,
      metadata_df,
      by = setNames(match_id, ".join_key")
    )
    
    # Clean up temporary join key
    result_df$.join_key <- NULL
    if (match_id == ".join_key") {
      metadata_df$.join_key <- NULL
    }
  }
  
  return(result_df)
}

#' Validate Dataframe for Plotting
#'
#' Checks dataframe integrity and required columns for plotting functions.
#'
#' @param df data.frame to validate
#' @param required_cols Character vector of required column names
#' @param check_na Logical. Check for NA values in required columns (default: TRUE)
#' @param check_duplicates Logical. Check for duplicate rownames (default: FALSE)
#' @param min_rows Integer. Minimum number of rows required (default: 1)
#'
#' @return Invisibly returns TRUE if validation passes, stops with error otherwise
#'
#' @examples
#' \dontrun{
#' .validate_df(df, required_cols = c("CD3D", "cell_type"), check_na = TRUE)
#' }
#'
#' @keywords internal
.validate_df <- function(df,
                         required_cols = NULL,
                         check_na = TRUE,
                         check_duplicates = FALSE,
                         min_rows = 1) {
  
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame")
  }
  
  # Check minimum rows
  if (nrow(df) < min_rows) {
    stop("Dataframe must have at least ", min_rows, " rows, but has ", nrow(df))
  }
  
  # Check required columns
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Check for NA values
    if (check_na) {
      na_cols <- character()
      for (col in required_cols) {
        if (any(is.na(df[[col]]))) {
          na_cols <- c(na_cols, col)
        }
      }
      if (length(na_cols) > 0) {
        warning("NA values found in columns: ", paste(na_cols, collapse = ", "))
      }
    }
  }
  
  # Check for duplicate rownames
  if (check_duplicates) {
    if (any(duplicated(rownames(df)))) {
      warning("Duplicate rownames found in dataframe")
    }
  }
  
  return(invisible(TRUE))
}

#' Prepare Plot Data from Seurat or Dataframe
#'
#' Unified function to prepare data for plotting from either Seurat object
#' or dataframe input. Handles feature extraction, metadata joining, and validation.
#'
#' @param data Seurat object or data.frame
#' @param features Character vector of features to extract
#' @param group.by Character. Column name for grouping (facet)
#' @param split.by Character. Column name for splitting (x-axis or color)
#' @param sample_col Character. Column name for sample/patient ID (for aggregation)
#' @param assay Character. Assay to use (Seurat only)
#' @param layer Character. Layer to use (Seurat only)
#' @param metadata_df Optional data.frame. Additional metadata to join
#' @param match_id Character. Column name in metadata_df for matching
#'
#' @return data.frame ready for plotting
#'
#' @examples
#' \dontrun{
#' plot_df <- .prepare_plot_data(
#'   data = sobj,
#'   features = c("CD3D", "CD4"),
#'   group.by = "cell_type",
#'   split.by = "condition"
#' )
#' }
#'
#' @keywords internal
.prepare_plot_data <- function(data,
                               features = NULL,
                               group.by = NULL,
                               split.by = NULL,
                               sample_col = NULL,
                               assay = NULL,
                               layer = "data",
                               metadata_df = NULL,
                               match_id = NULL) {
  
  `%||%` <- rlang::`%||%`
  
  # Convert Seurat to dataframe if needed
  if (inherits(data, "Seurat")) {
    plot_df <- .sobj_to_df(
      sobj = data,
      features = features,
      assay = assay,
      layer = layer,
      metadata_df = metadata_df,
      match_id = match_id
    )
  } else if (is.data.frame(data)) {
    plot_df <- data
    
    # Check required columns
    req_cols <- unique(c(features, group.by, split.by, sample_col))
    req_cols <- req_cols[!is.null(req_cols)]
    
    missing_cols <- setdiff(req_cols, names(plot_df))
    if (length(missing_cols) > 0) {
      stop("Missing columns in dataframe: ", paste(missing_cols, collapse = ", "))
    }
    
    # Join additional metadata if provided
    if (!is.null(metadata_df)) {
      if (is.null(match_id)) {
        if (is.null(rownames(plot_df)) || is.null(rownames(metadata_df))) {
          stop("Both plot_df and metadata_df must have rownames, or match_id must be specified")
        }
        plot_df$.join_key <- rownames(plot_df)
        metadata_df$.join_key <- rownames(metadata_df)
        match_id <- ".join_key"
      }
      
      plot_df <- dplyr::left_join(
        plot_df,
        metadata_df,
        by = setNames(match_id, ".join_key")
      )
      
      plot_df$.join_key <- NULL
      if (match_id == ".join_key") {
        metadata_df$.join_key <- NULL
      }
    }
  } else {
    stop("data must be a Seurat object or data.frame")
  }
  
  # Validate
  req_cols <- unique(c(features, group.by, split.by, sample_col))
  req_cols <- req_cols[!is.null(req_cols)]
  if (length(req_cols) > 0) {
    .validate_df(plot_df, required_cols = req_cols, check_na = FALSE)
  }
  
  return(plot_df)
}

