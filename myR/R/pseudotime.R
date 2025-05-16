#' Run Slingshot trajectory inference from a Seurat object
#'
#' This function takes a Seurat object (typically after PCA, UMAP, and clustering),
#' runs Slingshot to infer trajectories, and returns a SingleCellExperiment object
#' containing the Slingshot results (pseudotime and cell weights). Counts and
#' other relevant metadata are also included.
#'
#' @param seurat_obj A Seurat object. Must have dimensionality reduction (e.g., 'pca', 'umap')
#'   and cluster identities.
#' @param cluster_col A character string specifying the column name in Seurat object's
#'   metadata that contains the cluster labels to be used by Slingshot.
#' @param reduced_dim_name A character string specifying the name of the dimensionality
#'   reduction to use from the Seurat object (e.g., "UMAP", "PCA"). Default is "UMAP".
#' @param start_cluster A character string or numeric value specifying the identity of the
#'   starting cluster for trajectory inference. It will be coerced to character.
#'   If NULL, Slingshot will try to identify it automatically.
#' @param end_clusters Optional. A character vector or numeric vector specifying identities
#'   of known end clusters. These will be coerced to character.
#' @param counts_assay_name Character string, the name of the assay in the Seurat object
#'   containing counts suitable for downstream DE analysis (e.g., "RNA" or "SCT"). Default is "RNA".
#' @param main_trajectory_only Logical. Currently a placeholder; future versions might implement
#'   logic to select a primary trajectory if multiple disconnected ones are found. Default is FALSE.
#' @param ... Additional arguments to pass to `slingshot::slingshot()` (e.g., `omega`, `smoother`, `approx_points`).
#'
#' @return A SingleCellExperiment object with Slingshot's pseudotime and cell weights
#'   stored. Counts are also included. Returns NULL if a critical error occurs.
#'
#' @import Seurat
#' @import slingshot
#' @import SingleCellExperiment
#' @import SummarizedExperiment # For colData<-
#' @importFrom S4Vectors DataFrame # Explicitly import DataFrame
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'seu' is a Seurat object processed with PCA, UMAP, FindClusters
#' # sce_with_sling <- tryCatch({
#' #   run_slingshot_from_seurat(
#' #     seurat_obj = seu,
#' #     cluster_col = "seurat_clusters",
#' #     reduced_dim_name = "umap",
#' #     start_cluster = "0",
#' #     counts_assay_name = "RNA"
#' #   )
#' # }, error = function(e) {
#' #   message("Slingshot run failed: ", e$message)
#' #   return(NULL)
#' # })
#' #
#' # if (!is.null(sce_with_sling)) {
#' #   print("Slingshot analysis successful. SCE object created.")
#' # }
#' }
run_slingshot_from_seurat <- function(seurat_obj,
                                      cluster_col,
                                      reduced_dim_name = "UMAP",
                                      start_cluster = NULL,
                                      end_clusters = NULL,
                                      counts_assay_name = "RNA",
                                      main_trajectory_only = FALSE,
                                      ...) {
  # --- 0. Load required packages (or ensure they are loaded) ---
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("Package 'slingshot' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) { # For DataFrame
    stop("Package 'S4Vectors' is required. Please install it.", call. = FALSE)
  }
  
  
  # --- 1. Input Validation & Seurat Object Integrity ---
  message("--- Step 0: Validating inputs and Seurat object integrity ---")
  if (!is(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  base_cell_ids <- colnames(seurat_obj)
  if (is.null(base_cell_ids) || length(base_cell_ids) == 0) {
    stop("colnames(seurat_obj) are NULL or empty. Seurat object seems malformed.", call. = FALSE)
  }
  if (is.null(rownames(seurat_obj@meta.data))) {
    message("Warning: rownames(seurat_obj@meta.data) is NULL. Attempting to set from colnames(seurat_obj).")
    if (length(base_cell_ids) == nrow(seurat_obj@meta.data)) {
      rownames(seurat_obj@meta.data) <- base_cell_ids
      message("Successfully set rownames(seurat_obj@meta.data).")
    } else {
      stop("Cannot fix NULL rownames for meta.data: length mismatch with colnames(seurat_obj).", call. = FALSE)
    }
  }
  # Ensure meta.data rownames are perfectly aligned with Seurat object colnames
  if (!identical(base_cell_ids, rownames(seurat_obj@meta.data))) {
    message("Aligning Seurat object metadata rownames with colnames by reordering/subsetting meta.data.")
    # Check if all base_cell_ids are in meta.data rownames
    if(!all(base_cell_ids %in% rownames(seurat_obj@meta.data))){
      missing_in_meta <- base_cell_ids[!base_cell_ids %in% rownames(seurat_obj@meta.data)]
      warning("Some cells in colnames(seurat_obj) are missing from rownames(seurat_obj@meta.data): ",
              length(missing_in_meta), " cells. This should not happen. Object may be inconsistent.", call. = FALSE)
      # Forcibly subset seurat_obj to only cells present in meta.data if meta.data is smaller but has rownames
      # Or, more commonly, meta.data might have extra rows. We only care about cells in base_cell_ids.
    }
    seurat_obj@meta.data <- seurat_obj@meta.data[base_cell_ids, , drop = FALSE] # Reorder/subset meta.data
    if(!identical(colnames(seurat_obj), rownames(seurat_obj@meta.data))){
      stop("Failed to align cell names between Seurat object (colnames) and its meta.data (rownames) after attempting fix.", call. = FALSE)
    }
  }
  
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cluster column '", cluster_col, "' not found in Seurat object metadata. Available columns: ", paste(colnames(seurat_obj@meta.data), collapse=", ")), call. = FALSE)
  }
  
  available_reductions <- names(seurat_obj@reductions)
  if (!reduced_dim_name %in% available_reductions) {
    matched_rd <- available_reductions[tolower(available_reductions) == tolower(reduced_dim_name)]
    if(length(matched_rd) == 1){
      message("Note: Reduction '", reduced_dim_name, "' not found directly. Using case-insensitive match: '", matched_rd, "'.")
      reduced_dim_name <- matched_rd
    } else {
      stop(paste("Reduction '", reduced_dim_name, "' not found in Seurat object. Available reductions: ", paste(available_reductions, collapse=", ")), call. = FALSE)
    }
  }
  
  # Assay name validation
  message("--- Debug: Validating counts_assay_name (formal argument value) ---")
  # This 'counts_assay_name' is the formal argument, which should receive the user's value (e.g., "SCT")
  message("Formal argument 'counts_assay_name' inside function is: '", counts_assay_name, "' (class: ", class(counts_assay_name), ")")
  seurat_assay_names <- NULL
  tryCatch({
    seurat_assay_names <- Seurat::Assays(seurat_obj)
  }, error = function(e) {
    message("ERROR occurred when calling Seurat::Assays(seurat_obj): ", e$message)
    stop("Could not retrieve assay names using Seurat::Assays().", call. = FALSE)
  })
  if (is.null(seurat_assay_names) || !is.vector(seurat_assay_names) || !is.character(seurat_assay_names)) {
    stop(paste0("Seurat::Assays(seurat_obj) did NOT return a character vector. Returned class: '",
                paste(class(seurat_assay_names), collapse=", "), "'."), call. = FALSE)
  }
  if (!counts_assay_name %in% seurat_assay_names) {
    stop(paste("Assay '", counts_assay_name, "' (from formal argument) not found in Seurat object. Available: ", paste(seurat_assay_names, collapse=", ")), call. = FALSE)
  }
  message("Assay '", counts_assay_name, "' (from formal argument) successfully found in Seurat object.")
  
  
  message("--- Step 1: Extracting and validating core data components ---")
  dim_reduced_data <- Embeddings(seurat_obj, reduction = reduced_dim_name)
  cluster_labels_vector <- seurat_obj@meta.data[[cluster_col]]
  # Use the validated counts_assay_name (formal argument) for extracting counts
  counts_data <- GetAssayData(seurat_obj, assay = counts_assay_name, slot = "counts")
  
  if (is.null(dim_reduced_data)) stop(paste("Embeddings for reduction '", reduced_dim_name, "' returned NULL."), call. = FALSE)
  if (!is.matrix(dim_reduced_data)) {
    dim_reduced_data <- as.matrix(dim_reduced_data)
    if(!is.matrix(dim_reduced_data)) stop("Failed to coerce embeddings to matrix.")
  }
  if(is.null(rownames(dim_reduced_data))) {
    message("CRITICAL: rownames(Embeddings for '", reduced_dim_name, "') is NULL. Attempting to fix using Seurat object cell names.")
    if(nrow(dim_reduced_data) == length(base_cell_ids)) { # base_cell_ids from colnames(seurat_obj)
      rownames(dim_reduced_data) <- base_cell_ids
      if(is.null(rownames(dim_reduced_data))) stop("Failed to assign rownames to dim_reduced_data.")
      message("Successfully assigned rownames to dim_reduced_data.")
    } else {
      stop(paste("Cannot assign rownames to dim_reduced_data: nrow mismatch (", nrow(dim_reduced_data),
                 ") with Seurat cell count (", length(base_cell_ids), ")."), call. = FALSE)
    }
  }
  if (any(is.na(rownames(dim_reduced_data))) || anyDuplicated(rownames(dim_reduced_data))) {
    stop(paste("rownames(dim_reduced_data) for '", reduced_dim_name, "' contain NA or duplicated values."), call. = FALSE)
  }
  if (any(!is.finite(dim_reduced_data))) {
    num_non_finite <- sum(!is.finite(dim_reduced_data))
    warning(paste0("Non-finite values found in '", reduced_dim_name, "' coordinates (", num_non_finite, " instances)."), call. = FALSE)
  }
  
  if(is.null(colnames(counts_data))) stop("colnames(counts_data) is NULL.")
  if(!is.vector(cluster_labels_vector) && !is.factor(cluster_labels_vector)) stop("Metadata column '", cluster_col, "' is not a vector or factor.")
  
  if(any(counts_data@x < 0, na.rm=TRUE)) stop("Counts data for assay '", counts_assay_name, "' contain negative values.")
  if(!all(is.finite(counts_data@x))) warning("Non-finite values found in counts data for assay '", counts_assay_name, "'.", call. = FALSE)
  # Check for non-integers that are not NA before rounding
  non_int_idx <- counts_data@x %% 1 != 0 & !is.na(counts_data@x)
  if(any(non_int_idx)){
    message("Warning: Counts data for assay '", counts_assay_name, "' are not all integers. Rounding ", sum(non_int_idx) ," non-integer values.")
    counts_data@x[non_int_idx] <- round(counts_data@x[non_int_idx])
  }
  
  cluster_labels <- factor(cluster_labels_vector)
  names(cluster_labels) <- rownames(seurat_obj@meta.data) # Names from validated & aligned meta.data rownames
  
  start_cluster_char <- if (!is.null(start_cluster)) as.character(start_cluster) else NULL
  end_clusters_char <- if (!is.null(end_clusters)) as.character(end_clusters) else NULL
  
  if (!is.null(start_cluster_char) && !start_cluster_char %in% levels(cluster_labels)) {
    stop(paste("Start cluster '", start_cluster_char, "' not found in levels of '", cluster_col, "'. Available: ", paste(levels(cluster_labels), collapse=", ")), call. = FALSE)
  }
  
  
  message("--- Step 2: Aligning cell identifiers for Slingshot input ---")
  cells_in_reduction <- rownames(dim_reduced_data)
  cells_with_clusters <- names(cluster_labels)
  
  common_slingshot_input_cells <- intersect(cells_in_reduction, cells_with_clusters)
  message("Cells in reduction '", reduced_dim_name, "': ", length(cells_in_reduction),
          ". Cells with cluster labels in '", cluster_col, "': ", length(cells_with_clusters),
          ". Intersection for Slingshot input: ", length(common_slingshot_input_cells))
  
  if(length(common_slingshot_input_cells) < 2){
    stop(paste0("Too few common cells (", length(common_slingshot_input_cells), ") for Slingshot. Need at least 2."), call. = FALSE)
  }
  
  dim_reduced_data_filt <- dim_reduced_data[common_slingshot_input_cells, , drop = FALSE]
  cluster_labels_filt <- cluster_labels[common_slingshot_input_cells]
  
  if (!identical(rownames(dim_reduced_data_filt), names(cluster_labels_filt))) {
    stop("Internal error: Cell ID mismatch/order diff between filtered reduced_dims and cluster_labels pre-Slingshot.", call. = FALSE)
  }
  
  if(any(is.na(cluster_labels_filt))){
    message("WARNING: NA values found in cluster_labels_filt for ", sum(is.na(cluster_labels_filt)), " cells. These will be ignored by Slingshot.")
  }
  
  cluster_labels_filt <- droplevels(cluster_labels_filt)
  cluster_counts <- table(cluster_labels_filt)
  message("Cell counts per cluster for Slingshot input (total ", sum(cluster_counts) ," cells):")
  print(cluster_counts)
  if(any(cluster_counts == 0)) {
    warning("Some cluster levels have zero cells after filtering: ",
            paste(names(cluster_counts)[cluster_counts == 0], collapse=", "), call. = FALSE)
  }
  if(!is.null(start_cluster_char) && (!start_cluster_char %in% names(cluster_counts) || cluster_counts[start_cluster_char] == 0)){
    stop(paste("Start cluster '", start_cluster_char, "' has no cells in filtered data for Slingshot."), call.=FALSE)
  }
  
  
  message("--- Step 3: Running Slingshot ---")
  # Prepare base arguments for slingshot (controlled by the wrapper)
  slingshot_base_args <- list(
    data = dim_reduced_data_filt,
    clusterLabels = cluster_labels_filt
  )
  if(!is.null(start_cluster_char)) slingshot_base_args$start.clus <- start_cluster_char
  if(!is.null(end_clusters_char)) slingshot_base_args$end.clus <- end_clusters_char
  
  
  # Capture arguments passed in ... by the user, intended for slingshot::slingshot itself
  additional_user_args_for_slingshot <- list(...)
  
  message("--- Debug: Contents of '...' (ellipsis) passed to run_slingshot_from_seurat ---")
  if (length(additional_user_args_for_slingshot) == 0) {
    message("Ellipsis '...' is empty (as expected for a direct call with all main args specified).")
  } else {
    message("Ellipsis '...' contains the following named arguments that will be passed to slingshot:")
    print(names(additional_user_args_for_slingshot))
    # This is where we previously saw 'count_assay_name'. If it appears here again,
    # it means the call to run_slingshot_from_seurat is somehow including it in '...'
    # despite 'count_assay_name' being a formal argument.
    if ("count_assay_name" %in% names(additional_user_args_for_slingshot)) {
      warning("CRITICAL DEBUG: 'count_assay_name' was found in '...' (ellipsis_args). ",
              "This should NOT happen if 'run_slingshot_from_seurat' is called directly and 'count_assay_name' ",
              "is matched as a formal argument. Value from ellipsis: '",
              additional_user_args_for_slingshot$count_assay_name, "'.", call. = FALSE)
      # To prevent error, remove it from ... if it's there due to some R quirk
      # additional_user_args_for_slingshot[["count_assay_name"]] <- NULL
    }
  }
  
  final_slingshot_args <- c(slingshot_base_args, additional_user_args_for_slingshot)
  
  message("--- Final argument names being passed to slingshot::slingshot via do.call ---")
  message(paste(names(final_slingshot_args), collapse = ", "))
  
  
  slingshot_result <- tryCatch({
    do.call(slingshot::slingshot, final_slingshot_args)
  }, error = function(e) {
    message("ERROR during slingshot::slingshot() [called via do.call]: ", e$message)
    if(grepl("unique starting cluster", e$message, ignore.case = TRUE) && is.null(start_cluster_char)){
      message("Slingshot could not automatically identify a unique starting cluster. ",
              "Consider specifying 'start_cluster'. Available clusters with cells: ",
              paste(names(cluster_counts)[cluster_counts > 0], collapse=", "))
    }
    return(NULL)
  })
  
  if (is.null(slingshot_result)) return(NULL)
  message("Slingshot completed successfully.")
  
  # --- Step 4: Processing Slingshot output ---
  message("Processing Slingshot output...")
  pst_matrix <- slingshot::slingPseudotime(slingshot_result, na = FALSE)
  weights_matrix <- slingshot::slingCurveWeights(slingshot_result, na = FALSE)
  # Cell names for slingshot_result are based on the input dim_reduced_data_filt
  original_cell_names_for_slingshot_input <- rownames(dim_reduced_data_filt)
  
  if(nrow(pst_matrix) == 0 || ncol(pst_matrix) == 0){
    warning("Slingshot returned an empty pseudotime matrix. No trajectories found. Returning NULL.", call. = FALSE)
    return(NULL)
  }
  
  map_slingshot_rownames <- function(matrix_from_slingshot, original_names) {
    if (is.null(original_names) || length(original_names) == 0) {
      warning("map_slingshot_rownames: original_names is NULL or empty.", call. = FALSE)
      return(matrix_from_slingshot)
    }
    if (nrow(matrix_from_slingshot) > 0) {
      current_rn <- rownames(matrix_from_slingshot)
      if (!is.null(current_rn) && all(grepl("^[0-9]+$", current_rn)) && all(!is.na(as.integer(current_rn)))) {
        indices <- as.integer(current_rn)
        valid_idx_mask <- indices >= 1 & indices <= length(original_names) & !is.na(indices)
        if(sum(!valid_idx_mask) > 0) {
          warning("Some Slingshot output rowname indices are out of bounds or NA: ",
                  sum(!valid_idx_mask), " affected.", call.=FALSE)
        }
        matrix_from_slingshot <- matrix_from_slingshot[valid_idx_mask, , drop=FALSE]
        rownames(matrix_from_slingshot) <- original_names[indices[valid_idx_mask]]
      } else if (is.null(current_rn) && nrow(matrix_from_slingshot) == length(original_names)) {
        rownames(matrix_from_slingshot) <- original_names
      } else if (!is.null(current_rn) && !all(current_rn %in% original_names)) {
        # Names might already be character but a subset or different set
        # Filter to keep only those in original_names
        keep_rn <- current_rn[current_rn %in% original_names]
        if(length(keep_rn) < length(current_rn)) {
          warning("Some Slingshot output rownames were not in the original input cell names for Slingshot. Filtering.", call.=FALSE)
        }
        matrix_from_slingshot <- matrix_from_slingshot[keep_rn, , drop=FALSE]
      }
    }
    return(matrix_from_slingshot)
  }
  
  pst_matrix <- map_slingshot_rownames(pst_matrix, original_cell_names_for_slingshot_input)
  weights_matrix <- map_slingshot_rownames(weights_matrix, original_cell_names_for_slingshot_input)
  
  if(nrow(pst_matrix) == 0 || nrow(weights_matrix) == 0 || is.null(rownames(pst_matrix)) || is.null(rownames(weights_matrix)) ){
    warning("Pseudotime or weights matrix became empty or lost rownames after mapping. Returning NULL.", call.=FALSE)
    return(NULL)
  }
  
  cells_in_slingshot_output <- rownames(pst_matrix) # These should be character cell IDs from original input to slingshot
  if(length(cells_in_slingshot_output) == 0){
    warning("Slingshot output has no cells after mapping rownames. Returning NULL.", call. = FALSE)
    return(NULL)
  }
  
  
  # --- Step 5: Creating SingleCellExperiment object ---
  message("Creating SingleCellExperiment object...")
  # Final cells for SCE: must be in original Seurat object (base_cell_ids) AND in slingshot output
  final_common_cells <- intersect(cells_in_slingshot_output, base_cell_ids)
  
  if(length(final_common_cells) == 0){
    warning("No consistent cells to create SCE object (Slingshot output cells not in original Seurat object). Returning NULL.", call. = FALSE)
    return(NULL)
  }
  message("Finalizing SCE object with ", length(final_common_cells), " cells.")
  
  counts_data_sce <- counts_data[, final_common_cells, drop = FALSE]
  dim_reduced_data_sce <- dim_reduced_data[final_common_cells, , drop = FALSE]
  pst_matrix_sce <- pst_matrix[final_common_cells, , drop = FALSE]
  weights_matrix_sce <- weights_matrix[final_common_cells, , drop = FALSE]
  original_metadata_sce <- seurat_obj@meta.data[final_common_cells, , drop = FALSE]
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_data_sce),
    colData = S4Vectors::DataFrame(original_metadata_sce) # Start with original metadata for the final cells
  )
  
  tryCatch({
    reducedDim(sce, toupper(reduced_dim_name)) <- dim_reduced_data_sce
    reducedDim(sce, "slingshot") <- pst_matrix_sce
  }, error = function(e){
    warning("Failed to set reducedDims on SCE object: ", e$message, call.=FALSE)
  })
  
  weights_df <- S4Vectors::DataFrame(weights_matrix_sce)
  lineage_colnames_from_pst <- colnames(pst_matrix_sce)
  if (ncol(weights_df) > 0 && !is.null(lineage_colnames_from_pst) && ncol(weights_df) == ncol(pst_matrix_sce)) {
    colnames(weights_df) <- lineage_colnames_from_pst
  } else if (ncol(weights_df) > 0) {
    colnames(weights_df) <- paste0("Lineage", 1:ncol(weights_df))
  }
  
  # Carefully merge weights_df into existing colData(sce), which already has original_metadata_sce
  for(w_col in colnames(weights_df)){
    safe_w_col_name <- if(w_col %in% colnames(colData(sce))) paste0("slingWeight_", w_col) else w_col
    tryCatch(colData(sce)[[safe_w_col_name]] <- weights_df[[w_col]],
             error = function(e) warning("Could not add weight column '", safe_w_col_name, "' to colData: ", e$message, call. = FALSE))
  }
  
  # Ensure the specific cluster_col is present and is a factor with original levels
  if(cluster_col %in% colnames(original_metadata_sce)){
    colData(sce)[[cluster_col]] <- factor(original_metadata_sce[[cluster_col]],
                                          levels = levels(factor(seurat_obj@meta.data[[cluster_col]])))
  } else {
    warning("Original cluster column '", cluster_col, "' was lost from colData of SCE.", call. = FALSE)
  }
  
  message("SCE object created successfully with ", ncol(sce), " cells.")
  if(!is.null(slingshot_result) && "slingLineages" %in% ls(getNamespace("slingshot"))){ # Check if SlingshotDataSet was returned
    tryCatch({
      lineages_found <- slingshot::slingLineages(slingshot_result)
      message("Number of lineages found by Slingshot: ", length(lineages_found))
    }, error = function(e) message("Could not retrieve lineage count from Slingshot result."))
  }
  return(sce)
}


# 필요한 라이브러리 로드
# library(monocle3)
# library(mgcv)
# library(ggplot2)
# library(dplyr)
# library(purrr) # compact 함수 사용
# library(parallel) # mclapply 함수 사용 (Unix 계열)

# 플롯 저장 시 디렉토리 생성 및 파일명 중복 처리 헬퍼 함수
#' @export
save_plot_with_conflict_resolution <- function(plot_object, base_filename, output_dir, width, height, dpi) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Output directory created: ", output_dir)
  }
  
  filename_no_ext <- tools::file_path_sans_ext(base_filename)
  ext <- tools::file_ext(base_filename)
  if (ext == "") ext <- "png" # 기본 확장자
  
  full_path <- file.path(output_dir, paste0(filename_no_ext, ".", ext))
  counter <- 1
  while (file.exists(full_path)) {
    full_path <- file.path(output_dir, paste0(filename_no_ext, "_", counter, ".", ext))
    counter <- counter + 1
  }
  
  ggsave(filename = full_path, plot = plot_object, width = width, height = height, dpi = dpi)
  message("Plot saved: ", full_path)
  return(full_path)
}






# Roxygen2 주석 및 개선된 analyze_gene_dynamics 함수
# (이전 답변의 함수를 기반으로 수정)

#' Analyze gene expression dynamics along pseudotime using GAM
#'
#' This function fits a Generalized Additive Model (GAM) to model gene expression
#' changes along pseudotime, potentially accounting for different conditions.
#' It calculates various metrics and generates plots.
#'
#' @param gene_id A character string specifying the gene ID (e.g., gene symbol).
#' @param cds_obj A cell_data_set object from Monocle3.
#' @param condition_col_name A character string naming the column in `colData(cds_obj)`
#'   that contains the condition information (e.g., treatment group, patient prognosis).
#'   This column should ideally be a factor.
#' @param pseudotime_method A function or character string to specify how pseudotime is
#'   extracted. Defaults to `monocle3::pseudotime`. If your CDS object might have
#'   multiple disconnected trajectories or specific roots, you might need a custom
#'   function or ensure `order_cells` was run appropriately.
#' @param sample_col_name Optional. A character string naming the column in `colData(cds_obj)`
#'   that contains sample/patient IDs. If provided, metrics can be potentially
#'   summarized or checked at this level in future extensions, or used for pseudo-bulk aggregation.
#'   Currently, this argument is noted for future patient-level aggregation development.
#' @param output_dir A character string specifying the directory to save plots.
#'   If it doesn't exist, it will be created recursively.
#' @param k_val An integer, the number of basis functions for GAM splines (mgcv's k).
#' @param min_cells_for_fit An integer, the minimum number of cells with finite
#'   expression and pseudotime values required to fit the GAM.
#' @param plot_split Logical, if TRUE, plots for different conditions will be faceted.
#'   If FALSE (default), conditions will be overlaid on the same plot distinguished by color.
#' @param plot_width Width of the saved plot in inches.
#' @param plot_height Height of the saved plot in inches.
#' @param plot_dpi DPI of the saved plot.
#' @param scale_DR Logical, if TRUE, Dynamic Range (DR) will be scaled by the mean
#'   expression of the gene in the cells used for fitting.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `gene`: The input gene ID.
#'     \item `status`: Character string, "success", "skipped_insufficient_data", etc.
#'     \item `metrics`: A list or tibble of calculated metrics (e.g., deviance explained,
#'       interaction p-value, TV, DR per condition).
#'     \item `plot_path`: Path to the saved plot file, if successful.
#'     \item `plot_object`: The ggplot object itself.
#'   }
#'   Returns `NULL` if the gene is not found or critical errors occur early.
#'
#' @import monocle3
#' @import mgcv
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom methods is
#' @importFrom stats anova as.formula family formula lm predict residuals
#' @importFrom grDevices dev.off png
#' @importFrom tools file_path_sans_ext file_ext
#'
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'cds' is a valid Monocle3 object with pseudotime calculated
#' # and 'colData(cds)$condition' exists.
#' # gene_list <- c("geneA", "geneB")
#' # results <- lapply(gene_list, function(g) {
#' #   analyze_gene_dynamics(
#' #     gene_id = g,
#' #     cds_obj = cds,
#' #     condition_col_name = "condition",
#' #     output_dir = "my_gam_plots"
#' #   )
#' # })
#' # metrics_df <- dplyr::bind_rows(lapply(results, function(r) {
#' #   if(r$status == "success") {
#' #     # Further processing of r$metrics might be needed to flatten it
#' #     # This is a simplified example
#' #     df <- as.data.frame(r$metrics[!sapply(r$metrics, is.data.frame)])
#' #     if(!is.null(r$metrics$TV_per_condition)) {
#' #        tv_wide <- tidyr::pivot_wider(r$metrics$TV_per_condition,
#' #                                      names_from = cond, values_from = TV, names_prefix = "TV_")
#' #        df <- cbind(df, tv_wide)
#' #     }
#' #     if(!is.null(r$metrics$DR_per_condition)) {
#' #        dr_wide <- tidyr::pivot_wider(r$metrics$DR_per_condition,
#' #                                      names_from = cond, values_from = DR, names_prefix = "DR_")
#' #        df <- cbind(df, dr_wide)
#' #     }
#' #     return(df)
#' #   }
#' #   return(data.frame(gene = r$gene, status = r$status))
#' # }))
#' }
analyze_gene_dynamics <- function(gene_id,
                                  cds_obj,
                                  condition_col_name,
                                  pseudotime_method = monocle3::pseudotime,
                                  sample_col_name = NULL, # For future patient-level aggregation
                                  output_dir = "gam_analysis_plots",
                                  k_val = 6,
                                  min_cells_for_fit = 30,
                                  plot_split = FALSE, # Default to combined plot
                                  plot_width = 7,
                                  plot_height = 5,
                                  plot_dpi = 300,
                                  scale_DR = TRUE) {
  
  # --- 0. Input Validation and Data Extraction ---
  if (!is(cds_obj, "cell_data_set")) {
    stop("cds_obj must be a Monocle3 cell_data_set object.")
  }
  if (!gene_id %in% rownames(counts(cds_obj))) {
    warning("Gene '", gene_id, "' not found in CDS object. Skipping.")
    return(list(gene = gene_id, status = "skipped_gene_not_found", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  # Pseudotime extraction
  pt_values <- tryCatch({
    if (is.character(pseudotime_method)) {
      eval(parse(text = pseudotime_method))(cds_obj)
    } else if (is.function(pseudotime_method)) {
      pseudotime_method(cds_obj)
    } else {
      stop("pseudotime_method must be a function or a character string naming a function.")
    }
  }, error = function(e) {
    warning("Could not extract pseudotime for gene '", gene_id, "'. Error: ", e$message)
    return(NULL)
  })
  if (is.null(pt_values)) return(list(gene = gene_id, status = "skipped_pseudotime_extraction_failed", metrics = NULL, plot_path = NULL, plot_object = NULL))
  if (all(is.infinite(pt_values))) {
    warning("All pseudotime values are Inf for gene '", gene_id, "'. Skipping.")
    return(list(gene = gene_id, status = "skipped_all_pseudotime_inf", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  
  cell_metadata <- colData(cds_obj)
  if (!condition_col_name %in% colnames(cell_metadata)) {
    warning("Condition column '", condition_col_name, "' not found in colData(cds_obj). Skipping gene '", gene_id, "'.")
    return(list(gene = gene_id, status = "skipped_condition_col_not_found", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  conditions_original <- cell_metadata[[condition_col_name]]
  if (!is.factor(conditions_original)) {
    conditions_original <- factor(conditions_original)
  }
  names(conditions_original) <- rownames(cell_metadata)
  
  # Sample column (optional, for future use)
  # samples <- if (!is.null(sample_col_name) && sample_col_name %in% colnames(cell_metadata)) {
  #   cell_metadata[[sample_col_name]]
  # } else { NULL }
  # if (!is.null(samples)) names(samples) <- rownames(cell_metadata)
  
  gene_expr_matrix <- counts(cds_obj)[gene_id, , drop = FALSE]
  
  # --- 1. Data Preparation ---
  common_cells <- Reduce(intersect, list(names(pt_values), names(conditions_original), colnames(gene_expr_matrix)))
  
  if (length(common_cells) == 0) {
    warning("Gene '", gene_id, "': No common cells with pseudotime, condition, and expression data. Skipping.")
    return(list(gene = gene_id, status = "skipped_no_common_cells", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  num_inf_pt <- sum(is.infinite(pt_values[common_cells]))
  if (num_inf_pt > 0) {
    message("Gene '", gene_id, "': ", num_inf_pt, " out of ", length(common_cells),
            " common cells have Inf pseudotime. These will be filtered out.")
  }
  
  dat <- tibble(
    cell_id = common_cells,
    expr = as.numeric(gene_expr_matrix[gene_id, common_cells]),
    pseudotime = pt_values[common_cells],
    cond = conditions_original[common_cells]
    # sample = if(!is.null(samples)) samples[common_cells] else NA
  ) %>%
    filter(is.finite(expr) & is.finite(pseudotime))
  
  if (nrow(dat) < min_cells_for_fit) {
    warning("Gene '", gene_id, "': Fewer than ", min_cells_for_fit, " finite data points (actual: ", nrow(dat), "). Skipping GAM fit.")
    return(list(gene = gene_id, status = "skipped_insufficient_data", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  # Ensure 'cond' is a factor with potentially dropped levels
  dat$cond <- factor(dat$cond)
  if (nlevels(dat$cond) < 1) { # Should not happen if nrow(dat) > 0
    warning("Gene '", gene_id, "': No conditions left after filtering. Skipping.")
    return(list(gene = gene_id, status = "skipped_no_conditions_after_filter", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  if (nlevels(dat$cond) < 2) {
    message("Gene '", gene_id, "': Only one condition level present ('", levels(dat$cond),
            "') in the filtered data for GAM. Interaction term will not be meaningful/estimable in the standard way. Fitting a simpler model.")
    # For single condition, fit model without 'cond' and interaction.
    fit_formula_str <- "expr ~ s(pseudotime, k = k_val, bs = \"cr\")"
    # Or, if you want to fit for that single condition to get TV/DR for it:
    # fit_formula_str <- "expr ~ s(pseudotime, k = k_val, bs = \"cr\")"
    # And then patternTest related logic will be skipped.
  } else {
    fit_formula_str <- "expr ~ s(pseudotime, k = k_val, bs = \"cr\") + cond + s(pseudotime, by = cond, k = k_val, bs = \"cr\")"
  }
  fit_formula <- as.formula(fit_formula_str)
  
  
  # --- 2. GAM 모델 피팅 ---
  fit <- tryCatch({
    gam(fit_formula, family = nb(link = "log"), data = dat, method = "REML")
  }, error = function(e) {
    warning("Gene '", gene_id, "': GAM fitting failed. Error: ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(list(gene = gene_id, status = "gam_fit_failed", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  # --- 3. 모델 비교 (Interaction 항 유무) ---
  interaction_p_value <- NA_real_ # NA_real_ 로 명시적 초기화
  
  if (nlevels(dat$cond) >= 2 && "cond" %in% all.vars(fit_formula) && any(grepl("by = cond", deparse(fit_formula)))) {
    fit0_formula_str <- "expr ~ s(pseudotime, k = k_val, bs = \"cr\") + cond"
    fit0 <- tryCatch({
      gam(as.formula(fit0_formula_str), family = nb(link="log"), data = dat, method = "REML")
    }, error = function(e) {
      warning("Gene '", gene_id, "': GAM fitting for null model (fit0) failed. Error: ", e$message)
      return(NULL)
    })
    
    if (!is.null(fit0)) {
      interaction_anova_result <- tryCatch(anova(fit0, fit, test = "Chisq"), error = function(e) {
        warning("Gene '", gene_id, "': ANOVA test between models failed. Error: ", e$message)
        return(NULL)
      })
      if(!is.null(interaction_anova_result) && !is.null(interaction_anova_result$"P(>|Chi|)") && length(interaction_anova_result$"P(>|Chi|)") >= 2) {
        interaction_p_value <- interaction_anova_result$"P(>|Chi|)"[2]
      } #  else: interaction_p_value remains NA_real_
    } # else: interaction_p_value remains NA_real_
  } else {
    message("Gene '", gene_id, "': Interaction test skipped (e.g. single condition or model without interaction term).")
  }
  
  
  # --- 4. 예측 곡선 및 지표 계산 ---
  active_conditions <- levels(dat$cond)
  newd_list <- lapply(active_conditions, function(lvl) {
    pt_subset <- dat$pseudotime[dat$cond == lvl]
    if(length(pt_subset) < 2) return(NULL)
    data.frame(
      pseudotime = seq(min(pt_subset, na.rm = TRUE), max(pt_subset, na.rm = TRUE), length.out = 100),
      cond = factor(lvl, levels = active_conditions) # Ensure factor levels match for predict
    )
  })
  newd <- bind_rows(newd_list)
  
  if(nrow(newd) == 0) {
    warning("Gene '", gene_id, "': Could not generate prediction data.")
    gam_summary <- summary(fit) # Still get some summary
    metrics <- list(
      gene = gene_id, n_cells_fit = nrow(dat), deviance_explained = gam_summary$dev.expl,
      adj_r_squared = gam_summary$r.sq, interaction_p_value = interaction_p_value,
      smooth_term_summary = gam_summary$s.table, parametric_term_summary = gam_summary$p.table
    )
    return(list(gene = gene_id, status = "prediction_data_generation_failed", metrics = metrics, plot_path = NULL, plot_object = NULL))
  }
  
  newd$fit_response <- predict(fit, newdata = newd, type = "response")
  
  TV_per_cond <- newd %>%
    filter(is.finite(fit_response)) %>%
    group_by(cond) %>%
    summarize(TV = sum(abs(diff(fit_response))), .groups = "drop")
  
  DR_per_cond_calc <- newd %>%
    filter(is.finite(fit_response)) %>%
    group_by(cond) %>%
    summarize(DR_val = if(n() > 1 && (max(fit_response) - min(fit_response) > 0) ) max(fit_response) - min(fit_response) else 0,
              mean_expr_for_scaling = mean(dat$expr[dat$cond == first(cond)], na.rm = TRUE), # Mean of raw expr for that cond
              .groups = "drop")
  
  if(scale_DR){
    DR_per_cond <- DR_per_cond_calc %>%
      mutate(DR = ifelse(mean_expr_for_scaling > 1e-6, DR_val / mean_expr_for_scaling, DR_val)) %>% # Avoid division by zero
      select(cond, DR)
  } else {
    DR_per_cond <- DR_per_cond_calc %>% select(cond, DR = DR_val)
  }
  
  
  gam_summary <- summary(fit)
  dev_explained <- gam_summary$dev.expl
  adj_r_squared <- gam_summary$r.sq
  
  metrics <- list(
    gene = gene_id,
    n_cells_fit = nrow(dat),
    deviance_explained = dev_explained,
    adj_r_squared = adj_r_squared,
    interaction_p_value = interaction_p_value,
    TV_per_condition = if(nrow(TV_per_cond) > 0) TV_per_cond else NA,
    DR_per_condition = if(nrow(DR_per_cond) > 0) DR_per_cond else NA,
    smooth_term_summary = gam_summary$s.table, #DataFrame: edf, Ref.df, F, p-value for s(pseudotime), s(pseudotime, by=cond)
    parametric_term_summary = gam_summary$p.table #DataFrame: Estimate, Std. Error, t/z value, p-value for cond
  )
  
  # --- 5. 시각화 ---
  subtitle_parts <- c()
  if (is.data.frame(TV_per_cond) && nrow(TV_per_cond) > 0) {
    tv_strings <- apply(TV_per_cond, 1, function(r) paste0("TV(", r["cond"], "):", round(as.numeric(r["TV"]), 2)))
    subtitle_parts <- c(subtitle_parts, paste(tv_strings, collapse=", "))
  }
  if (is.data.frame(DR_per_cond) && nrow(DR_per_cond) > 0) {
    dr_strings <- apply(DR_per_cond, 1, function(r) paste0("DR(", r["cond"], "):", round(as.numeric(r["DR"]), 2)))
    subtitle_parts <- c(subtitle_parts, paste(dr_strings, collapse=", "))
  }
  plot_subtitle <- paste0(
    "Dev.Expl: ", round(dev_explained, 2),
    if(!is.na(interaction_p_value)) paste0(", P(int): ", format.pval(interaction_p_value, digits = 2, eps = 0.001)) else "",
    "\n", paste(subtitle_parts, collapse=" | ")
  )
  
  p <- ggplot() +
    geom_point(data = dat, aes(x = pseudotime, y = expr, color = cond), alpha = 0.2, size = 0.7, stroke=0) +
    geom_line(data = newd %>% filter(is.finite(fit_response)), aes(x = pseudotime, y = fit_response, color = cond), linewidth = 1) +
    labs(title = paste0("Gene: ", gene_id),
         subtitle = plot_subtitle,
         x = "Pseudotime", y = "Expression Level") +
    theme_classic(base_size = 10) +
    scale_color_discrete(name = if(nlevels(dat$cond) > 1) "Condition" else "") # Use actual colname or generic
  
  if(plot_split && nlevels(dat$cond) > 1){
    p <- p + facet_wrap(~cond, scales = "free_y", ncol = 1)
  }
  
  sanitized_gene_id <- gsub("[^a-zA-Z0-9_.-]", "_", gene_id)
  plot_base_filename <- paste0("GAM_dynamics_", sanitized_gene_id, ".png")
  
  # Using the helper function from your previous context, assuming it's defined
  # save_plot_with_conflict_resolution <- function(...) { ggsave(...); return(full_path)}
  plot_filepath <- tryCatch({
    save_plot_with_conflict_resolution( # Ensure this helper is available
      plot_object = p,
      base_filename = plot_base_filename,
      output_dir = output_dir,
      width = plot_width, height = plot_height, dpi = plot_dpi
    )}, error = function(e) {
      warning("Plot saving failed for gene '", gene_id, "': ", e$message)
      return(NA_character_)
    })
  
  
  return(list(gene = gene_id, status = "success", metrics = metrics, plot_path = plot_filepath, plot_object = p))
}


#' Helper to process a list of genes using analyze_gene_dynamics
#'
#' @param gene_list A character vector of gene IDs.
#' @param cds_obj A cell_data_set object from Monocle3.
#' @param condition_col_name Column name for conditions in `colData(cds_obj)`.
#' @param output_dir Directory to save plots.
#' @param num_cores Number of cores for parallel processing (uses `mclapply`).
#' @param ... Other arguments to pass to `analyze_gene_dynamics`.
#' @return A data.frame summarizing metrics for all successfully processed genes.
#'   Plots are saved to `output_dir`.
#' @export
#' @examples
#' \dontrun{
#' # summary_df <- process_gene_list_dynamics(
#' #   gene_list = c("geneA", "geneB"),
#' #   cds_obj = cds,
#' #   condition_col_name = "condition",
#' #   output_dir = "all_gam_plots",
#' #   num_cores = 2
#' # )
#' # print(summary_df)
#' }
process_gene_list_dynamics <- function(gene_list, cds_obj, condition_col_name, output_dir, num_cores = 1, ...) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  run_analysis_for_gene <- function(g) {
    analyze_gene_dynamics(
      gene_id = g,
      cds_obj = cds_obj,
      condition_col_name = condition_col_name,
      output_dir = output_dir,
      ...
    )
  }
  
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    all_results <- parallel::mclapply(gene_list, run_analysis_for_gene, mc.cores = num_cores)
  } else {
    all_results <- lapply(gene_list, run_analysis_for_gene)
  }
  
  # Filter out NULLs or complete failures if any
  valid_results <- purrr::compact(all_results) # Removes NULLs
  valid_results <- Filter(function(x) !is.null(x$status) && x$status == "success", valid_results)
  
  
  metrics_summary_list <- lapply(valid_results, function(res) {
    m <- res$metrics
    # Basic metrics
    base_df <- data.frame(
      gene = m$gene,
      n_cells_fit = m$n_cells_fit,
      deviance_explained = m$deviance_explained,
      adj_r_squared = m$adj_r_squared,
      interaction_p_value = m$interaction_p_value,
      stringsAsFactors = FALSE
    )
    
    # Flatten TV_per_condition
    if (is.data.frame(m$TV_per_condition) && nrow(m$TV_per_condition) > 0) {
      tv_wide <- tryCatch(tidyr::pivot_wider(m$TV_per_condition, names_from = cond, values_from = TV, names_prefix = "TV_"), error = function(e) NULL)
      if(!is.null(tv_wide) && ncol(tv_wide)>0) base_df <- cbind(base_df, tv_wide) else base_df$TV_data <- NA # Add placeholder if error
    } else {
      base_df$TV_data <- NA
    }
    
    
    # Flatten DR_per_condition
    if (is.data.frame(m$DR_per_condition) && nrow(m$DR_per_condition) > 0) {
      dr_wide <- tryCatch(tidyr::pivot_wider(m$DR_per_condition, names_from = cond, values_from = DR, names_prefix = "DR_"), error = function(e) NULL)
      if(!is.null(dr_wide) && ncol(dr_wide)>0) base_df <- cbind(base_df, dr_wide) else base_df$DR_data <- NA
    } else {
      base_df$DR_data <- NA
    }
    
    # Extract key p-values from smooth_term_summary if available
    # Example for s(pseudotime) overall effect (first row often)
    # And s(pseudotime, by=condition) (often multiple rows, one for each non-reference condition level)
    if(is.data.frame(m$smooth_term_summary) && nrow(m$smooth_term_summary) > 0){
      # Overall pseudotime effect (assuming first s(pseudotime) term is the global one)
      s_pt_row <- m$smooth_term_summary[grepl("^s\\(pseudotime\\)$", rownames(m$smooth_term_summary), ignore.case = TRUE), , drop = FALSE]
      if(nrow(s_pt_row) == 1) base_df$s_pseudotime_pval <- s_pt_row[1, "p-value"] else base_df$s_pseudotime_pval <- NA
      
      # Interaction p-values (example for one specific interaction term if name is known)
      # This can be complex if there are many 'by' variables or levels.
      # For simplicity, the overall interaction_p_value from ANOVA is often more direct.
    } else {
      base_df$s_pseudotime_pval <- NA
    }
    
    
    return(base_df)
  })
  
  # Combine all data frames
  if (length(metrics_summary_list) > 0) {
    final_summary_df <- dplyr::bind_rows(metrics_summary_list)
    return(final_summary_df)
  } else {
    return(data.frame()) # Return empty data frame if no successful results
  }
}


#' Analyze gene expression dynamics using tradeSeq
#'
#' Fits GAMs using the tradeSeq framework for a specific gene and performs tests
#' for differential progression or expression patterns along pseudotime trajectories.
#'
#' @param gene_id A character string, the ID of the gene to analyze.
#' @param sce_obj A SingleCellExperiment object, typically output from
#'   `run_slingshot_from_seurat` or similarly prepared. It must contain:
#'   - counts in `assays(sce_obj)$counts`
#'   - Slingshot pseudotime matrix (cells x lineages) in `reducedDim(sce_obj, "slingshot")`
#'   - Slingshot cell weights matrix (cells x lineages) as DataFrame columns in `colData(sce_obj)`
#'     (e.g., colData(sce_obj)$slingWeight1, colData(sce_obj)$slingWeight2, ...) or a matrix named 'slingshot_weights'.
#' @param condition_col A character string specifying the column name in `colData(sce_obj)`
#'   that contains the condition factor (e.g., treatment group, patient prognosis).
#' @param lineage_names A character vector specifying which lineages to test (e.g., "Lineage1", "Lineage2").
#'   If NULL (default), all lineages found in pseudotime matrix are considered if possible, or the first one.
#'   It's often better to specify. `tradeSeq` functions like `patternTest` can accept multiple lineages.
#' @param nknots An integer, the number of knots for GAM splines in `fitGAM`.
#' @param test_to_perform A character string. Currently supports "patternTest".
#'   Future versions could add "conditionTest", "diffEndTest", etc.
#' @param pseudotime_assay_name Character. Name of the reducedDim in SCE containing pseudotime. Default "slingshot".
#' @param weights_col_prefix Character. Prefix for colData columns containing cell weights. Default "slingWeight".
#'   Alternatively, if `colData(sce_obj)$slingshot_weights` is a matrix, it will be used.
#' @param output_dir A character string specifying the directory to save plots.
#' @param plot_split Logical, if TRUE, plots for different conditions will be faceted. Default FALSE.
#' @param scale_DR Logical, if TRUE, Dynamic Range (DR) will be scaled by the mean
#'   expression of the gene in the cells used for fitting for each lineage and condition.
#' @param fitGAM_args A list of additional arguments to pass to `tradeSeq::fitGAM`.
#' @param test_args A list of additional arguments to pass to the chosen test function (e.g., `patternTest`).
#'
#' @return A list containing:
#'   \itemize{
#'     \item `gene`: The input gene ID.
#'     \item `status`: Character string, "success" or an error/skip message.
#'     \item `metrics`: A list or tibble of calculated metrics (e.g., p-value from the test).
#'       DR and TV might be added if prediction extraction is robustly implemented.
#'     \item `plot_path`: Path to the saved plot file, if successful.
#'     \item `plot_object`: The ggplot object itself.
#'   }
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import tradeSeq
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'sce_sling' is an SCE object from run_slingshot_from_seurat
#' # and colData(sce_sling)$my_condition exists.
#' # Make sure lineage names in pseudotime matrix are "Lineage1", "Lineage2", etc.
#' # if (ncol(reducedDim(sce_sling, "slingshot")) > 0) {
#' #   lineage_to_test <- colnames(reducedDim(sce_sling, "slingshot"))[1]
#' #   tradeSeq_result <- analyze_gene_dynamics_tradeSeq(
#' #     gene_id = "GeneA", # Replace with an actual gene in your data
#' #     sce_obj = sce_sling,
#' #     condition_col = "my_condition",
#' #     lineage_names = lineage_to_test, # Test a specific lineage
#' #     nknots = 5,
#' #     output_dir = "tradeSeq_gene_plots"
#' #   )
#' #   if (!is.null(tradeSeq_result$plot_object)) print(tradeSeq_result$plot_object)
#' # }
#' }
analyze_gene_dynamics_tradeSeq <- function(gene_id,
                                           sce_obj,
                                           condition_col,
                                           lineage_names = NULL,
                                           nknots = 6,
                                           test_to_perform = "patternTest",
                                           pseudotime_assay_name = "slingshot",
                                           weights_col_prefix = "slingWeight",
                                           output_dir = "tradeSeq_gene_plots",
                                           plot_split = FALSE,
                                           scale_DR = TRUE,
                                           fitGAM_args = list(),
                                           test_args = list()) {
  
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop("Package 'tradeSeq' is required. Please install it.", call. = FALSE)
  }
  
  # --- 0. Input Validation and Data Preparation ---
  if (!is(sce_obj, "SingleCellExperiment")) {
    stop("sce_obj must be a SingleCellExperiment object.", call. = FALSE)
  }
  if (!gene_id %in% rownames(sce_obj)) {
    warning("Gene '", gene_id, "' not found in rownames(sce_obj). Skipping.")
    return(list(gene = gene_id, status = "skipped_gene_not_found", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  if (!condition_col %in% colnames(colData(sce_obj))) {
    warning("Condition column '", condition_col, "' not found in colData(sce_obj). Skipping.")
    return(list(gene = gene_id, status = "skipped_condition_col_not_found", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  if (!pseudotime_assay_name %in% reducedDimNames(sce_obj)) {
    stop("Pseudotime assay '", pseudotime_assay_name, "' not found in reducedDims(sce_obj).", call. = FALSE)
  }
  
  counts_matrix <- assays(sce_obj)$counts[gene_id, , drop = FALSE]
  if(!all(counts_matrix >= 0, na.rm = TRUE)){
    warning("Counts for gene '", gene_id, "' contain negative values. tradeSeq expects non-negative counts.")
    return(list(gene = gene_id, status = "skipped_negative_counts", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  if(!all(counts_matrix %% 1 == 0, na.rm = TRUE)){
    message("Warning: Counts for gene '", gene_id, "' are not all integers. Rounding for tradeSeq.")
    counts_matrix <- round(counts_matrix)
  }
  
  
  conditions <- factor(colData(sce_obj)[[condition_col]])
  names(conditions) <- colnames(sce_obj)
  
  # Pseudotime and weights
  pst_matrix_full <- reducedDim(sce_obj, pseudotime_assay_name) # Cells x Lineages
  
  # Try to get weights matrix
  cell_weights_matrix_full <- NULL
  if ("slingshot_weights" %in% names(colData(sce_obj)) && is.matrix(colData(sce_obj)$slingshot_weights)) {
    cell_weights_matrix_full <- colData(sce_obj)$slingshot_weights
  } else {
    weight_cols <- grep(paste0("^", weights_col_prefix), colnames(colData(sce_obj)), value = TRUE)
    if (length(weight_cols) > 0) {
      cell_weights_matrix_full <- as.matrix(colData(sce_obj)[, weight_cols, drop = FALSE])
      # Ensure colnames match lineage names if possible, or are standardized
      if(ncol(cell_weights_matrix_full) == ncol(pst_matrix_full) && is.null(colnames(cell_weights_matrix_full)) && !is.null(colnames(pst_matrix_full))){
        colnames(cell_weights_matrix_full) <- colnames(pst_matrix_full)
      }
    }
  }
  if(is.null(cell_weights_matrix_full)){
    stop("Cell weights matrix could not be found or constructed using prefix '", weights_col_prefix,
         "' or as 'slingshot_weights' in colData(sce_obj).", call. = FALSE)
  }
  if(nrow(pst_matrix_full) != nrow(cell_weights_matrix_full) || !all(rownames(pst_matrix_full) %in% rownames(cell_weights_matrix_full))){
    stop("Mismatch in cells between pseudotime and cell weights matrices.", call.=FALSE)
  }
  
  # Align cells (pst, weights, counts, conditions)
  common_cells <- Reduce(intersect, list(rownames(pst_matrix_full), rownames(cell_weights_matrix_full), colnames(counts_matrix), names(conditions)))
  if(length(common_cells) < 10){ # Arbitrary minimum
    warning("Gene '", gene_id, "': Too few common cells (", length(common_cells), ") after aligning all data. Skipping.")
    return(list(gene = gene_id, status = "skipped_too_few_common_cells", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  counts_matrix <- counts_matrix[, common_cells, drop = FALSE]
  conditions <- conditions[common_cells]
  pst_matrix <- pst_matrix_full[common_cells, , drop = FALSE]
  cell_weights_matrix <- cell_weights_matrix_full[common_cells, , drop = FALSE]
  
  # Handle lineage_names
  available_lineages <- colnames(pst_matrix)
  if (is.null(lineage_names)) {
    lineage_names <- available_lineages
    if(length(lineage_names) == 0) {
      warning("Gene '", gene_id, "': No lineages found in pseudotime matrix. Skipping.")
      return(list(gene = gene_id, status = "skipped_no_lineages", metrics = NULL, plot_path = NULL, plot_object = NULL))
    }
    message("Analyzing all available lineages: ", paste(lineage_names, collapse=", "))
  } else {
    lineage_names <- intersect(lineage_names, available_lineages)
    if (length(lineage_names) == 0) {
      warning("Gene '", gene_id, "': Specified lineage_names not found in pseudotime matrix. Available: ",
              paste(available_lineages, collapse=", "), ". Skipping.")
      return(list(gene = gene_id, status = "skipped_invalid_lineages", metrics = NULL, plot_path = NULL, plot_object = NULL))
    }
  }
  pst_matrix <- pst_matrix[, lineage_names, drop = FALSE]
  cell_weights_matrix <- cell_weights_matrix[, lineage_names, drop = FALSE]
  
  
  # Filter cells with NA/Inf pseudotime or all-zero weights for the selected lineages
  # A cell is valid if its pseudotime is finite for AT LEAST ONE lineage it's weighted to,
  # and its weights are not all zero.
  valid_cell_filter <- apply(pst_matrix, 1, function(row_pt) any(is.finite(row_pt))) &
    apply(cell_weights_matrix, 1, function(row_w) sum(row_w, na.rm = TRUE) > 1e-6) # Check for non-zero weight sum
  
  if(sum(valid_cell_filter) < 10){ # Arbitrary minimum
    warning("Gene '", gene_id, "': Too few cells (", sum(valid_cell_filter), ") after filtering for finite pseudotime/weights. Skipping.")
    return(list(gene = gene_id, status = "skipped_too_few_valid_cells_for_fitGAM", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  counts_matrix_filt <- counts_matrix[, valid_cell_filter, drop = FALSE]
  conditions_filt <- conditions[valid_cell_filter]
  pst_matrix_filt <- pst_matrix[valid_cell_filter, , drop = FALSE]
  cell_weights_matrix_filt <- cell_weights_matrix[valid_cell_filter, , drop = FALSE]
  
  if(nlevels(conditions_filt) < 2 && test_to_perform %in% c("patternTest", "conditionTest")){
    message("Gene '", gene_id, "': Only one condition level present for tradeSeq test '", test_to_perform,
            "'. Test will likely not be meaningful for condition comparison. Consider `startVsEndTest` or fitting without conditions.")
    # Proceeding, but test results might be trivial or error out.
  }
  
  
  # --- 1. Fit GAM using tradeSeq::fitGAM ---
  message("Fitting GAM for gene '", gene_id, "' with tradeSeq...")
  # fitGAM expects genes in rows, cells in columns
  default_fitGAM_args <- list(
    counts = counts_matrix_filt,
    pseudotime = pst_matrix_filt,
    cellWeights = cell_weights_matrix_filt,
    conditions = conditions_filt,
    nknots = nknots,
    verbose = FALSE,
    parallel = FALSE # For single gene, parallel=TRUE might add overhead
  )
  current_fitGAM_args <- utils::modifyList(default_fitGAM_args, fitGAM_args)
  
  # Ensure 'conditions' argument is only passed if there are multiple condition levels
  # or if the test explicitly needs it. fitGAM can run without conditions.
  if (nlevels(conditions_filt) < 2 && !is.null(current_fitGAM_args$conditions)) {
    message("Note: Only one condition level. Running fitGAM without explicit conditions parameter for this model, unless test requires it.")
    # current_fitGAM_args$conditions <- NULL # This might be too aggressive, let tradeSeq handle it.
  }
  
  
  fit_res_sce <- tryCatch({
    do.call(tradeSeq::fitGAM, current_fitGAM_args)
  }, error = function(e) {
    message("fitGAM failed for gene '", gene_id, "': ", e$message)
    return(NULL)
  })
  
  if (is.null(fit_res_sce)) {
    return(list(gene = gene_id, status = "fitGAM_failed", metrics = NULL, plot_path = NULL, plot_object = NULL))
  }
  
  # --- 2. Perform requested test ---
  test_df <- NULL
  p_value_col_name <- "pvalue" # Default
  message("Performing test: '", test_to_perform, "' for gene '", gene_id, "'...")
  
  if (test_to_perform == "patternTest") {
    default_test_args <- list(object = fit_res_sce, conditions = conditions_filt, global = TRUE, pairwise = TRUE) # Pass conditions for pairwise if desired
    current_test_args <- utils::modifyList(default_test_args, test_args)
    if(nlevels(conditions_filt) < 2 && !is.null(current_test_args$conditions)){
      # If only one condition, patternTest against a flat line might be implied or conditions removed
      # current_test_args$conditions <- NULL # Let tradeSeq handle it
    }
    test_df <- tryCatch(do.call(tradeSeq::patternTest, current_test_args), error = function(e){ message(e); NULL})
    p_value_col_name <- "pvalue_lineage1" # patternTest often gives p-values per lineage if not global
    # If global=TRUE and pairwise=FALSE, it's often 'pvalue' or similar.
    # If pairwise=TRUE, it's more complex, e.g. 'pvalue_condA_vs_condB_lineageX'
    # For simplicity, we might need to pick a primary p-value or report all.
    # Let's assume we are interested in the first lineage's p-value or a global one.
    # This part needs careful handling based on exact test_args.
    if(!is.null(test_df) && "waldStat" %in% colnames(test_df)) p_value_col_name <- "pvalue" # A common output
    else if(!is.null(test_df) && paste0("pvalue_", lineage_names[1]) %in% colnames(test_df)) p_value_col_name <- paste0("pvalue_", lineage_names[1])
    
    
  } else if (test_to_perform == "conditionTest") { # conditionTest often uses conditions within fitGAM
    default_test_args <- list(object = fit_res_sce, l2fc = 0) # Example arg
    current_test_args <- utils::modifyList(default_test_args, test_args)
    test_df <- tryCatch(do.call(tradeSeq::conditionTest, current_test_args), error = function(e){ message(e); NULL})
    p_value_col_name <- "pvalue" # Check tradeSeq output for exact name
  } else if (test_to_perform == "startVsEndTest") {
    default_test_args <- list(object = fit_res_sce, lineages=TRUE) # lineages=TRUE for per-lineage test
    current_test_args <- utils::modifyList(default_test_args, test_args)
    test_df <- tryCatch(do.call(tradeSeq::startVsEndTest, current_test_args), error = function(e){ message(e); NULL})
    if(!is.null(test_df) && paste0("pvalue_",lineage_names[1]) %in% colnames(test_df)) p_value_col_name <- paste0("pvalue_",lineage_names[1]) else p_value_col_name <- "pvalue"
  } else {
    warning("Test '", test_to_perform, "' is not currently implemented in this wrapper. Skipping test.")
  }
  
  main_p_value <- NA_real_
  if (!is.null(test_df) && gene_id %in% rownames(test_df) && p_value_col_name %in% colnames(test_df)) {
    main_p_value <- test_df[gene_id, p_value_col_name]
  } else if (!is.null(test_df) && nrow(test_df)==1 && rownames(test_df)[1] == "1" && p_value_col_name %in% colnames(test_df)){ # if only one gene was fit
    main_p_value <- test_df[1, p_value_col_name]
    rownames(test_df)[1] <- gene_id # rename for consistency
  }
  
  
  # --- 3. Calculate DR, TV (from predicted values) ---
  # This is more complex with tradeSeq as it doesn't directly provide a simple predict interface for arbitrary new data for one gene easily.
  # For DR/TV, one might need to extract GAM objects or use specific tradeSeq prediction functions.
  # For now, we'll focus on p-values and plotting.
  # If predictSmooth or predictCells can be adapted:
  # pred_vals_list <- list()
  # for (lin in lineage_names) {
  #   for (cond_lvl in levels(conditions_filt)) {
  #       # Define pseudotime range for this lineage & condition
  #       idx_filt <- conditions_filt == cond_lvl & cell_weights_matrix_filt[, lin] > 0.5 # Example filter
  #       if(sum(idx_filt) < 2) next
  #       pt_for_pred <- pst_matrix_filt[idx_filt, lin]
  #       if(length(unique(pt_for_pred)) < 2) next
  #       
  #       # This requires the actual GAM object for the gene, which is within fit_res_sce
  #       # Or use predictCells if appropriate
  #       # yhat <- predictSmooth(models = fit_res_sce, gene = gene_id, pseudotime = pred_pt_seq, condition = cond_lvl, lineage = lin)
  #       # Placeholder for actual prediction logic
  #   }
  # }
  # TV_per_cond_lin <- ...
  # DR_per_cond_lin_scaled <- ...
  # For now, returning NA for these
  TV_metric <- NA
  DR_metric <- NA
  
  
  # --- 4. Plotting ---
  plot_object <- NULL
  plot_filepath <- NA_character_
  message("Generating plot for gene '", gene_id, "'...")
  tryCatch({
    p_smooth <- plotSmoothers(models = fit_res_sce,
                              counts = counts_matrix_filt, # Original counts for the gene
                              gene = gene_id,
                              conditions = conditions_filt,
                              alpha = 0.6, pointSize = 0.8) +
      ggtitle(paste0(gene_id, " (Test P: ", format.pval(main_p_value, digits=3, eps=0.001), ")")) +
      theme_classic(base_size=10)
    
    if(plot_split && nlevels(conditions_filt) > 1){
      p_smooth <- p_smooth + facet_wrap(~condition, scales="free_y")
    }
    plot_object <- p_smooth
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    sanitized_gene_id <- gsub("[^a-zA-Z0-9_.-]", "_", gene_id)
    base_fn <- paste0("tradeSeq_dynamics_", sanitized_gene_id, ".png")
    # Assuming save_plot_with_conflict_resolution is defined globally from previous context
    plot_filepath <- save_plot_with_conflict_resolution(
      plot_object = plot_object,
      base_filename = base_fn,
      output_dir = output_dir,
      width = 7, height = 5, dpi = 300
    )
  }, error = function(e) {
    message("Plotting failed for gene '", gene_id, "': ", e$message)
  })
  
  
  # --- 5. Return results ---
  metrics_out <- list(
    gene = gene_id,
    test_performed = test_to_perform,
    main_p_value = main_p_value,
    # TV = TV_metric, # Add when implemented
    # DR_scaled = DR_metric, # Add when implemented
    full_test_result = if(!is.null(test_df) && gene_id %in% rownames(test_df)) test_df[gene_id, , drop=FALSE] else if (!is.null(test_df) && nrow(test_df)==1) test_df[1,,drop=FALSE] else NA
  )
  
  return(list(gene = gene_id, status = "success", metrics = metrics_out,
              plot_path = plot_filepath, plot_object = plot_object))
}