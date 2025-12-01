#' Demultiplexing Functions for Pipeline
#'
#' Functions for handling SNP and HTO demultiplexing
#'
#' @name pipe_demulti
NULL

#' Process Demuxalot Output (SNP-based)
#'
#' Processes demuxalot posterior CSV file and creates barcode mapping
#'
#' @param demux_file Path to demuxalot posterior CSV file
#' @param sample_names Vector of sample names (will be used as column names)
#' @param barcode_col Column name for barcodes (default: "BARCODE")
#' @param singlet_threshold Minimum probability for singlet (default: 0.5)
#' @param doublet_threshold Minimum probability for doublet (default: 0.3)
#' @param gem_name GEM name to add to metadata
#' @param gem_index Index for creating unique barcodes (default: 1)
#'
#' @return Data frame with barcode mapping including:
#'   \item{Barcode}{Unique barcode ID}
#'   \item{Best_Sample}{Assigned sample(s)}
#'   \item{droplet_demulti}{singlet_demulti or doublet_demulti}
#'   \item{GEM}{GEM name}
#'
#' @export
process_demuxalot <- function(demux_file, sample_names, barcode_col = "BARCODE",
                              singlet_threshold = 0.5, doublet_threshold = 0.3,
                              gem_name, gem_index = 1) {
  
  if (!requireNamespace("myR", quietly = TRUE)) {
    stop("myR package is required. Please load it with devtools::load_all()")
  }
  
  # Load demux data
  demux_data <- read.csv(demux_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Verify that sample_names match the columns in demux file
  # sample_names should already be the column names from the demux file
  all_cols <- colnames(demux_data)
  expected_cols <- c(barcode_col, sample_names)
  
  if (!all(expected_cols %in% all_cols)) {
    missing <- setdiff(expected_cols, all_cols)
    stop(sprintf("Mismatch: demux file columns don't match expected. Missing: %s. Found columns: %s",
                 paste(missing, collapse = ", "), paste(all_cols, collapse = ", ")))
  }
  
  # Use existing column names (don't rename)
  # The sample_names are already the correct column names from the demux file
  
  # Convert to matrix for get_barcode_mapping
  prob_matrix <- as.matrix(demux_data[, sample_names, drop = FALSE])
  rownames(prob_matrix) <- demux_data[[barcode_col]]
  
  # Get barcode mapping using myR function
  barcode_map <- myR::get_barcode_mapping(
    prob_matrix,
    singlet_threshold = singlet_threshold,
    doublet_threshold = doublet_threshold,
    return_probs = FALSE
  )
  
  # Rename assignment to Best_Sample (demuxalot convention)
  barcode_map$Best_Sample <- barcode_map$assignment
  barcode_map$assignment <- NULL
  
  # Add droplet_demulti tag
  barcode_map$droplet_demulti <- ifelse(
    myR::is_doublet(barcode_map$Best_Sample),
    "doublet_demulti",
    "singlet_demulti"
  )
  
  # Add GEM and create unique barcode
  barcode_map$GEM <- gem_name
  barcode_map$Barcode <- paste0(barcode_map$barcode, "_", gem_index)
  rownames(barcode_map) <- barcode_map$Barcode
  
  barcode_map
}

#' Process HTO Demultiplexing
#'
#' Performs HTO demultiplexing using Seurat's HTODemux or MULTIseqDemux
#'
#' @param obj Seurat object with HTO assay
#' @param method Demultiplexing method: "HTODemux" or "MULTIseqDemux"
#' @param assay_name Name of HTO assay (default: "Multiplexing Capture")
#' @param positive_quantile Quantile for positive HTO signal (default: 0.99)
#' @param sample_name Sample name to add to metadata
#' @param gem_name GEM name to add to metadata
#'
#' @return Seurat object with demultiplexing results in metadata
#'
#' @export
process_hto_demux <- function(obj, method = "HTODemux", assay_name = "Multiplexing Capture",
                              positive_quantile = 0.99, sample_name, gem_name) {
  
  if (!method %in% c("HTODemux", "MULTIseqDemux")) {
    stop(sprintf("Unknown HTO demultiplexing method: %s. Use 'HTODemux' or 'MULTIseqDemux'", method))
  }
  
  # Check if HTO assay exists
  if (!assay_name %in% names(obj@assays)) {
    stop(sprintf("HTO assay '%s' not found in Seurat object", assay_name))
  }
  
  # Perform demultiplexing
  if (method == "HTODemux") {
    # Filter out cells with zero HTO counts before demultiplexing
    hto_assay <- obj[[assay_name]]
    hto_counts <- Seurat::GetAssayData(hto_assay, slot = "counts")
    cells_with_hto <- colnames(hto_counts)[colSums(hto_counts) > 0]
    
    if (length(cells_with_hto) == 0) {
      stop("No cells with HTO counts found")
    }
    
    if (length(cells_with_hto) < ncol(obj)) {
      warning(sprintf("Filtering out %d cells with zero HTO counts", ncol(obj) - length(cells_with_hto)))
      obj <- obj[, cells_with_hto]
    }
    
    # Try HTODemux with error handling
    tryCatch({
      obj <- Seurat::HTODemux(obj, assay = assay_name, positive.quantile = positive_quantile)
      # HTODemux creates HTO_classification and HTO_maxID
      obj$sample_id <- obj$HTO_maxID
      obj$droplet_demulti <- ifelse(obj$HTO_classification == "Doublet", "doublet_demulti", "singlet_demulti")
    }, error = function(e) {
      # If HTODemux fails, try with nstarts parameter or fallback to simple thresholding
      warning(sprintf("HTODemux failed: %s. Trying with nstarts=100", e$message))
      tryCatch({
        obj <<- Seurat::HTODemux(obj, assay = assay_name, positive.quantile = positive_quantile, nstarts = 100)
        obj$sample_id <<- obj$HTO_maxID
        obj$droplet_demulti <<- ifelse(obj$HTO_classification == "Doublet", "doublet_demulti", "singlet_demulti")
      }, error = function(e2) {
        # Final fallback: assign all cells as unknown
        warning(sprintf("HTODemux failed even with nstarts: %s. Assigning all cells as unknown", e2$message))
        obj$sample_id <<- "unknown"
        obj$droplet_demulti <<- "unknown"
      })
    })
  } else if (method == "MULTIseqDemux") {
    if (!requireNamespace("MULTIseq", quietly = TRUE)) {
      stop("MULTIseq package is required for MULTIseqDemux method")
    }
    obj <- MULTIseq::MULTIseqDemux(obj, assay = assay_name, quantile = positive_quantile)
    # MULTIseq creates MULTI_ID and MULTI_classification
    obj$sample_id <- obj$MULTI_ID
    obj$droplet_demulti <- ifelse(obj$MULTI_classification == "Doublet", "doublet_demulti", "singlet_demulti")
  }
  
  # Add sample and GEM metadata
  obj$sample_id <- sample_name
  obj$GEM <- gem_name
  
  obj
}

#' Generate Sample Names
#'
#' Generates sample names in the format used by the pipeline
#'
#' @param n_samples Number of samples
#' @param start_num Starting number (default: 1)
#' @param format Format string (default: "Sample %d")
#'
#' @return Character vector of sample names
#'
#' @export
generate_sample_names <- function(n_samples, start_num = 1, format = "Sample %d") {
  sprintf(format, start_num:(start_num + n_samples - 1))
}

