#' Demultiplexing Utility Functions
#'
#' This module provides utilities for demultiplexing single-cell data,
#' including barcode assignment and doublet detection.
#'
#' @name demulti_utils
NULL

#' Get Best Two Probabilities
#'
#' Finds the two highest probabilities from a vector.
#'
#' @param probs Numeric vector of probabilities
#'
#' @return List containing:
#'   \item{best}{Value of the highest probability}
#'   \item{best_idx}{Index of the highest probability}
#'   \item{second}{Value of the second highest probability}
#'   \item{second_idx}{Index of the second highest probability}
#'
#' @examples
#' get_best_two(c(0.1, 0.7, 0.2))
#' # Returns: list(best=0.7, best_idx=2, second=0.2, second_idx=3)
#'
#' @export
get_best_two <- function(probs) {
  sorted_idx <- order(probs, decreasing = TRUE)
  
  list(
    best = probs[sorted_idx[1]],
    best_idx = sorted_idx[1],
    second = probs[sorted_idx[2]],
    second_idx = sorted_idx[2]
  )
}

#' Get Barcode Mapping
#'
#' Assigns barcodes (cell identities) to samples based on probability matrix.
#'
#' @param prob_matrix Matrix where rows are cells and columns are samples,
#'   containing assignment probabilities
#' @param singlet_threshold Minimum probability for singlet assignment (default: 0.5)
#' @param doublet_threshold Minimum probability for doublet assignment (default: 0.3)
#' @param return_probs Whether to return probability values (default: FALSE)
#'
#' @return Data frame with columns:
#'   \item{barcode}{Cell barcode}
#'   \item{assignment}{Assigned sample(s)}
#'   \item{type}{"singlet" or "doublet"}
#'   \item{prob1}{(if return_probs=TRUE) Top probability}
#'   \item{prob2}{(if return_probs=TRUE) Second probability}
#'
#' @examples
#' \dontrun{
#' prob_matrix <- matrix(runif(100), nrow = 10, ncol = 10)
#' colnames(prob_matrix) <- paste0("Sample", 1:10)
#' rownames(prob_matrix) <- paste0("Cell", 1:10)
#' assignments <- get_barcode_mapping(prob_matrix)
#' }
#'
#' @export
get_barcode_mapping <- function(prob_matrix, 
                                singlet_threshold = 0.5, 
                                doublet_threshold = 0.3,
                                return_probs = FALSE) {
  
  if (!is.matrix(prob_matrix)) {
    stop("prob_matrix must be a matrix")
  }

  if (is.null(rownames(prob_matrix))) {
    rownames(prob_matrix) <- paste0("Cell_", seq_len(nrow(prob_matrix)))
  }
  
  if (is.null(colnames(prob_matrix))) {
    colnames(prob_matrix) <- paste0("Sample_", seq_len(ncol(prob_matrix)))
  }
  
  results <- data.frame(
    barcode = rownames(prob_matrix),
    assignment = character(nrow(prob_matrix)),
    type = character(nrow(prob_matrix)),
    stringsAsFactors = FALSE
  )
  
  if (return_probs) {
    results$prob1 <- numeric(nrow(prob_matrix))
    results$prob2 <- numeric(nrow(prob_matrix))
  }
  
  for (i in seq_len(nrow(prob_matrix))) {
    best_two <- get_best_two(prob_matrix[i, ])
    
    # Check for singlet
    if (best_two$best >= singlet_threshold) {
      results$assignment[i] <- colnames(prob_matrix)[best_two$best_idx]
      results$type[i] <- "singlet"
    } 
    # Check for doublet
    else if (best_two$best >= doublet_threshold && 
             best_two$second >= doublet_threshold) {
      sample1 <- colnames(prob_matrix)[best_two$best_idx]
      sample2 <- colnames(prob_matrix)[best_two$second_idx]
      results$assignment[i] <- paste(sort(c(sample1, sample2)), collapse = "+")
      results$type[i] <- "doublet"
    } 
    # Unassigned/negative
    else {
      results$assignment[i] <- "Negative"
      results$type[i] <- "negative"
    }
    
    if (return_probs) {
      results$prob1[i] <- best_two$best
      results$prob2[i] <- best_two$second
    }
  }
  
  return(results)
}

#' Check if Sample is Doublet
#'
#' Determines if a sample identifier represents a doublet (contains "+").
#'
#' @param sample_name Character vector of sample names
#'
#' @return Logical vector indicating doublet status
#'
#' @examples
#' is_doublet(c("Sample1", "Sample1+Sample2", "Sample3"))
#' # Returns: c(FALSE, TRUE, FALSE)
#'
#' @export
is_doublet <- function(sample_name) {
  grepl("+", sample_name, fixed = TRUE)
}

#' Generate Sample Values
#'
#' Creates a list of sample identifiers including both singlets and doublets.
#'
#' @param n_samples Number of samples
#' @param include_doublets Whether to include doublet combinations (default: TRUE)
#' @param prefix Prefix for sample names (default: "Sample")
#'
#' @return Character vector of sample identifiers
#'
#' @examples
#' generate_sample_values(3)
#' # Returns: c("Sample1", "Sample2", "Sample3", "Sample1+Sample2", 
#' #            "Sample1+Sample3", "Sample2+Sample3")
#'
#' @export
generate_sample_values <- function(n_samples, 
                                   include_doublets = TRUE,
                                   prefix = "Sample") {
  
  # Generate singlet names
  singlets <- paste0(prefix, seq_len(n_samples))
  
  if (!include_doublets) {
    return(singlets)
  }
  
  # Generate doublet combinations
  doublets <- character(0)
  if (n_samples >= 2) {
    for (i in seq_len(n_samples - 1)) {
      for (j in (i + 1):n_samples) {
        doublet_name <- paste0(prefix, i, "+", prefix, j)
        doublets <- c(doublets, doublet_name)
      }
    }
  }
  
  return(c(singlets, doublets))
}

#' Generate Sample Names (Formatted)
#'
#' Creates formatted sample name lists for display/plotting.
#' Supports both old-style (vector input) and new-style (n_samples/n_total_cols) usage.
#'
#' @param n_samples Number of samples (if NULL, will be inferred from n_total_cols or vector)
#' @param format Format string with \%d placeholder (default: "Sample \%d")
#' @param include_doublets Whether to include doublet combinations (default: TRUE)
#' @param n_total_cols Total number of columns (singlets + doublets). If provided, n_samples will be calculated.
#' @param start_num Starting number for sample names (default: 1)
#' @param vector Vector of sample numbers/names (for backward compatibility with old code)
#'
#' @return Character vector of formatted sample names
#'
#' @examples
#' generate_sample_names(3, format = "Patient %d")
#' # Returns: c("Patient 1", "Patient 2", "Patient 3", 
#' #            "Patient 1+Patient 2", "Patient 1+Patient 3", "Patient 2+Patient 3")
#' 
#' # Infer n_samples from total columns (e.g., from demux file)
#' generate_sample_names(n_total_cols = 21)  # 6 singlets + 15 doublets = 21
#' # Returns: c("Sample 1", "Sample 2", ..., "Sample 6", "Sample 1+Sample 2", ...)
#' 
#' # Old-style usage (backward compatibility)
#' generate_sample_names(vector = 1:6)
#' # Returns: c("1", "2", ..., "6", "1+2", "1+3", ...)
#'
#' @export
generate_sample_names <- function(n_samples = NULL, 
                                  format = "Sample %d",
                                  include_doublets = TRUE,
                                  n_total_cols = NULL,
                                  start_num = 1,
                                  vector = NULL) {
  
  # Backward compatibility: if vector is provided, use it
  if (!is.null(vector)) {
    if (is.numeric(vector)) {
      # Vector of numbers: generate names directly
      singlets <- as.character(vector)
      if (!include_doublets) {
        return(singlets)
      }
      # Generate doublet combinations
      doublets <- character(0)
      if (length(vector) >= 2) {
        for (i in seq_len(length(vector) - 1)) {
          for (j in (i + 1):length(vector)) {
            doublet_name <- paste0(vector[i], "+", vector[j])
            doublets <- c(doublets, doublet_name)
          }
        }
      }
      return(c(singlets, doublets))
    } else {
      # Vector of names: use as singlets
      singlets <- as.character(vector)
      if (!include_doublets) {
        return(singlets)
      }
      # Generate doublet combinations
      doublets <- character(0)
      if (length(vector) >= 2) {
        for (i in seq_len(length(vector) - 1)) {
          for (j in (i + 1):length(vector)) {
            doublet_name <- paste0(vector[i], "+", vector[j])
            doublets <- c(doublets, doublet_name)
          }
        }
      }
      return(c(singlets, doublets))
    }
  }
  
  # If n_total_cols is provided, calculate n_samples
  if (!is.null(n_total_cols) && is.null(n_samples)) {
    # n_total_cols = n_samples + nC2
    # n_total_cols = n_samples + n_samples * (n_samples - 1) / 2
    # Solve: n_total_cols = n + n*(n-1)/2
    # n_total_cols = n + (n^2 - n)/2 = (2n + n^2 - n)/2 = (n^2 + n)/2
    # n^2 + n - 2*n_total_cols = 0
    # n = (-1 + sqrt(1 + 8*n_total_cols)) / 2
    
    n_samples <- round((-1 + sqrt(1 + 8 * n_total_cols)) / 2)
    
    # Verify: n_samples + nC2 should equal n_total_cols
    expected_total <- n_samples + choose(n_samples, 2)
    if (abs(expected_total - n_total_cols) > 1) {
      warning(sprintf("Calculated n_samples=%d from n_total_cols=%d, but expected total=%d. Using calculated value anyway.",
                     n_samples, n_total_cols, expected_total))
    }
  }
  
  if (is.null(n_samples)) {
    stop("Either n_samples, n_total_cols, or vector must be provided")
  }
  
  # Generate singlet names
  singlets <- sprintf(format, start_num:(start_num + n_samples - 1))
  
  if (!include_doublets) {
    return(singlets)
  }
  
  # Generate doublet combinations
  doublets <- character(0)
  if (n_samples >= 2) {
    for (i in seq_len(n_samples - 1)) {
      for (j in (i + 1):n_samples) {
        doublet_name <- paste0(sprintf(format, start_num + i - 1), "+", 
                               sprintf(format, start_num + j - 1))
        doublets <- c(doublets, doublet_name)
      }
    }
  }
  
  return(c(singlets, doublets))
}

#' Parse Doublet Name
#'
#' Extracts the two sample identifiers from a doublet name.
#'
#' @param doublet_name Character string representing a doublet (e.g., "Sample1+Sample2")
#' @param separator Character used to separate samples (default: "+")
#'
#' @return Character vector of length 2 with the two sample names, or NULL if not a doublet
#'
#' @examples
#' parse_doublet_name("Sample1+Sample2")
#' # Returns: c("Sample1", "Sample2")
#'
#' @export
parse_doublet_name <- function(doublet_name, separator = "+") {
  if (!is_doublet(doublet_name)) {
    return(NULL)
  }
  
  parts <- strsplit(doublet_name, separator, fixed = TRUE)[[1]]
  
  if (length(parts) != 2) {
    warning("Unexpected doublet format: ", doublet_name)
    return(NULL)
  }
  
  return(trimws(parts))
}

#' Filter Demultiplexing Results
#'
#' Filters barcode assignments based on type and quality criteria.
#'
#' @param assignments Data frame from get_barcode_mapping
#' @param keep_singlets Keep singlet assignments (default: TRUE)
#' @param keep_doublets Keep doublet assignments (default: FALSE)
#' @param keep_negative Keep negative/unassigned cells (default: FALSE)
#' @param min_prob Minimum probability threshold (if prob columns exist)
#'
#' @return Filtered data frame
#'
#' @export
filter_demulti_results <- function(assignments,
                                   keep_singlets = TRUE,
                                   keep_doublets = FALSE,
                                   keep_negative = FALSE,
                                   min_prob = NULL) {
  
  # Filter by type
  keep_types <- character(0)
  if (keep_singlets) keep_types <- c(keep_types, "singlet")
  if (keep_doublets) keep_types <- c(keep_types, "doublet")
  if (keep_negative) keep_types <- c(keep_types, "negative")
  
  filtered <- assignments[assignments$type %in% keep_types, ]
  
  # Filter by probability if applicable
  if (!is.null(min_prob) && "prob1" %in% colnames(filtered)) {
    filtered <- filtered[filtered$prob1 >= min_prob, ]
  }
  
  return(filtered)
}

#' Summarize Demultiplexing Results
#'
#' Generates summary statistics for demultiplexing results.
#'
#' @param assignments Data frame from get_barcode_mapping
#'
#' @return List containing:
#'   \item{n_total}{Total number of cells}
#'   \item{n_singlets}{Number of singlet assignments}
#'   \item{n_doublets}{Number of doublet assignments}
#'   \item{n_negative}{Number of negative/unassigned cells}
#'   \item{doublet_rate}{Proportion of doublets}
#'   \item{assignment_rate}{Proportion successfully assigned}
#'
#' @export
summarize_demulti_results <- function(assignments) {
  n_total <- nrow(assignments)
  n_singlets <- sum(assignments$type == "singlet")
  n_doublets <- sum(assignments$type == "doublet")
  n_negative <- sum(assignments$type == "negative")
  
  list(
    n_total = n_total,
    n_singlets = n_singlets,
    n_doublets = n_doublets,
    n_negative = n_negative,
    doublet_rate = n_doublets / n_total,
    assignment_rate = (n_singlets + n_doublets) / n_total
  )
}

#' Demultiplex Demuxalot Output
#'
#' Processes demuxalot posterior CSV file and creates barcode mapping.
#' Automatically calculates n_samples from column count (demuxalot format: 1 BARCODE + n singlets + nC2 doublets).
#'
#' @param demuxalot_posterior Path to demuxalot posterior CSV file, or data frame (already loaded)
#' @param barcode_col Column name for barcodes (default: "BARCODE")
#' @param singlet_threshold Minimum probability for singlet (default: 0.5)
#' @param doublet_threshold Minimum probability for doublet (default: 0.3)
#' @param gem_name GEM name to add to metadata (optional)
#' @param gem_col Column name for GEM in output (default: "GEM")
#' @param return_probs Whether to return probability values (default: FALSE)
#'
#' @return Data frame with barcode mapping including:
#'   \item{Barcode}{Original barcode (with GEM suffix if gem_name provided)}
#'   \item{Best_Sample}{Assigned sample(s)}
#'   \item{droplet_demulti}{singlet_demulti or doublet_demulti}
#'   \item{GEM}{(if gem_name provided) GEM name}
#'   \item{Best_Probability}{(if return_probs=TRUE) Top probability}
#'   \item{Second_Best_Sample}{(if return_probs=TRUE) Second best sample}
#'   \item{Second_Best_Probability}{(if return_probs=TRUE) Second probability}
#'
#' @export
demultiplex_demuxalot <- function(demuxalot_posterior,
                                  barcode_col = "BARCODE",
                                  singlet_threshold = 0.5,
                                  doublet_threshold = 0.3,
                                  gem_name = NULL,
                                  gem_col = "GEM",
                                  return_probs = FALSE) {
  
  # Load data if path provided
  if (is.character(demuxalot_posterior)) {
    demux_data <- read.csv(demuxalot_posterior, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(demuxalot_posterior)) {
    demux_data <- demuxalot_posterior
  } else {
    stop("demuxalot_posterior must be a file path (character) or data frame")
  }
  
  # Calculate n_samples from column count
  # demuxalot format: 1 (BARCODE) + n (singlets) + nC2 (doublets) = 1 + n + n*(n-1)/2
  # n_total_cols = ncols - 1 (excluding BARCODE)
  n_total_cols <- ncol(demux_data) - 1
  n_samples <- round((-1 + sqrt(1 + 8 * n_total_cols)) / 2)
  
  # Verify calculation
  expected_total <- n_samples + choose(n_samples, 2)
  if (abs(expected_total - n_total_cols) > 1) {
    warning(sprintf("Calculated n_samples=%d from n_total_cols=%d, but expected total=%d. Using calculated value anyway.",
                   n_samples, n_total_cols, expected_total))
  }
  
  # Get all columns except barcode
  all_cols <- colnames(demux_data)
  if (!barcode_col %in% all_cols) {
    stop(sprintf("Barcode column '%s' not found in demux data. Found columns: %s",
                 barcode_col, paste(all_cols, collapse = ", ")))
  }
  
  # Get sample columns (all columns except barcode)
  sample_cols <- setdiff(all_cols, barcode_col)
  
  # Convert to matrix for get_barcode_mapping
  prob_matrix <- as.matrix(demux_data[, sample_cols, drop = FALSE])
  rownames(prob_matrix) <- demux_data[[barcode_col]]
  
  # Get barcode mapping
  barcode_map <- get_barcode_mapping(
    prob_matrix,
    singlet_threshold = singlet_threshold,
    doublet_threshold = doublet_threshold,
    return_probs = return_probs
  )
  
  # Rename assignment to Best_Sample (demuxalot convention)
  barcode_map$Best_Sample <- barcode_map$assignment
  barcode_map$assignment <- NULL
  
  # Add droplet_demulti tag
  barcode_map$droplet_demulti <- ifelse(
    is_doublet(barcode_map$Best_Sample),
    "doublet_demulti",
    "singlet_demulti"
  )
  
  # Add GEM if provided
  if (!is.null(gem_name)) {
    barcode_map[[gem_col]] <- gem_name
    # Create unique barcode with GEM suffix
    barcode_map$Barcode <- paste0(barcode_map$barcode, "_", gem_name)
  } else {
    barcode_map$Barcode <- barcode_map$barcode
  }
  
  # Set rownames to Barcode
  rownames(barcode_map) <- barcode_map$Barcode
  
  # Rename prob columns if return_probs
  if (return_probs) {
    # Get second best sample name (need to recalculate)
    prob_matrix_with_names <- prob_matrix
    colnames(prob_matrix_with_names) <- sample_cols
    
    second_best_samples <- character(nrow(barcode_map))
    for (i in seq_len(nrow(barcode_map))) {
      orig_barcode <- barcode_map$barcode[i]
      probs <- prob_matrix_with_names[orig_barcode, ]
      sorted_idx <- order(probs, decreasing = TRUE)
      if (length(sorted_idx) >= 2) {
        second_best_samples[i] <- sample_cols[sorted_idx[2]]
      } else {
        second_best_samples[i] <- NA
      }
    }
    barcode_map$Second_Best_Sample <- second_best_samples
    names(barcode_map)[names(barcode_map) == "prob1"] <- "Best_Probability"
    names(barcode_map)[names(barcode_map) == "prob2"] <- "Second_Best_Probability"
  }
  
  barcode_map
}

#' Demultiplex HTO Data
#'
#' Processes HTO data from filtered barcode matrix and creates barcode mapping.
#'
#' @param filtered_barcode_matrix Path to filtered barcode matrix directory, or Seurat object (already loaded)
#' @param hto_assay_name Name of HTO assay (default: "Multiplexing Capture")
#' @param method Demultiplexing method: "HTODemux" or "MULTIseqDemux" (default: "HTODemux")
#' @param positive_quantile Quantile for positive HTO signal (default: 0.99)
#' @param gem_name GEM name to add to metadata (optional)
#' @param gem_col Column name for GEM in output (default: "GEM")
#'
#' @return Data frame with barcode mapping including:
#'   \item{Barcode}{Cell barcode}
#'   \item{Best_Sample}{Assigned sample (HTO_maxID or MULTI_ID)}
#'   \item{droplet_demulti}{singlet_demulti or doublet_demulti}
#'   \item{GEM}{(if gem_name provided) GEM name}
#'   \item{Probability}{(if available) Assignment probability}
#'
#' @export
demultiplex_HTODemux <- function(filtered_barcode_matrix,
                                 hto_assay_name = "Multiplexing Capture",
                                 method = "HTODemux",
                                 positive_quantile = 0.99,
                                 gem_name = NULL,
                                 gem_col = "GEM") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for HTO demultiplexing")
  }
  
  # Load data if path provided
  if (is.character(filtered_barcode_matrix)) {
    filtered_counts_list <- Seurat::Read10X(filtered_barcode_matrix)
    
    # Extract Gene Expression and HTO
    if (is.list(filtered_counts_list)) {
      filtered_counts <- filtered_counts_list$`Gene Expression`
      hto_counts <- filtered_counts_list[[hto_assay_name]]
    } else {
      stop(sprintf("HTO assay '%s' not found in filtered barcode matrix", hto_assay_name))
    }
    
    # Create Seurat object
    obj <- Seurat::CreateSeuratObject(counts = filtered_counts)
    obj[[hto_assay_name]] <- Seurat::CreateAssayObject(counts = hto_counts)
    
    # Seurat converts spaces to dots
    actual_assay_name <- ifelse(hto_assay_name %in% names(obj@assays),
                               hto_assay_name,
                               gsub(" ", ".", hto_assay_name))
  } else if (inherits(filtered_barcode_matrix, "Seurat")) {
    obj <- filtered_barcode_matrix
    actual_assay_name <- ifelse(hto_assay_name %in% names(obj@assays),
                               hto_assay_name,
                               gsub(" ", ".", hto_assay_name))
  } else {
    stop("filtered_barcode_matrix must be a file path (character) or Seurat object")
  }
  
  # Perform demultiplexing
  if (method == "HTODemux") {
    # Filter out cells with zero HTO counts
    hto_assay <- obj[[actual_assay_name]]
    hto_counts <- Seurat::GetAssayData(hto_assay, slot = "counts")
    cells_with_hto <- colnames(hto_counts)[colSums(hto_counts) > 0]
    
    if (length(cells_with_hto) == 0) {
      stop("No cells with HTO counts found")
    }
    
    if (length(cells_with_hto) < ncol(obj)) {
      obj <- obj[, cells_with_hto]
    }
    
    # Try HTODemux
    tryCatch({
      obj <- Seurat::HTODemux(obj, assay = actual_assay_name, positive.quantile = positive_quantile)
    }, error = function(e) {
      warning(sprintf("HTODemux failed: %s. Trying with nstarts=100", e$message))
      tryCatch({
        obj <<- Seurat::HTODemux(obj, assay = actual_assay_name, positive.quantile = positive_quantile, nstarts = 100)
      }, error = function(e2) {
        warning(sprintf("HTODemux failed even with nstarts: %s. Assigning all cells as unknown", e2$message))
        obj$HTO_classification <<- "Unknown"
        obj$HTO_maxID <<- "unknown"
      })
    })
    
    # Create barcode mapping
    barcode_map <- data.frame(
      Barcode = colnames(obj),
      Best_Sample = obj$HTO_maxID,
      droplet_demulti = ifelse(obj$HTO_classification == "Doublet", "doublet_demulti", "singlet_demulti"),
      stringsAsFactors = FALSE
    )
    
    # Add probability if available
    if ("HTO_classification.global" %in% colnames(obj@meta.data)) {
      # HTODemux doesn't provide explicit probabilities, but we can use classification
      barcode_map$Probability <- ifelse(obj$HTO_classification == "Negative", 0, 1)
    }
    
  } else if (method == "MULTIseqDemux") {
    if (!requireNamespace("MULTIseq", quietly = TRUE)) {
      stop("MULTIseq package is required for MULTIseqDemux method")
    }
    obj <- MULTIseq::MULTIseqDemux(obj, assay = actual_assay_name, quantile = positive_quantile)
    
    # Create barcode mapping
    barcode_map <- data.frame(
      Barcode = colnames(obj),
      Best_Sample = obj$MULTI_ID,
      droplet_demulti = ifelse(obj$MULTI_classification == "Doublet", "doublet_demulti", "singlet_demulti"),
      stringsAsFactors = FALSE
    )
  } else {
    stop(sprintf("Unknown HTO demultiplexing method: %s. Use 'HTODemux' or 'MULTIseqDemux'", method))
  }
  
  # Add GEM if provided
  if (!is.null(gem_name)) {
    barcode_map[[gem_col]] <- gem_name
  }
  
  # Set rownames to Barcode
  rownames(barcode_map) <- barcode_map$Barcode
  
  barcode_map
}


