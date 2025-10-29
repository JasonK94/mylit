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
#'
#' @param n_samples Number of samples
#' @param format Format string with \%d placeholder (default: "Sample \%d")
#' @param include_doublets Whether to include doublet combinations (default: TRUE)
#'
#' @return Character vector of formatted sample names
#'
#' @examples
#' generate_sample_names(3, format = "Patient %d")
#' # Returns: c("Patient 1", "Patient 2", "Patient 3", 
#' #            "Patient 1+Patient 2", "Patient 1+Patient 3", "Patient 2+Patient 3")
#'
#' @export
generate_sample_names <- function(n_samples, 
                                  format = "Sample %d",
                                  include_doublets = TRUE) {
  
  # Generate singlet names
  singlets <- sprintf(format, seq_len(n_samples))
  
  if (!include_doublets) {
    return(singlets)
  }
  
  # Generate doublet combinations
  doublets <- character(0)
  if (n_samples >= 2) {
    for (i in seq_len(n_samples - 1)) {
      for (j in (i + 1):n_samples) {
        doublet_name <- paste0(sprintf(format, i), "+", sprintf(format, j))
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


