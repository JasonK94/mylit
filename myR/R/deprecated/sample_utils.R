#' Sample Naming and Sorting Utilities
#'
#' This module provides utilities for handling sample identifiers and sorting.
#'
#' @name sample_utils
NULL

#' Sort Sample Identifiers
#'
#' Sorts character sample identifiers, handling both single-number and 
#' dual-number strings (e.g., "1", "2+3", "10"). Numeric samples are sorted
#' numerically, with doublets (containing "+") sorted after singlets.
#'
#' @param samples Character vector of sample identifiers
#'
#' @return Sorted character vector of sample IDs
#'
#' @examples
#' \dontrun{
#' samples <- c("10", "1+2", "2", "1", "11+12", "3")
#' sort_samples(samples)
#' # Returns: "1" "2" "3" "10" "1+2" "11+12"
#' }
#'
#' @export
sort_samples <- function(samples) {
  # Helper function to extract numbers from a string
  extract_numbers <- function(x) {
    parts <- strsplit(x, "\\+")[[1]]
    nums  <- as.numeric(parts)
    if (any(is.na(nums))) return(c(Inf, Inf))
    nums
  }
  
  # Check if all samples are non-numeric
  all_non_numeric <- all(sapply(samples, function(x) {
    !grepl("\\+", x) && is.na(suppressWarnings(as.numeric(x)))
  }))
  
  if (all_non_numeric) {
    return(sort(samples))  # If all non-numeric, return alphabetically sorted
  }
  
  # Separate samples into those with and without "+"
  is_doublet <- grepl("\\+", samples)
  single_samples <- samples[!is_doublet]
  doublet_samples <- samples[is_doublet]
  
  # Sort single samples numerically where possible
  single_sorted <- single_samples[order(
    suppressWarnings(as.numeric(single_samples)), 
    single_samples
  )]
  
  # Sort doublet samples by first number, then second number
  if (length(doublet_samples) > 0) {
    doublet_sorted <- doublet_samples[order(
      sapply(doublet_samples, function(x) extract_numbers(x)[1]),  # First number
      sapply(doublet_samples, function(x) extract_numbers(x)[2])   # Second number
    )]
  } else {
    doublet_sorted <- character(0)
  }
  
  # Combine and return
  return(c(single_sorted, doublet_sorted))
}

#' Generate Sample Names
#'
#' Helper function to generate standardized sample names.
#'
#' @param n_samples Number of samples to generate names for
#' @param prefix Prefix for sample names (default: "Sample")
#' @param start_index Starting index (default: 1)
#'
#' @return Character vector of sample names
#'
#' @examples
#' generate_sample_names(5)
#' # Returns: "Sample1" "Sample2" "Sample3" "Sample4" "Sample5"
#'
#' @export
generate_sample_names <- function(n_samples, prefix = "Sample", start_index = 1) {
  paste0(prefix, seq(start_index, start_index + n_samples - 1))
}

#' Parse Sample Identifiers
#'
#' Extracts numeric components from sample identifiers.
#'
#' @param samples Character vector of sample identifiers
#' @param return_doublets Whether to return doublet information
#'
#' @return Data frame with columns: sample, numeric1, numeric2 (if doublet), is_doublet
#'
#' @examples
#' \dontrun{
#' samples <- c("1", "2+3", "10")
#' parse_sample_ids(samples, return_doublets = TRUE)
#' }
#'
#' @export
parse_sample_ids <- function(samples, return_doublets = TRUE) {
  result <- data.frame(
    sample = samples,
    numeric1 = NA_real_,
    numeric2 = NA_real_,
    is_doublet = grepl("\\+", samples),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(samples)) {
    parts <- strsplit(samples[i], "\\+")[[1]]
    nums <- suppressWarnings(as.numeric(parts))
    
    if (!is.na(nums[1])) {
      result$numeric1[i] <- nums[1]
    }
    
    if (length(nums) > 1 && !is.na(nums[2])) {
      result$numeric2[i] <- nums[2]
    }
  }
  
  if (!return_doublets) {
    result <- result[!result$is_doublet, ]
  }
  
  return(result)
}


