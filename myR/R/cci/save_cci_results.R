#' Save CCI Results
#'
#' This module provides functions for saving and loading CCI analysis results.
#'
#' @name save_cci_results
NULL

#' Save CCI Intermediate Results
#'
#' Saves intermediate results from CCI analysis preparation.
#'
#' @param prepared_data List with prepared data (from prepare_cci_data functions)
#' @param output_dir Character string, output directory (default: "/data/user3/sobj")
#' @param prefix Character string, filename prefix (default: "cci_prepared_data")
#'
#' @return Character string, path to saved file
#' @export
save_cci_intermediate <- function(prepared_data, 
                                   output_dir = "/data/user3/sobj",
                                   prefix = "cci_prepared_data") {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for saving results")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(prefix, "_", timestamp, ".qs")
  filepath <- file.path(output_dir, filename)
  
  tryCatch({
    qs::qsave(prepared_data, filepath)
    message("Intermediate results saved to: ", filepath)
    return(filepath)
  }, error = function(e) {
    warning("Failed to save intermediate results: ", e$message)
    return(NULL)
  })
}

#' Save CCI Final Results
#'
#' Saves final CCI analysis results.
#'
#' @param results List with CCI analysis results
#' @param output_dir Character string, output directory (default: "/data/user3/sobj")
#' @param prefix Character string, filename prefix (default: "cci_analysis_results")
#'
#' @return Character string, path to saved file
#' @export
save_cci_final <- function(results,
                            output_dir = "/data/user3/sobj",
                            prefix = "cci_analysis_results") {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for saving results")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(prefix, "_", timestamp, ".qs")
  filepath <- file.path(output_dir, filename)
  
  tryCatch({
    qs::qsave(results, filepath)
    message("Final results saved to: ", filepath)
    return(filepath)
  }, error = function(e) {
    warning("Failed to save final results: ", e$message)
    return(NULL)
  })
}

#' Load CCI Results
#'
#' Loads previously saved CCI analysis results.
#'
#' @param filepath Character string, path to saved results file
#'
#' @return List with loaded results
#' @export
load_cci_results <- function(filepath) {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for loading results")
  }
  
  if (!file.exists(filepath)) {
    stop("Results file not found: ", filepath)
  }
  
  tryCatch({
    results <- qs::qread(filepath)
    message("Results loaded from: ", filepath)
    return(results)
  }, error = function(e) {
    stop("Failed to load results: ", e$message)
  })
}

