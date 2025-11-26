#' Save CCI Results
#'
#' This module provides functions for saving and loading CCI analysis results.
#'
#' @importFrom dplyr %>%
#' @name save_cci_results
NULL

#' Save CCI Intermediate Results
#'
#' Saves intermediate results from CCI analysis preparation.
#' If output_dir contains run folders, saves to the latest run folder's prepared_data subdirectory.
#'
#' @param prepared_data List with prepared data (from prepare_cci_data functions)
#' @param output_dir Character string, output directory. If NULL, uses default "/data/user3/sobj/cci".
#'                   If provided and contains run folders, saves to run*/prepared_data/.
#' @param prefix Character string, filename prefix (default: "prepared_data")
#'
#' @return Character string, path to saved file
#' @export
save_cci_intermediate <- function(prepared_data, 
                                   output_dir = NULL,
                                   prefix = "prepared_data") {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for saving results")
  }
  
  if (is.null(output_dir)) {
    output_dir <- "/data/user3/sobj/cci"
  }
  
  # Check if output_dir contains run folders
  if (dir.exists(output_dir)) {
    existing_runs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
    run_numbers <- grep("^run[0-9]+$", existing_runs, value = TRUE) %>%
      sub("run", "", .) %>%
      as.numeric()
    
    if (length(run_numbers) > 0 && !all(is.na(run_numbers))) {
      # Save to the latest run folder's prepared_data subdirectory
      latest_run <- max(run_numbers, na.rm = TRUE)
      run_folder <- file.path(output_dir, paste0("run", latest_run))
      prepared_data_dir <- file.path(run_folder, "prepared_data")
      if (dir.exists(run_folder)) {
        output_dir <- prepared_data_dir
      }
    }
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
#' If output_dir is provided and contains run folders, saves to the latest run folder.
#' Otherwise saves to the specified output_dir or default cci directory.
#'
#' @param results List with CCI analysis results
#' @param output_dir Character string, output directory. If NULL, uses default "/data/user3/sobj/cci".
#'                   If provided and contains run folders, saves to the latest run folder.
#' @param prefix Character string, filename prefix (default: "nichenet_results")
#'
#' @return Character string, path to saved file
#' @export
save_cci_final <- function(results,
                            output_dir = NULL,
                            prefix = "nichenet_results") {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for saving results")
  }
  
  # Determine where to save
  # Priority: 1) output_dir from results$output_path, 2) provided output_dir, 3) default
  if (is.null(output_dir) && !is.null(results$output_path)) {
    output_dir <- results$output_path
  }
  
  if (is.null(output_dir)) {
    output_dir <- "/data/user3/sobj/cci"
  }
  
  # Check if output_dir contains run folders (from run_nichenet_analysis)
  if (dir.exists(output_dir)) {
    existing_runs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
    run_numbers <- grep("^run[0-9]+$", existing_runs, value = TRUE) %>%
      sub("run", "", .) %>%
      as.numeric()
    
    if (length(run_numbers) > 0 && !all(is.na(run_numbers))) {
      # Save to the latest run folder (where plots are already saved)
      latest_run <- max(run_numbers, na.rm = TRUE)
      run_folder <- file.path(output_dir, paste0("run", latest_run))
      if (dir.exists(run_folder)) {
        output_dir <- run_folder
      }
    }
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  filename <- paste0(prefix, ".qs")
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

