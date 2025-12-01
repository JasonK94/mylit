#' Pipeline Utility Functions
#'
#' Helper functions for the single-cell RNA-seq processing pipeline
#'
#' @name pipe_utils
NULL

#' Load Configuration Files
#'
#' Loads config.csv, config_default.csv, and methods.csv into a list
#'
#' @param config_path Path to config.csv file
#' @param config_dir Directory containing config files (default: dirname of config_path)
#'
#' @return List containing:
#'   \item{config}{Data frame with sample-specific configuration}
#'   \item{defaults}{Named vector of default parameter values}
#'   \item{methods}{Data frame with method-specific parameters}
#'
#' @export
load_config <- function(config_path, config_dir = NULL) {
  if (is.null(config_dir)) {
    config_dir <- dirname(config_path)
  }
  
  # Load config.csv
  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s", config_path))
  }
  
  # Use readr for better CSV handling
  if (requireNamespace("readr", quietly = TRUE)) {
    config <- readr::read_csv(config_path, show_col_types = FALSE, comment = "", 
                              na = c("", "NA", "na"), trim_ws = TRUE)
    # Convert to data.frame for compatibility
    config <- as.data.frame(config)
  } else {
    # Fallback to base read.csv with fill=TRUE to handle missing columns
    config <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE, 
                       comment.char = "", quote = "\"", sep = ",", fill = TRUE,
                       na.strings = c("", "NA", "na"))
  }
  
  # Load config_default.csv
  default_path <- file.path(config_dir, "config_default.csv")
  if (!file.exists(default_path)) {
    warning(sprintf("Default config not found: %s. Using empty defaults.", default_path))
    defaults <- character(0)
  } else {
    if (requireNamespace("readr", quietly = TRUE)) {
      default_df <- readr::read_csv(default_path, show_col_types = FALSE, comment = "#", 
                                   na = c("", "NA", "na"), trim_ws = TRUE)
      default_df <- as.data.frame(default_df)
    } else {
      default_df <- read.csv(default_path, stringsAsFactors = FALSE, comment.char = "#",
                            na.strings = c("", "NA", "na"))
    }
    defaults <- setNames(default_df$value, default_df$parameter)
  }
  
  # Load methods.csv
  methods_path <- file.path(config_dir, "methods.csv")
  if (!file.exists(methods_path)) {
    warning(sprintf("Methods config not found: %s. Using empty methods.", methods_path))
    methods <- data.frame()
  } else {
    methods <- read.csv(methods_path, stringsAsFactors = FALSE)
  }
  
  list(config = config, defaults = defaults, methods = methods)
}

#' Get Parameter Value
#'
#' Gets a parameter value from config, with fallback to defaults
#'
#' @param param_name Name of the parameter
#' @param config_list List returned by load_config()
#' @param default_value Default value if not found in config or defaults
#'
#' @return Parameter value (character, numeric, or logical as appropriate)
#'
#' @export
get_param <- function(param_name, config_list, default_value = NULL) {
  # Try config defaults first
  if (param_name %in% names(config_list$defaults)) {
    val <- config_list$defaults[[param_name]]
    # Try to convert to appropriate type
    if (val == "TRUE") return(TRUE)
    if (val == "FALSE") return(FALSE)
    if (grepl("^\\d+$", val)) return(as.numeric(val))
    if (grepl("^\\d+:\\d+$", val)) {
      parts <- strsplit(val, ":")[[1]]
      return(as.numeric(parts[1]):as.numeric(parts[2]))
    }
    return(val)
  }
  
  # Fallback to provided default
  if (!is.null(default_value)) {
    return(default_value)
  }
  
  stop(sprintf("Parameter '%s' not found in config and no default provided", param_name))
}

#' Setup Logging
#'
#' Sets up logging for a pipeline run
#'
#' @param run_id Run identifier (e.g., "run1")
#' @param log_dir Base directory for logs (default: "logs")
#' @param master_log_file Name of master log file (default: "total.log")
#'
#' @return List with log file connections and paths
#'
#' @export
setup_logging <- function(run_id, log_dir = "logs", master_log_file = "total.log") {
  # Create directories
  run_log_dir <- file.path(log_dir, run_id)
  dir.create(run_log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Master log
  master_log_path <- file.path(log_dir, master_log_file)
  master_log <- file(master_log_path, open = "a")
  
  # Run log
  run_log_path <- file.path(run_log_dir, sprintf("%s_log.log", run_id))
  run_log <- file(run_log_path, open = "a")
  
  # Log initial message
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  writeLines(sprintf("[%s] Starting run: %s", timestamp, run_id), master_log)
  writeLines(sprintf("[%s] Starting run: %s", timestamp, run_id), run_log)
  
  list(
    master_log = master_log,
    run_log = run_log,
    master_log_path = master_log_path,
    run_log_path = run_log_path,
    run_log_dir = run_log_dir
  )
}

#' Log Message
#'
#' Writes a message to log files
#'
#' @param message Message to log
#' @param log_list List returned by setup_logging()
#' @param level Log level (INFO, WARNING, ERROR)
#'
#' @export
log_message <- function(message, log_list, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- sprintf("[%s] [%s] %s", timestamp, level, message)
  
  writeLines(formatted_msg, log_list$master_log)
  writeLines(formatted_msg, log_list$run_log)
  
  # Also print to console
  cat(formatted_msg, "\n")
}

#' Close Logging
#'
#' Closes log file connections
#'
#' @param log_list List returned by setup_logging()
#'
#' @export
close_logging <- function(log_list) {
  close(log_list$master_log)
  close(log_list$run_log)
}

#' Get Output Path
#'
#' Constructs output file path based on run_id and step
#'
#' @param run_id Run identifier
#' @param step Step number or name
#' @param filename Filename
#' @param output_base_dir Base directory for outputs (default: "/data/user3/sobj/pipe")
#'
#' @return Full path to output file
#'
#' @export
get_output_path <- function(run_id, step, filename, output_base_dir = "/data/user3/sobj/pipe") {
  step_dir <- file.path(output_base_dir, run_id, sprintf("step%d", step))
  dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(step_dir, filename)
}

#' Save Intermediate Object
#'
#' Saves an R object using qs for fast serialization
#'
#' @param obj Object to save
#' @param filepath Full path to output file
#' @param log_list Optional log list for logging
#'
#' @export
save_intermediate <- function(obj, filepath, log_list = NULL) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for fast serialization")
  }
  
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(obj, filepath)
  
  if (!is.null(log_list)) {
    log_message(sprintf("Saved intermediate object to: %s", filepath), log_list)
  }
}

#' Load Intermediate Object
#'
#' Loads an R object saved with qs
#'
#' @param filepath Full path to input file
#' @param log_list Optional log list for logging
#'
#' @return Loaded object
#'
#' @export
load_intermediate <- function(filepath, log_list = NULL) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for fast serialization")
  }
  
  if (!file.exists(filepath)) {
    stop(sprintf("Intermediate file not found: %s", filepath))
  }
  
  obj <- qs::qread(filepath)
  
  if (!is.null(log_list)) {
    log_message(sprintf("Loaded intermediate object from: %s", filepath), log_list)
  }
  
  obj
}

#' Extract Metadata from Filename
#'
#' Extracts metadata values from sample names using regex patterns
#'
#' @param sample_name Sample name
#' @param patterns Named list of regex patterns (names are metadata column names)
#' @param mappings Optional named list of mappings (e.g., list(day = list("24" = 1, "72" = 3)))
#'
#' @return Named vector of extracted metadata values
#'
#' @export
extract_metadata_from_filename <- function(sample_name, patterns, mappings = NULL) {
  result <- character(length(patterns))
  names(result) <- names(patterns)
  
  for (col_name in names(patterns)) {
    pattern <- patterns[[col_name]]
    match <- regmatches(sample_name, regexpr(pattern, sample_name))
    
    if (length(match) > 0) {
      value <- match[1]
      
      # Apply mapping if provided
      if (!is.null(mappings) && col_name %in% names(mappings)) {
        mapping <- mappings[[col_name]]
        if (value %in% names(mapping)) {
          value <- mapping[[value]]
        }
      }
      
      result[col_name] <- as.character(value)
    } else {
      result[col_name] <- NA_character_
    }
  }
  
  result
}

#' Calculate Ribosomal Gene Percentage
#'
#' Calculates percentage of ribosomal genes, excluding known non-ribosomal genes
#'
#' @param obj Seurat object
#' @param rps_pattern Pattern for RPS genes (default: "^RPS")
#' @param rpl_pattern Pattern for RPL genes (default: "^RPL")
#' @param assay Assay to use (default: "RNA")
#'
#' @return Seurat object with percent.ribo added to metadata
#'
#' @export
calculate_ribosomal_percentage <- function(obj, rps_pattern = "^RPS", rpl_pattern = "^RPL", assay = "RNA") {
  # Get all genes
  all_genes <- rownames(obj[[assay]])
  
  # Find RPS and RPL genes
  rps_genes <- grep(rps_pattern, all_genes, value = TRUE)
  rpl_genes <- grep(rpl_pattern, all_genes, value = TRUE)
  
  # Known non-ribosomal genes that match patterns (add more as needed)
  non_ribo <- c("RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6",
                 "RPL10A", "RPL10L")  # Add more known exceptions
  
  # Filter out non-ribosomal
  rps_genes <- setdiff(rps_genes, non_ribo)
  rpl_genes <- setdiff(rpl_genes, non_ribo)
  
  ribo_genes <- c(rps_genes, rpl_genes)
  
  if (length(ribo_genes) > 0) {
    obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, features = ribo_genes, assay = assay)
  } else {
    obj[["percent.ribo"]] <- 0
    warning("No ribosomal genes found matching patterns")
  }
  
  obj
}

