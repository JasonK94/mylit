#' Pipeline Utility Functions
#'
#' Helper functions for the single-cell RNA-seq processing pipeline
#'
#' @name pipe_utils
NULL

#' Load Configuration Files
#'
#' Loads Manifest, Pipeline Config, and Execution Config into a list
#'
#' @param config_path Path to the Manifest file (e.g., manifest_stroke.csv) or legacy config.csv
#' @param config_dir Directory containing config files (default: dirname of config_path)
#' @param execution_config_path Path to execution config (optional)
#'
#' @return List containing:
#'   \item{config}{Data frame with sample-specific configuration (Manifest)}
#'   \item{defaults}{Named vector of pipeline parameters (Pipeline Config)}
#'   \item{execution}{List of execution-specific settings}
#'   \item{methods}{Data frame with method-specific parameters}
#'
#' @export
load_config <- function(config_path, config_dir = NULL, execution_config_path = NULL) {
  if (is.null(config_dir)) {
    config_dir <- dirname(config_path)
  }

  # 1. Load Manifest (formerly config.csv)
  if (!file.exists(config_path)) {
    stop(sprintf("Manifest file not found: %s", config_path))
  }

  if (requireNamespace("readr", quietly = TRUE)) {
    config <- readr::read_csv(config_path,
      show_col_types = FALSE, comment = "",
      na = c("", "NA", "na"), trim_ws = TRUE
    )
    config <- as.data.frame(config)
  } else {
    config <- read.csv(config_path,
      stringsAsFactors = FALSE, check.names = FALSE,
      comment.char = "", quote = "\"", sep = ",", fill = TRUE,
      na.strings = c("", "NA", "na")
    )
  }

  # Clean up file paths and quotes
  char_cols <- sapply(config, is.character)
  for (col in names(char_cols)[char_cols]) {
    config[[col]] <- gsub('^["\']+|["\']+$', "", config[[col]])
    config[[col]] <- trimws(config[[col]])
    config[[col]][config[[col]] == ""] <- NA
  }

  # 2. Load Pipeline Config (config_default.csv)
  # TODO: Migrate to pipeline_config.yaml
  default_path <- file.path(config_dir, "config_default.csv")
  if (!file.exists(default_path)) {
    warning(sprintf("Pipeline config not found: %s. Using empty defaults.", default_path))
    defaults <- character(0)
  } else {
    if (requireNamespace("readr", quietly = TRUE)) {
      default_df <- readr::read_csv(default_path,
        show_col_types = FALSE, comment = "#",
        na = c("", "NA", "na"), trim_ws = TRUE
      )
      default_df <- as.data.frame(default_df)
    } else {
      default_df <- read.csv(default_path,
        stringsAsFactors = FALSE, comment.char = "#",
        na.strings = c("", "NA", "na")
      )
    }
    defaults <- setNames(default_df$value, default_df$parameter)
  }

  # 3. Load Execution Config (Optional)
  execution <- list()
  if (!is.null(execution_config_path) && file.exists(execution_config_path)) {
    # Assuming JSON for execution config
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      execution <- jsonlite::read_json(execution_config_path)

      # Flatten 'steps' if present to allow grouping in JSON
      if ("steps" %in% names(execution)) {
        steps <- execution$steps
        execution$steps <- NULL # Remove the nested structure after flattening
        for (step_name in names(steps)) {
          step_params <- steps[[step_name]]
          for (param in names(step_params)) {
            # Check for collision
            if (param %in% names(execution)) {
              warning(sprintf("Parameter '%s' in step '%s' overrides global parameter", param, step_name))
            }
            execution[[param]] <- step_params[[param]]
          }
        }
      }
    } else {
      warning("jsonlite not installed. Skipping execution config.")
    }
  }

  # 4. Load Methods Config (Legacy/Optional)
  methods_path <- file.path(config_dir, "methods.csv")
  if (!file.exists(methods_path)) {
    methods <- data.frame()
  } else {
    methods <- read.csv(methods_path, stringsAsFactors = FALSE)
  }

  list(config = config, defaults = defaults, execution = execution, methods = methods)
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
  # 1. Try execution config (highest priority)
  if (!is.null(config_list$execution) && param_name %in% names(config_list$execution)) {
    return(config_list$execution[[param_name]])
  }

  # 2. Try config defaults (Pipeline Config)
  if (param_name %in% names(config_list$defaults)) {
    val <- config_list$defaults[[param_name]]
    # Try to convert to appropriate type
    if (val == "TRUE") {
      return(TRUE)
    }
    if (val == "FALSE") {
      return(FALSE)
    }
    # Check for numeric but be careful not to convert strings that look like numbers but should be strings if context implies
    # But for now, simple heuristic:
    if (grepl("^-?\\d+\\.?\\d*$", val)) {
      return(as.numeric(val))
    }
    if (grepl("^\\d+:\\d+$", val)) {
      parts <- strsplit(val, ":")[[1]]
      return(as.numeric(parts[1]):as.numeric(parts[2]))
    }
    return(val)
  }

  # 3. Fallback to provided default
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

  # Add separator if log file already exists (same run_id reused)
  if (file.exists(run_log_path) && file.info(run_log_path)$size > 0) {
    existing_log <- file(run_log_path, open = "a")
    writeLines("", existing_log)
    writeLines("========================================", existing_log)
    writeLines(sprintf("NEW RUN STARTED (same run_id: %s)", run_id), existing_log)
    writeLines("========================================", existing_log)
    close(existing_log)
  }

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
#' If output_dir is provided, uses it directly; otherwise constructs from run_id
#'
#' @param run_id Run identifier
#' @param step Step number or name
#' @param filename Filename
#' @param output_base_dir Base directory for outputs (default: "/data/user3/sobj/pipe")
#' @param output_dir Optional: direct output directory path (overrides run_id-based path)
#'
#' @return Full path to output file
#'
#' @export
get_output_path <- function(run_id, step, filename, output_base_dir = "/data/user3/sobj/pipe", output_dir = NULL) {
  if (!is.null(output_dir)) {
    # Use provided output_dir directly
    step_dir <- file.path(output_dir, sprintf("step%d", step))
  } else {
    # Use run_id-based path
    step_dir <- file.path(output_base_dir, run_id, sprintf("step%d", step))
  }
  dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(step_dir, filename)
}

#' Save Intermediate Object
#'
#' Saves an R object using qs for fast serialization
#' Automatically appends _1, _2, ... if file already exists to prevent overwriting
#'
#' @param obj Object to save
#' @param filepath Full path to output file
#' @param log_list Optional log list for logging
#' @param prevent_overwrite If TRUE, appends _1, _2, ... if file exists (default: TRUE)
#'
#' @return Actual filepath used (may differ if prevent_overwrite is TRUE)
#'
#' @export
save_intermediate <- function(obj, filepath, log_list = NULL, prevent_overwrite = TRUE) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for fast serialization")
  }

  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  # Prevent overwriting by appending _1, _2, ... if file exists
  if (prevent_overwrite && file.exists(filepath)) {
    base_path <- tools::file_path_sans_ext(filepath)
    ext <- tools::file_ext(filepath)
    if (ext != "") {
      ext <- paste0(".", ext)
    }

    counter <- 1
    new_filepath <- paste0(base_path, "_", counter, ext)
    while (file.exists(new_filepath)) {
      counter <- counter + 1
      new_filepath <- paste0(base_path, "_", counter, ext)
    }
    filepath <- new_filepath

    if (!is.null(log_list)) {
      log_message(sprintf("File already exists, saving as: %s", filepath), log_list, level = "WARNING")
    }
  }

  qs::qsave(obj, filepath)

  if (!is.null(log_list)) {
    log_message(sprintf("Saved intermediate object to: %s", filepath), log_list)
  }

  invisible(filepath)
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
load_intermediate <- function(filepath, log_list = NULL, input_dir = NULL) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for fast serialization")
  }

  # If input_dir is provided and filepath is relative, construct full path
  if (!is.null(input_dir)) {
    # Check if filepath is absolute (starts with /)
    if (!startsWith(filepath, "/")) {
      filepath <- file.path(input_dir, filepath)
    }
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
  non_ribo <- c(
    "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6",
    "RPL10A", "RPL10L"
  ) # Add more known exceptions

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
