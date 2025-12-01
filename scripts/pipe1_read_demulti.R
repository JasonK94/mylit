#!/usr/bin/env Rscript
# Step 1: Read data and perform demultiplexing
# Usage: Rscript pipe1_read_demulti.R --config <config_path> --run_id <run_id> [--input_step <step>] [--output_step <step>]

# Load lightweight pipeline environment FIRST (before any library calls)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grepl("--file=", cmd_args)]
if (length(file_arg) > 0) {
  script_path <- sub("--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  # Fallback: try to find script directory
  script_dir <- tryCatch({
    dirname(normalizePath(sys.frame(1)$ofile))
  }, error = function(e) {
    getwd()
  })
}
source(file.path(script_dir, "start_pipe.R"))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(stringr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type = "character", default = NULL,
              help = "Path to config.csv file", metavar = "character"),
  make_option(c("--run_id", "-r"), type = "character", default = "run1",
              help = "Run ID", metavar = "character"),
  make_option(c("--input_step", "-i"), type = "integer", default = 0,
              help = "Input step number (0 = start from scratch)", metavar = "integer"),
  make_option(c("--output_step", "-o"), type = "integer", default = 1,
              help = "Output step number", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

# Source utility functions
pipe_dir <- dirname(script_dir)

# Source pipeline utility functions
source(file.path(pipe_dir, "myR/R/pipe_utils.R"))
source(file.path(pipe_dir, "myR/R/pipe_demulti.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting read and demultiplexing", opt$output_step), log_list)

# Load config
config_list <- load_config(opt$config)
config <- config_list$config
defaults <- config_list$defaults

# Get output directory
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")
output_path <- get_output_path(opt$run_id, opt$output_step, 
                               get_param("output_step1_demulti", config_list, "step1_demulti_list.qs"),
                               output_base_dir)

# Set random seed
seed <- as.numeric(get_param("seed", config_list, 1234))
set.seed(seed)

# Initialize result list
sl <- list()

# Process each sample
for (i in seq_len(nrow(config))) {
  row <- config[i, ]
  sample_name <- row$sample_name
  gem_name <- row$gem_name
  multiplex_method <- row$multiplex_method
  
  log_message(sprintf("Processing sample %d/%d: %s (GEM: %s, Method: %s)", 
                     i, nrow(config), sample_name, gem_name, multiplex_method), log_list)
  
  processed_obj <- tryCatch({
    
    if (multiplex_method == "SNP") {
      # SNP-based demultiplexing
      log_message(sprintf("  Step 1.1: Loading demuxalot output for %s", sample_name), log_list)
      
      # Check if demux file exists
      if (is.na(row$dir_demultiplex_output) || row$dir_demultiplex_output == "") {
        stop(sprintf("dir_demultiplex_output is required for SNP samples but is missing for %s", sample_name))
      }
      
      # Load demux file to get column names (sample names from demultiplex_id)
      demux_data <- read.csv(row$dir_demultiplex_output, nrows = 1, stringsAsFactors = FALSE, check.names = FALSE)
      barcode_col <- get_param("demuxalot_barcode_col", config_list, "BARCODE")
      doublet_separator <- get_param("demuxalot_doublet_separator", config_list, "+")
      
      # Get sample names from demux file columns (excluding barcode column)
      # These should match demultiplex_id values from config
      all_cols <- colnames(demux_data)
      all_sample_cols <- setdiff(all_cols, barcode_col)
      
      # Filter to only singlet samples (exclude doublet combinations which contain separator)
      sample_names <- all_sample_cols[!grepl(doublet_separator, all_sample_cols, fixed = TRUE)]
      
      if (length(sample_names) == 0) {
        stop(sprintf("No singlet sample columns found in demux file (expected columns without '%s' separator): %s", 
                     doublet_separator, row$dir_demultiplex_output))
      }
      
      log_message(sprintf("  Found %d singlet sample(s) in demux file: %s", 
                         length(sample_names), paste(sample_names, collapse = ", ")), log_list)
      
      # Process demuxalot output
      # Note: sample_names are now the actual column names from demux file
      # get_barcode_mapping will use these column names to create Best_Sample
      barcode_map <- process_demuxalot(
        demux_file = row$dir_demultiplex_output,
        sample_names = sample_names,
        barcode_col = barcode_col,
        singlet_threshold = as.numeric(get_param("demuxalot_singlet_threshold", config_list, 0.5)),
        doublet_threshold = as.numeric(get_param("demuxalot_doublet_threshold", config_list, 0.3)),
        gem_name = gem_name,
        gem_index = i
      )
      
      log_message(sprintf("  Step 1.2: Loading count matrices for %s", sample_name), log_list)
      
      # Load count matrices
      filtered_counts <- Seurat::Read10X(row$dir_input_filtered_barcode_matrix)
      
      # Create Seurat object with demultiplexing metadata
      # Match barcodes (need to add GEM suffix to barcodes in counts)
      colnames(filtered_counts) <- paste0(colnames(filtered_counts), "_", i)
      
      # Get cells that are in both counts and barcode_map
      common_cells <- intersect(colnames(filtered_counts), rownames(barcode_map))
      
      if (length(common_cells) == 0) {
        stop(sprintf("No common cells found between count matrix and demux results for %s", sample_name))
      }
      
      # Create Seurat object
      obj <- Seurat::CreateSeuratObject(
        counts = filtered_counts[, common_cells],
        meta.data = barcode_map[common_cells, , drop = FALSE]
      )
      
      # Add GEM metadata
      obj$GEM <- gem_name
      
    } else if (multiplex_method %in% c("HTO", "CMO")) {
      # HTO-based demultiplexing
      log_message(sprintf("  Step 1.1: Loading HTO data for %s", sample_name), log_list)
      
      # Load filtered counts
      filtered_counts_list <- Seurat::Read10X(row$dir_input_filtered_barcode_matrix)
      
      # Extract Gene Expression matrix
      if (is.list(filtered_counts_list)) {
        filtered_counts <- filtered_counts_list$`Gene Expression`
        hto_counts <- filtered_counts_list[[get_param("hto_assay_name", config_list, "Multiplexing Capture")]]
      } else {
        filtered_counts <- filtered_counts_list
        hto_counts <- NULL
      }
      
      # Create initial Seurat object
      obj <- Seurat::CreateSeuratObject(counts = filtered_counts)
      
      # Add HTO assay if available
      hto_assay_name <- get_param("hto_assay_name", config_list, "Multiplexing Capture")
      if (!is.null(hto_counts)) {
        obj[[hto_assay_name]] <- Seurat::CreateAssayObject(counts = hto_counts)
        # Seurat converts spaces to dots in assay names, so check for both
        actual_assay_name <- ifelse(hto_assay_name %in% names(obj@assays), 
                                   hto_assay_name,
                                   gsub(" ", ".", hto_assay_name))
        log_message(sprintf("  Added HTO assay '%s' (stored as '%s') with %d features", 
                           hto_assay_name, actual_assay_name, nrow(hto_counts)), log_list)
      } else {
        actual_assay_name <- NULL
      }
      
      # Perform HTO demultiplexing
      demux_method <- ifelse(is.na(row$demultiplex_method) || row$demultiplex_method == "",
                            get_param("hto_demux_method", config_list, "HTODemux"),
                            row$demultiplex_method)
      
      # Check if HTO assay exists in object (try both original and dot-separated names)
      if (!is.null(actual_assay_name) && actual_assay_name %in% names(obj@assays)) {
        log_message(sprintf("  Step 1.2: Performing %s demultiplexing for %s", demux_method, sample_name), log_list)
        obj <- process_hto_demux(
          obj = obj,
          method = demux_method,
          assay_name = actual_assay_name,  # Use the actual stored name
          positive_quantile = as.numeric(get_param("hto_positive_quantile", config_list, 0.99)),
          sample_name = sample_name,
          gem_name = gem_name
        )
      } else {
        # No HTO assay found, just add metadata
        obj$sample_id <- sample_name
        obj$GEM <- gem_name
        obj$droplet_demulti <- "unknown"
        log_message(sprintf("  Step 1.2: No HTO assay found for %s (tried '%s'), skipping demultiplexing", 
                           sample_name, hto_assay_name), log_list, level = "WARNING")
      }
      
    } else {
      stop(sprintf("Unknown multiplex_method: %s for sample %s", multiplex_method, sample_name))
    }
    
    # Extract metadata from filename if patterns are provided
    timepoint_regex <- get_param("metadata_timepoint_regex", config_list, NULL)
    if (!is.null(timepoint_regex) && timepoint_regex != "") {
      timepoint <- str_extract(sample_name, timepoint_regex)
      if (!is.na(timepoint)) {
        # Map timepoint to day if mapping exists
        day_map_24 <- get_param("metadata_day_map_24", config_list, NULL)
        day_map_72 <- get_param("metadata_day_map_72", config_list, NULL)
        
        if (timepoint == "24" && !is.null(day_map_24)) {
          obj$day <- as.numeric(day_map_24)
        } else if (timepoint == "72" && !is.null(day_map_72)) {
          obj$day <- as.numeric(day_map_72)
        }
        obj$time_point <- timepoint
      }
    }
    
    log_message(sprintf("  Completed %s: %d cells", sample_name, ncol(obj)), log_list)
    obj
    
  }, error = function(e) {
    log_message(sprintf("!!! ERROR in step %d processing %s: %s !!!", 
                       opt$output_step, sample_name, e$message), log_list, level = "ERROR")
    log_message(sprintf("!!! Skipping sample %s !!!", sample_name), log_list, level = "ERROR")
    return(NULL)
  })
  
  if (!is.null(processed_obj)) {
    sl[[sample_name]] <- processed_obj
  }
}

# Save results
if (length(sl) > 0) {
  log_message(sprintf("Saving results to: %s", output_path), log_list)
  save_intermediate(sl, output_path, log_list)
  log_message(sprintf("Step %d completed: %d samples processed", opt$output_step, length(sl)), log_list)
} else {
  log_message("No samples were successfully processed!", log_list, level = "ERROR")
  close_logging(log_list)
  stop("Step 1 failed: No samples processed")
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Processed %d samples.\n", opt$output_step, length(sl)))

