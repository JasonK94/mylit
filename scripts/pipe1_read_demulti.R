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
              help = "Output step number", metavar = "integer"),
  make_option(c("--downsample", "-d"), type = "numeric", default = NULL,
              help = "Downsample fraction (0-1) for testing, e.g., 0.1 for 10%% of cells", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

# Source utility functions
# Use consistent approach: pipe_dir is parent of script_dir
pipe_dir <- dirname(script_dir)

# Source pipeline utility functions
source(file.path(pipe_dir, "myR/R/pipe_utils.R"))
source(file.path(pipe_dir, "myR/R/pipe_demulti.R"))
source(file.path(pipe_dir, "myR/R/utils_demulti.R"))

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

# Check if output already exists and skip if requested
if (file.exists(output_path) && opt$input_step == opt$output_step) {
  log_message(sprintf("Output file already exists: %s. Loading existing results...", output_path), log_list)
  sl <- load_intermediate(output_path, log_list)
  log_message(sprintf("Loaded %d samples from existing file", length(sl)), log_list)
  close_logging(log_list)
  cat(sprintf("Step %d completed: Loaded %d samples from existing file.\n", opt$output_step, length(sl)))
  quit(save = "no", status = 0)
}

# Create plots directory
plots_dir <- file.path(output_base_dir, opt$run_id, "plots", "step1_demultiplex")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Set random seed
seed <- as.numeric(get_param("seed", config_list, 1234))
set.seed(seed)

# Separate SNP and HTO samples
snp_samples <- config[config$multiplex_method == "SNP", ]
hto_samples <- config[config$multiplex_method %in% c("HTO", "CMO"), ]

log_message(sprintf("Found %d SNP samples and %d HTO samples", 
                   nrow(snp_samples), nrow(hto_samples)), log_list)

# Initialize result lists
snp_barcode_mappings <- list()
hto_objects <- list()
sl <- list()

# ============================================================================
# Process SNP samples
# ============================================================================
if (nrow(snp_samples) > 0) {
  log_message("Processing SNP samples...", log_list)
  
  for (i in seq_len(nrow(snp_samples))) {
    row <- snp_samples[i, ]
    sample_name <- row$sample_name
    gem_name <- row$gem_name
    
    log_message(sprintf("Processing SNP sample %d/%d: %s (GEM: %s)", 
                       i, nrow(snp_samples), sample_name, gem_name), log_list)
    
    tryCatch({
      # Get demux file
      demux_file <- as.character(row$dir_demultiplex_output)
      demux_file <- gsub('^["\']+|["\']+$', '', demux_file)  # Remove quotes
      demux_file <- trimws(demux_file)
      
      if (is.na(demux_file) || demux_file == "" || !file.exists(demux_file)) {
        stop(sprintf("dir_demultiplex_output is required for SNP samples but is missing or file not found for %s: %s", 
                     sample_name, row$dir_demultiplex_output))
      }
      
      # Use new demultiplex_demuxalot function
      log_message(sprintf("  Demultiplexing %s using demuxalot", sample_name), log_list)
      barcode_map <- demultiplex_demuxalot(
        demuxalot_posterior = demux_file,
        barcode_col = get_param("demuxalot_barcode_col", config_list, "BARCODE"),
        singlet_threshold = as.numeric(get_param("demuxalot_singlet_threshold", config_list, 0.5)),
        doublet_threshold = as.numeric(get_param("demuxalot_doublet_threshold", config_list, 0.3)),
        gem_name = gem_name,
        gem_col = get_param("metadata_gem_col", config_list, "GEM"),
        return_probs = TRUE
      )
      
      # Add config metadata columns
      config_cols_to_add <- c("name", "patient_name", "contamination_risk", "sample_name", "gem_name", "demultiplex_id")
      for (col in config_cols_to_add) {
        if (col %in% colnames(row) && !is.na(row[[col]])) {
          barcode_map[[col]] <- as.character(row[[col]])
        }
      }
      
      snp_barcode_mappings[[sample_name]] <- barcode_map
      log_message(sprintf("  Completed %s: %d barcodes mapped", sample_name, nrow(barcode_map)), log_list)
      
    }, error = function(e) {
      log_message(sprintf("!!! ERROR processing SNP sample %s: %s !!!", sample_name, e$message), log_list, level = "ERROR")
    })
  }
  
  # Combine all SNP barcode mappings
  if (length(snp_barcode_mappings) > 0) {
    log_message("Combining SNP barcode mappings...", log_list)
    # Remove rownames before binding to avoid conflicts
    for (i in seq_along(snp_barcode_mappings)) {
      rownames(snp_barcode_mappings[[i]]) <- NULL
    }
    all_snp_mappings <- bind_rows(snp_barcode_mappings)
    # Reset rownames to Barcode column (which already has suffix)
    rownames(all_snp_mappings) <- all_snp_mappings$Barcode
    log_message(sprintf("Total SNP barcodes: %d", nrow(all_snp_mappings)), log_list)
    
    # Save SNP mappings separately
    snp_output_path <- get_output_path(opt$run_id, opt$output_step, 
                                      "step1_snp_barcode_mappings.qs",
                                      output_base_dir)
    save_intermediate(all_snp_mappings, snp_output_path, log_list)
  }
}

# ============================================================================
# Process HTO samples
# ============================================================================
if (nrow(hto_samples) > 0) {
  log_message("Processing HTO samples...", log_list)
  
  for (i in seq_len(nrow(hto_samples))) {
    row <- hto_samples[i, ]
    sample_name <- row$sample_name
    gem_name <- row$gem_name
    
    log_message(sprintf("Processing HTO sample %d/%d: %s (GEM: %s)", 
                       i, nrow(hto_samples), sample_name, gem_name), log_list)
    
    processed_obj <- tryCatch({
      # Load filtered counts
      log_message(sprintf("  Loading HTO data for %s", sample_name), log_list)
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
        log_message(sprintf("  Performing %s demultiplexing for %s", demux_method, sample_name), log_list)
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
        log_message(sprintf("  No HTO assay found for %s (tried '%s'), skipping demultiplexing", 
                           sample_name, hto_assay_name), log_list, level = "WARNING")
      }
      
      # Add config metadata columns
      config_cols_to_add <- c("name", "patient_name", "contamination_risk", "sample_name", "gem_name", "demultiplex_id")
      for (col in config_cols_to_add) {
        if (col %in% colnames(row) && !is.na(row[[col]])) {
          obj[[col]] <- as.character(row[[col]])
        }
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
      log_message(sprintf("!!! ERROR in step %d processing HTO sample %s: %s !!!", 
                         opt$output_step, sample_name, e$message), log_list, level = "ERROR")
      log_message(sprintf("!!! Skipping HTO sample %s !!!", sample_name), log_list, level = "ERROR")
      return(NULL)
    })
    
    if (!is.null(processed_obj)) {
      hto_objects[[sample_name]] <- processed_obj
    }
  }
  
  # Save HTO objects separately
  if (length(hto_objects) > 0) {
    log_message("Saving HTO objects...", log_list)
    hto_output_path <- get_output_path(opt$run_id, opt$output_step, 
                                      "step1_hto_list.qs",
                                      output_base_dir)
    save_intermediate(hto_objects, hto_output_path, log_list)
  }
}

# ============================================================================
# Create Seurat objects for SNP samples (using barcode mappings)
# ============================================================================
if (length(snp_barcode_mappings) > 0 && exists("all_snp_mappings")) {
  log_message("Creating Seurat objects for SNP samples...", log_list)
  
  # Group by GEM (same GEM can have multiple samples)
  gem_groups <- split(snp_samples, snp_samples$gem_name)
  
  for (gem_name in names(gem_groups)) {
    gem_samples <- gem_groups[[gem_name]]
    
    for (i in seq_len(nrow(gem_samples))) {
      row <- gem_samples[i, ]
      sample_name <- row$sample_name
      
      tryCatch({
        log_message(sprintf("  Creating Seurat object for %s (GEM: %s)", sample_name, gem_name), log_list)
        
        # Load count matrices
        filtered_counts <- Seurat::Read10X(row$dir_input_filtered_barcode_matrix)
        
        # Get barcode mapping for this sample
        sample_mapping <- all_snp_mappings[all_snp_mappings$sample_name == sample_name, ]
        
        if (nrow(sample_mapping) == 0) {
          stop(sprintf("No barcode mapping found for sample %s", sample_name))
        }
        
        # Match barcodes using original barcode (before suffix)
        # demultiplex_demuxalot already added suffix to Barcode column
        # We need to match original barcodes first, then add suffix to counts
        original_barcodes_in_mapping <- sample_mapping$barcode  # Original barcode (no suffix)
        barcodes_with_suffix <- sample_mapping$Barcode  # Barcode with suffix (rownames)
        
        # Find common original barcodes
        common_original_barcodes <- intersect(colnames(filtered_counts), original_barcodes_in_mapping)
        
        if (length(common_original_barcodes) == 0) {
          # Debug: show some examples
          log_message(sprintf("  DEBUG: Count matrix barcodes (first 5): %s", 
                            paste(head(colnames(filtered_counts), 5), collapse = ", ")), log_list)
          log_message(sprintf("  DEBUG: Mapping barcodes (first 5): %s", 
                            paste(head(original_barcodes_in_mapping, 5), collapse = ", ")), log_list)
          stop(sprintf("No common original barcodes found between count matrix and demux results for %s", sample_name))
        }
        
        # Get corresponding barcodes with suffix from mapping
        mapping_subset <- sample_mapping[sample_mapping$barcode %in% common_original_barcodes, , drop = FALSE]
        
        # Create a named vector mapping original -> suffix barcode
        orig_to_suffix <- setNames(mapping_subset$Barcode, mapping_subset$barcode)
        
        # Update count matrix barcode names: add suffix to match mapping rownames
        count_barcodes_original <- colnames(filtered_counts)
        count_barcodes_to_update <- count_barcodes_original %in% names(orig_to_suffix)
        colnames(filtered_counts)[count_barcodes_to_update] <- orig_to_suffix[count_barcodes_original[count_barcodes_to_update]]
        
        # Now match using barcodes with suffix (rownames of sample_mapping)
        common_cells <- intersect(colnames(filtered_counts), rownames(sample_mapping))
        
        if (length(common_cells) == 0) {
          # Additional debug info
          log_message(sprintf("  DEBUG: After suffix addition, count barcodes (first 5): %s", 
                            paste(head(colnames(filtered_counts), 5), collapse = ", ")), log_list)
          log_message(sprintf("  DEBUG: Mapping rownames (first 5): %s", 
                            paste(head(rownames(sample_mapping), 5), collapse = ", ")), log_list)
          log_message(sprintf("  DEBUG: Common original barcodes: %d", length(common_original_barcodes)), log_list)
          stop(sprintf("No common cells found between count matrix and demux results for %s", sample_name))
        }
        
        # Create Seurat object
        obj <- Seurat::CreateSeuratObject(
          counts = filtered_counts[, common_cells],
          meta.data = sample_mapping[common_cells, , drop = FALSE]
        )
        
        # Add config metadata columns
        config_cols_to_add <- c("name", "patient_name", "contamination_risk", "sample_name", "gem_name", "demultiplex_id")
        for (col in config_cols_to_add) {
          if (col %in% colnames(row) && !is.na(row[[col]])) {
            obj[[col]] <- as.character(row[[col]])
          }
        }
        
        # Extract metadata from filename if patterns are provided
        timepoint_regex <- get_param("metadata_timepoint_regex", config_list, NULL)
        if (!is.null(timepoint_regex) && timepoint_regex != "") {
          timepoint <- str_extract(sample_name, timepoint_regex)
          if (!is.na(timepoint)) {
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
        
        # Downsample if requested (for testing)
        if (!is.null(opt$downsample) && opt$downsample > 0 && opt$downsample < 1) {
          n_cells_before <- ncol(obj)
          n_cells_after <- max(1, round(n_cells_before * opt$downsample))
          if (n_cells_after < n_cells_before) {
            set.seed(seed + i)  # Reproducible sampling
            cells_to_keep <- sample(colnames(obj), n_cells_after)
            obj <- obj[, cells_to_keep]
            log_message(sprintf("  Downsampled %s: %d -> %d cells (%.1f%%)", 
                              sample_name, n_cells_before, n_cells_after, 
                              opt$downsample * 100), log_list)
          }
        }
        
        sl[[sample_name]] <- obj
        log_message(sprintf("  Completed %s: %d cells", sample_name, ncol(obj)), log_list)
        
      }, error = function(e) {
        log_message(sprintf("!!! ERROR creating Seurat object for SNP sample %s: %s !!!", 
                           sample_name, e$message), log_list, level = "ERROR")
      })
    }
  }
}

# Combine SNP and HTO objects
if (length(hto_objects) > 0) {
  sl <- c(sl, hto_objects)
}

# Save results
if (length(sl) > 0) {
  log_message(sprintf("Saving results to: %s", output_path), log_list)
  save_intermediate(sl, output_path, log_list)
  log_message(sprintf("Step %d completed: %d samples processed (%d SNP, %d HTO)", 
                     opt$output_step, length(sl), 
                     length(snp_barcode_mappings), length(hto_objects)), log_list)
} else {
  log_message("No samples were successfully processed!", log_list, level = "ERROR")
  close_logging(log_list)
  stop("Step 1 failed: No samples processed")
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Processed %d samples.\n", opt$output_step, length(sl)))
