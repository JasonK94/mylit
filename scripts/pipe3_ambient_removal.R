#!/usr/bin/env Rscript
# Step 3: Ambient RNA removal using SoupX
# Usage: Rscript pipe3_ambient_removal.R --config <config_path> --run_id <run_id> --input_step <step> --output_step <step>

# Load lightweight pipeline environment first
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grepl("--file=", cmd_args)]
if (length(file_arg) > 0) {
  script_path <- sub("--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "start_pipe.R"))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SoupX)
  library(ggplot2)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type = "character", default = NULL,
              help = "Path to config.csv file", metavar = "character"),
  make_option(c("--run_id", "-r"), type = "character", default = "run1",
              help = "Run ID", metavar = "character"),
  make_option(c("--input_step", "-i"), type = "integer", default = 2,
              help = "Input step number", metavar = "integer"),
  make_option(c("--output_step", "-o"), type = "integer", default = 3,
              help = "Output step number", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

# Source utility functions
# Use consistent approach: pipe_dir is parent of scripts directory
script_dir_for_source <- dirname(normalizePath(opt$config))
pipe_dir <- dirname(script_dir_for_source)
source(file.path(pipe_dir, "myR/R/pipe_utils.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting SoupX ambient RNA removal", opt$output_step), log_list)

# Load config
config_list <- load_config(opt$config)
config <- config_list$config
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")

# SoupX tuning parameters
soupx_tfidf_start <- as.numeric(get_param("soupx_tfidf_start", config_list, 1))
soupx_soup_quantile_start <- as.numeric(get_param("soupx_soup_quantile_start", config_list, 0.9))
soupx_tfidf_min_floor <- as.numeric(get_param("soupx_tfidf_floor", config_list, 0.3))
soupx_soup_quantile_floor <- as.numeric(get_param("soupx_soup_quantile_floor", config_list, 0.1))
soupx_param_step <- as.numeric(get_param("soupx_param_step", config_list, 0.1))
if (soupx_param_step <= 0) {
  soupx_param_step <- 0.1
}

# Load input from previous step
input_path <- get_output_path(opt$run_id, opt$input_step, 
                              get_param("output_step2_nmz", config_list, "step2_nmz_list.qs"),
                              output_base_dir)
log_message(sprintf("Loading input from: %s", input_path), log_list)
sl <- load_intermediate(input_path, log_list)

# Create plots directory
plots_dir <- file.path(output_base_dir, opt$run_id, "plots", "step3_soupx")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Create log file for plot generation
plot_log_file <- file.path(plots_dir, "plot_generation.log")
if (file.exists(plot_log_file)) {
  file.remove(plot_log_file)
}
write_log_entry <- function(sample_name, plot_type, status, message = "") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s - %s: %s %s\n", 
                       timestamp, sample_name, plot_type, status, message)
  cat(log_entry, file = plot_log_file, append = TRUE)
}

auto_estimate_with_retry <- function(sc, sample_name, plot_path) {
  current_tfidf <- soupx_tfidf_start
  current_soup_q <- soupx_soup_quantile_start
  attempt <- 1
  last_error <- NULL
  
  repeat {
    log_message(sprintf("  autoEstCont attempt %d for %s (tfidfMin=%.2f, soupQuantile=%.2f)", 
                       attempt, sample_name, current_tfidf, current_soup_q), log_list)
    temp_plot_path <- paste0(tools::file_path_sans_ext(plot_path), sprintf("_attempt%d.pdf", attempt))
    success <- FALSE
    sc_result <- NULL
    
    tryCatch({
      pdf(temp_plot_path, width = plot_width, height = plot_height)
      sc_result <- autoEstCont(
        sc,
        doPlot = do_plot,
        forceAccept = force_accept,
        tfidfMin = current_tfidf,
        soupQuantile = current_soup_q
      )
      dev.off()
      file.rename(temp_plot_path, plot_path)
      success <- TRUE
    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
      unlink(temp_plot_path)
      last_error <<- e$message
      write_log_entry(
        sample_name,
        "Rho_Distribution",
        "FAILED",
        sprintf("tfidfMin=%.2f, soupQuantile=%.2f, %s", current_tfidf, current_soup_q, e$message)
      )
      log_message(sprintf("  ERROR: autoEstCont attempt %d failed for %s: %s", 
                          attempt, sample_name, e$message), log_list, level = "ERROR")
    })
    
    if (success) {
      write_log_entry(
        sample_name,
        "Rho_Distribution",
        "SUCCESS",
        sprintf("tfidfMin=%.2f, soupQuantile=%.2f", current_tfidf, current_soup_q)
      )
      if (attempt > 1) {
        log_message(sprintf("  autoEstCont succeeded after lowering thresholds to tfidfMin=%.2f, soupQuantile=%.2f", 
                            current_tfidf, current_soup_q), log_list)
      }
      return(list(sc = sc_result, attempts = attempt, tfidf = current_tfidf, soup = current_soup_q))
    }
    
    if (current_tfidf <= soupx_tfidf_min_floor && current_soup_q <= soupx_soup_quantile_floor) {
      stop(sprintf(
        paste(
          "autoEstCont failed for %s even at minimum thresholds",
          "(tfidfMin=%.2f, soupQuantile=%.2f). Last error: %s",
          "Consider lowering soupx_tfidf_floor or soupx_soup_quantile_floor in config."
        ),
        sample_name,
        current_tfidf,
        current_soup_q,
        ifelse(is.null(last_error), "Unknown error", last_error)
      ))
    }
    
    new_tfidf <- max(soupx_tfidf_min_floor, current_tfidf - soupx_param_step)
    new_soup_q <- max(soupx_soup_quantile_floor, current_soup_q - soupx_param_step)
    
    if (new_tfidf == current_tfidf && new_soup_q == current_soup_q) {
      stop(sprintf(
        paste(
          "autoEstCont failed for %s and parameters can no longer be reduced.",
          "Last error: %s. Adjust soupx_tfidf_start/soupx_tfidf_floor or soupx_soup_quantile_start/soupx_soup_quantile_floor in config."
        ),
        sample_name,
        ifelse(is.null(last_error), "Unknown error", last_error)
      ))
    }
    
    current_tfidf <- new_tfidf
    current_soup_q <- new_soup_q
    attempt <- attempt + 1
  }
}

# Process each sample
for (i in seq_len(nrow(config))) {
  row <- config[i, ]
  sample_name <- row$sample_name
  
  if (!sample_name %in% names(sl)) {
    log_message(sprintf("Skipping %s (not in input)", sample_name), log_list, level = "WARNING")
    next
  }
  
  log_message(sprintf("Processing %s", sample_name), log_list)
  
  obj <- sl[[sample_name]]
  
  tryCatch({
    # Check if raw counts are available
    if (is.na(row$dir_input_raw_barcode_matrix) || row$dir_input_raw_barcode_matrix == "") {
      log_message(sprintf("  No raw counts available for %s, skipping SoupX", sample_name), log_list, level = "WARNING")
      next
    }
    
    log_message(sprintf("  Loading raw counts for %s", sample_name), log_list)
    
    # Load raw counts
    raw_counts_list <- Seurat::Read10X(row$dir_input_raw_barcode_matrix)
    
    # Extract Gene Expression if it's a list
    if (is.list(raw_counts_list)) {
      raw_counts <- raw_counts_list$`Gene Expression`
    } else {
      raw_counts <- raw_counts_list
    }
    
    # Match barcode suffixes if needed (for SNP samples)
    # Use gem_name instead of index for consistency with Step 1
    if (row$multiplex_method == "SNP") {
      gem_name <- row$gem_name
      colnames(raw_counts) <- paste0(colnames(raw_counts), "_", gem_name)
    }
    
    log_message(sprintf("  Running SoupX for %s", sample_name), log_list)
    
    # Get common features
    common_features <- intersect(rownames(raw_counts), rownames(obj))
    if (length(common_features) == 0) {
      stop("No common features between raw and filtered counts")
    }
    
    # Create SoupChannel
    sc <- SoupChannel(raw_counts[common_features, ], 
                     obj@assays$RNA$counts[common_features, ])
    
    # Set clusters
    if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
      stop("seurat_clusters not found in metadata. Run clustering first.")
    }
    sc <- setClusters(sc, obj$seurat_clusters)
    
    # Estimate contamination
    force_accept <- get_param("soupx_force_accept", config_list, TRUE)
    do_plot <- get_param("soupx_do_plot", config_list, TRUE)
    
    # Save SoupX plots
    plot_width <- as.numeric(get_param("soupx_plot_width", config_list, 7))
    plot_height <- as.numeric(get_param("soupx_plot_height", config_list, 6))
    
    rho_plot_path <- file.path(plots_dir, sprintf("SoupX_Rho_Distribution_%s.pdf", sample_name))
    est_result <- auto_estimate_with_retry(sc, sample_name, rho_plot_path)
    sc <- est_result$sc
    log_message(sprintf("  Saved Rho Distribution plot: %s", rho_plot_path), log_list)
    
    # Marker Distribution plot
    marker_plot_path <- file.path(plots_dir, sprintf("SoupX_Marker_Distribution_%s.pdf", sample_name))
    tryCatch({
      pdf(marker_plot_path, width = plot_width, height = plot_height)
      print(plotMarkerDistribution(sc))
      dev.off()
      write_log_entry(sample_name, "Marker_Distribution", "SUCCESS", "")
      log_message(sprintf("  Saved Marker Distribution plot: %s", marker_plot_path), log_list)
    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
      suggestion <- "Try lowering soupx_tfidf_min or soupx_soup_quantile in config if markers are sparse."
      error_msg <- sprintf("Error: %s. %s", e$message, suggestion)
      write_log_entry(sample_name, "Marker_Distribution", "FAILED", error_msg)
      log_message(sprintf("  WARNING: Failed to generate Marker Distribution plot for %s: %s. %s", 
                         sample_name, e$message, suggestion), log_list, level = "WARNING")
      # Continue processing even if plot fails
    })
    
    # Adjust counts with guard for missing contamination fractions
    log_message(sprintf("  Adjusting counts for %s", sample_name), log_list)
    corrected_counts <- tryCatch({
      adjustCounts(sc)
    }, error = function(e) {
      if (grepl("Contamination fractions must have already been calculated/set", e$message, fixed = TRUE)) {
        suggestion <- "autoEstCont produced no contamination fraction. Consider lowering soupx_tfidf_start/floor or soupx_soup_quantile_start/floor."
        log_message(sprintf("  WARNING: %s - %s", e$message, suggestion), log_list, level = "WARNING")
        write_log_entry(sample_name, "AdjustCounts", "SKIPPED", sprintf("%s %s", e$message, suggestion))
        return(NULL)
      }
      stop(e)
    })
    
    if (is.null(corrected_counts)) {
      log_message(sprintf("  Keeping original counts for %s due to missing contamination fraction", sample_name), log_list, level = "WARNING")
      next
    }
    
    # Update object with corrected counts
    # Match cells
    common_cells <- intersect(colnames(corrected_counts), colnames(obj))
    if (length(common_cells) == 0) {
      stop("No common cells between corrected counts and object")
    }
    
    # Create new object with corrected counts
    obj <- Seurat::CreateSeuratObject(
      counts = corrected_counts[, common_cells],
      meta.data = obj@meta.data[common_cells, , drop = FALSE]
    )
    
    # Preserve other assays if they exist
    for (assay_name in names(obj@assays)) {
      if (assay_name != "RNA") {
        obj[[assay_name]] <- obj[[assay_name]]
      }
    }
    
    sl[[sample_name]] <- obj
    log_message(sprintf("  Completed %s: %d cells", sample_name, ncol(obj)), log_list)
    
  }, error = function(e) {
    log_message(sprintf("!!! ERROR in step %d processing %s: %s !!!", 
                       opt$output_step, sample_name, e$message), log_list, level = "ERROR")
    # Keep original object if SoupX fails
    log_message(sprintf("  Keeping original counts for %s", sample_name), log_list, level = "WARNING")
  })
}

# Save results
if (length(sl) > 0) {
  output_path <- get_output_path(opt$run_id, opt$output_step, 
                                 get_param("output_step3_soupx", config_list, "step3_soupx_list.qs"),
                                 output_base_dir)
  
  log_message(sprintf("Saving results to: %s", output_path), log_list)
  save_intermediate(sl, output_path, log_list)
  log_message(sprintf("Step %d completed: %d samples processed", opt$output_step, length(sl)), log_list)
} else {
  log_message("No samples were successfully processed!", log_list, level = "ERROR")
  close_logging(log_list)
  stop("Step 3 failed: No samples processed")
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Processed %d samples.\n", opt$output_step, length(sl)))

