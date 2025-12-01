#!/usr/bin/env Rscript
# Step 4: Doublet detection using scDblFinder
# Usage: Rscript pipe4_doubletfinder.R --config <config_path> --run_id <run_id> --input_step <step> --output_step <step>

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
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type = "character", default = NULL,
              help = "Path to config.csv file", metavar = "character"),
  make_option(c("--run_id", "-r"), type = "character", default = "run1",
              help = "Run ID", metavar = "character"),
  make_option(c("--input_step", "-i"), type = "integer", default = 3,
              help = "Input step number", metavar = "integer"),
  make_option(c("--output_step", "-o"), type = "integer", default = 4,
              help = "Output step number", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

# Source utility functions
script_dir <- dirname(normalizePath(opt$config))
source(file.path(script_dir, "../myR/R/pipe_utils.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting doublet detection (scDblFinder)", opt$output_step), log_list)

# Load config
config_list <- load_config(opt$config)
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")

# Load input from previous step
input_path <- get_output_path(opt$run_id, opt$input_step, 
                              get_param("output_step3_soupx", config_list, "step3_soupx_list.qs"),
                              output_base_dir)
log_message(sprintf("Loading input from: %s", input_path), log_list)
sl <- load_intermediate(input_path, log_list)

# Create plots directory
plots_dir <- file.path(output_base_dir, opt$run_id, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Process each sample
for (sample_name in names(sl)) {
  log_message(sprintf("Processing %s", sample_name), log_list)
  
  obj <- sl[[sample_name]]
  
  tryCatch({
    # Need to run SCTransform first if not already done
    if (!"SCT" %in% names(obj@assays)) {
      log_message(sprintf("  Running SCTransform for %s (required for scDblFinder)", sample_name), log_list)
      
      method <- get_param("sct_method", config_list, "glmGamPoi")
      vst_flavor <- get_param("sct_vst_flavor", config_list, "v2")
      variable_features_n <- as.numeric(get_param("sct_variable_features_n", config_list, 3000))
      conserve_memory <- get_param("sct_conserve_memory", config_list, TRUE)
      
      obj <- SCTransform(obj, 
                        method = method,
                        vst.flavor = vst_flavor,
                        variable.features.n = variable_features_n,
                        conserve.memory = conserve_memory,
                        verbose = FALSE)
      
      # Run PCA for visualization
      obj <- RunPCA(obj, verbose = FALSE)
    }
    
    log_message(sprintf("  Running scDblFinder for %s", sample_name), log_list)
    
    # Convert to SingleCellExperiment
    sce <- as.SingleCellExperiment(obj)
    
    # Run scDblFinder
    use_clusters <- get_param("scdblfinder_clusters", config_list, TRUE)
    sce <- scDblFinder(sce, clusters = use_clusters)
    
    # Add results back to Seurat object
    obj$scDblFinder.class <- sce$scDblFinder.class
    obj$scDblFinder.score <- sce$scDblFinder.score
    
    # Create plots
    plot_width <- as.numeric(get_param("scdblfinder_plot_width", config_list, 7))
    plot_height <- as.numeric(get_param("scdblfinder_plot_height", config_list, 6))
    
    # Doublet score plot
    if ("pca" %in% names(obj@reductions)) {
      p_dbl_score <- FeaturePlot(obj, features = "scDblFinder.score", reduction = "pca") + 
        ggplot2::labs(title = paste("Doublet Scores -", sample_name))
      ggplot2::ggsave(file.path(plots_dir, sprintf("Doublet_Score_Plot_%s.png", sample_name)), 
                     plot = p_dbl_score, width = plot_width, height = plot_height)
      
      # Doublet classification plot
      p_dbl_class <- DimPlot(obj, group.by = "scDblFinder.class", reduction = "pca") + 
        ggplot2::labs(title = paste("Doublet Classification -", sample_name))
      ggplot2::ggsave(file.path(plots_dir, sprintf("Doublet_Class_Plot_%s.png", sample_name)), 
                     plot = p_dbl_class, width = plot_width, height = plot_height)
    }
    
    # Optional: Remove doublets (default: keep all)
    remove_doublets <- get_param("scdblfinder_remove_doublets", config_list, FALSE)
    if (remove_doublets) {
      n_before <- ncol(obj)
      obj <- subset(obj, subset = scDblFinder.class == "singlet")
      n_after <- ncol(obj)
      log_message(sprintf("  Removed doublets: %d -> %d cells", n_before, n_after), log_list)
    } else {
      log_message(sprintf("  Keeping all cells (including doublets) for %s", sample_name), log_list)
    }
    
    sl[[sample_name]] <- obj
    log_message(sprintf("  Completed %s: %d cells", sample_name, ncol(obj)), log_list)
    
  }, error = function(e) {
    log_message(sprintf("!!! ERROR in step %d processing %s: %s !!!", 
                       opt$output_step, sample_name, e$message), log_list, level = "ERROR")
    # Keep original object if scDblFinder fails
    log_message(sprintf("  Keeping original object for %s", sample_name), log_list, level = "WARNING")
  })
}

# Save results
if (length(sl) > 0) {
  output_path <- get_output_path(opt$run_id, opt$output_step, 
                                 get_param("output_step4_doublet", config_list, "step4_doublet_list.qs"),
                                 output_base_dir)
  
  log_message(sprintf("Saving results to: %s", output_path), log_list)
  save_intermediate(sl, output_path, log_list)
  log_message(sprintf("Step %d completed: %d samples processed", opt$output_step, length(sl)), log_list)
} else {
  log_message("No samples were successfully processed!", log_list, level = "ERROR")
  close_logging(log_list)
  stop("Step 4 failed: No samples processed")
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Processed %d samples.\n", opt$output_step, length(sl)))

