#!/usr/bin/env Rscript
# Step 2: Normalization and clustering
# Usage: Rscript pipe2_nmz_clustering.R --config <config_path> --run_id <run_id> --input_step <step> --output_step <step> --nmz <LogNormalize|SCTransform>

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
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type = "character", default = NULL,
              help = "Path to config.csv file", metavar = "character"),
  make_option(c("--run_id", "-r"), type = "character", default = "run1",
              help = "Run ID", metavar = "character"),
  make_option(c("--input_step", "-i"), type = "integer", default = 1,
              help = "Input step number", metavar = "integer"),
  make_option(c("--output_step", "-o"), type = "integer", default = 2,
              help = "Output step number", metavar = "integer"),
  make_option(c("--nmz", "-n"), type = "character", default = "LogNormalize",
              help = "Normalization method: LogNormalize or SCTransform", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

if (!opt$nmz %in% c("LogNormalize", "SCTransform")) {
  stop("--nmz must be either 'LogNormalize' or 'SCTransform'")
}

# Source utility functions
script_dir <- dirname(normalizePath(opt$config))
source(file.path(script_dir, "../myR/R/pipe_utils.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting normalization (%s) and clustering", 
                   opt$output_step, opt$nmz), log_list)

# Load config
config_list <- load_config(opt$config)
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")

# Load input from previous step
input_path <- get_output_path(opt$run_id, opt$input_step, 
                              get_param("output_step1_demulti", config_list, "step1_demulti_list.qs"),
                              output_base_dir)
log_message(sprintf("Loading input from: %s", input_path), log_list)
sl <- load_intermediate(input_path, log_list)

# Process each sample
for (sample_name in names(sl)) {
  log_message(sprintf("Processing %s", sample_name), log_list)
  
  obj <- sl[[sample_name]]
  
  tryCatch({
    # Calculate QC metrics
    log_message(sprintf("  Calculating QC metrics for %s", sample_name), log_list)
    
    # Mitochondrial percentage
    mt_pattern <- get_param("qc_mt_pattern", config_list, "^MT-")
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
    
    # Ribosomal percentage (optional)
    rps_pattern <- get_param("qc_ribo_pattern_rps", config_list, "^RPS")
    rpl_pattern <- get_param("qc_ribo_pattern_rpl", config_list, "^RPL")
    obj <- calculate_ribosomal_percentage(obj, rps_pattern, rpl_pattern)
    
    # Optional QC filtering (if enabled in config)
    min_nfeature <- as.numeric(get_param("qc_min_nfeature", config_list, 200))
    max_nfeature <- as.numeric(get_param("qc_max_nfeature", config_list, 5000))
    max_percent_mt <- as.numeric(get_param("qc_max_percent_mt", config_list, 20))
    min_cells <- as.numeric(get_param("qc_min_cells", config_list, 20))
    
    # Apply QC filtering
    n_before <- ncol(obj)
    obj <- subset(obj, subset = nFeature_RNA > min_nfeature & 
                             nFeature_RNA < max_nfeature & 
                             percent.mt < max_percent_mt)
    n_after <- ncol(obj)
    
    if (n_after < min_cells) {
      stop(sprintf("Fewer than %d cells after QC (%d cells), skipping sample", min_cells, n_after))
    }
    
    log_message(sprintf("  QC filtering: %d -> %d cells", n_before, n_after), log_list)
    
    # Normalization
    if (opt$nmz == "LogNormalize") {
      log_message(sprintf("  Performing LogNormalize for %s", sample_name), log_list)
      
      scale_factor <- as.numeric(get_param("lognorm_scale_factor", config_list, 10000))
      nfeatures <- as.numeric(get_param("lognorm_nfeatures", config_list, 2000))
      selection_method <- get_param("lognorm_selection_method", config_list, "disp")
      
      obj <- NormalizeData(obj, scale.factor = scale_factor, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = selection_method, 
                                  nfeatures = nfeatures, verbose = FALSE)
      obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
      
    } else if (opt$nmz == "SCTransform") {
      log_message(sprintf("  Performing SCTransform for %s (this may take a while)", sample_name), log_list)
      
      method <- get_param("sct_method", config_list, "glmGamPoi")
      vst_flavor <- get_param("sct_vst_flavor", config_list, "v2")
      variable_features_n <- as.numeric(get_param("sct_variable_features_n", config_list, 3000))
      conserve_memory <- get_param("sct_conserve_memory", config_list, TRUE)
      
      obj <- SCTransform(obj, 
                        method = method,
                        vst.flavor = vst_flavor,
                        variable.features.n = variable_features_n,
                        conserve.memory = conserve_memory,
                        verbose = TRUE)
    }
    
    # PCA
    log_message(sprintf("  Running PCA for %s", sample_name), log_list)
    npcs <- as.numeric(get_param("pre_soupx_npcs", config_list, 20))
    obj <- RunPCA(obj, npcs = min(npcs, ncol(obj) - 1), verbose = FALSE)
    
    # Clustering (for SoupX or downstream analysis)
    log_message(sprintf("  Clustering for %s", sample_name), log_list)
    dims_str <- get_param("pre_soupx_dims", config_list, "1:10")
    dims <- eval(parse(text = dims_str))
    resolution <- as.numeric(get_param("pre_soupx_resolution", config_list, 0.5))
    
    obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
    
    sl[[sample_name]] <- obj
    log_message(sprintf("  Completed %s: %d cells", sample_name, ncol(obj)), log_list)
    
  }, error = function(e) {
    log_message(sprintf("!!! ERROR in step %d processing %s: %s !!!", 
                       opt$output_step, sample_name, e$message), log_list, level = "ERROR")
    # Remove failed sample from list
    sl[[sample_name]] <<- NULL
  })
}

# Save results
if (length(sl) > 0) {
  output_filename <- ifelse(opt$nmz == "LogNormalize", 
                           get_param("output_step2_nmz", config_list, "step2_nmz_list.qs"),
                           "step2_sct_list.qs")
  output_path <- get_output_path(opt$run_id, opt$output_step, output_filename, output_base_dir)
  
  log_message(sprintf("Saving results to: %s", output_path), log_list)
  save_intermediate(sl, output_path, log_list)
  log_message(sprintf("Step %d completed: %d samples processed", opt$output_step, length(sl)), log_list)
} else {
  log_message("No samples were successfully processed!", log_list, level = "ERROR")
  close_logging(log_list)
  stop("Step 2 failed: No samples processed")
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Processed %d samples.\n", opt$output_step, length(sl)))

