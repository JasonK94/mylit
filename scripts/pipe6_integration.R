#!/usr/bin/env Rscript
# Step 6: Integration using RPCA or scVI
# Usage: Rscript pipe6_integration.R --config <config_path> --run_id <run_id> --input_step <step> --output_step <step> --method <RPCA|scVI>

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
  library(future)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type = "character", default = NULL,
              help = "Path to config.csv file", metavar = "character"),
  make_option(c("--run_id", "-r"), type = "character", default = "run1",
              help = "Run ID", metavar = "character"),
  make_option(c("--input_step", "-i"), type = "integer", default = 5,
              help = "Input step number", metavar = "integer"),
  make_option(c("--output_step", "-o"), type = "integer", default = 6,
              help = "Output step number", metavar = "integer"),
  make_option(c("--method", "-m"), type = "character", default = "RPCA",
              help = "Integration method: RPCA or scVI", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

if (!opt$method %in% c("RPCA", "scVI")) {
  stop("--method must be either 'RPCA' or 'scVI'")
}

# Source utility functions
script_dir <- dirname(normalizePath(opt$config))
source(file.path(script_dir, "../myR/R/pipe_utils.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting integration (%s)", opt$output_step, opt$method), log_list)

# Load config
config_list <- load_config(opt$config)
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")

# Setup future for parallelization
log_message("Setting up parallel processing...", log_list)
# Use more workers for faster processing (Rscript allows more freedom than RStudio)
# 4-6 workers is a good balance between speed and overhead
n_workers <- as.numeric(get_param("integration_workers", config_list, 6))
log_message(sprintf("Using %d workers for parallel processing", n_workers), log_list)
future::plan(future::multisession, workers = n_workers)
options(future.globals.maxSize = 200 * 1024^3)

# Load input from previous step (Step 5: Doublet Detection)
input_path <- get_output_path(opt$run_id, opt$input_step, 
                              get_param("output_step5_doublet", config_list, "step5_doublet_list.qs"),
                              output_base_dir)
log_message(sprintf("Loading input from: %s", input_path), log_list)
sl <- load_intermediate(input_path, log_list)

# Merge all samples into one object
log_message("Merging all samples...", log_list)
log_message(sprintf("Merging %d samples: %s", length(sl), paste(names(sl), collapse = ", ")), log_list)
for (i in seq_along(sl)) {
  log_message(sprintf("  Sample %d: %s has %d cells", i, names(sl)[i], ncol(sl[[i]])), log_list)
}
merged <- merge(sl[[1]], y = sl[-1], add.cell.ids = names(sl))
log_message(sprintf("Merged object has %d cells", ncol(merged)), log_list)

# Set default assay
DefaultAssay(merged) <- "RNA"

if (opt$method == "RPCA") {
  log_message("Performing RPCA integration using FindIntegrationAnchors/IntegrateData method...", log_list)
  
  # Join layers if Seurat v5
  if (packageVersion("Seurat") >= "5.0.0") {
    merged <- JoinLayers(merged, assay = "RNA")
    merged <- DietSeurat(merged, assays = "RNA", dimreducs = NULL, graphs = NULL)
  }
  
  # Check for GEM column
  batch_col <- get_param("metadata_gem_col", config_list, "GEM")
  if (!batch_col %in% colnames(merged@meta.data)) {
    stop(sprintf("Batch column '%s' not found in metadata", batch_col))
  }
  
  # Split by batch
  log_message(sprintf("Splitting by batch column: %s", batch_col), log_list)
  obj.list <- SplitObject(merged, split.by = batch_col)
  
  # Filter out batches with too few cells
  min_cells_per_batch <- 50
  obj.list <- obj.list[vapply(obj.list, function(x) as.integer(ncol(x)), integer(1)) >= min_cells_per_batch]
  
  if (length(obj.list) < 2) {
    stop(sprintf("Need at least 2 batches for integration, but only %d batch(es) with >= %d cells", 
                 length(obj.list), min_cells_per_batch))
  }
  
  log_message(sprintf("Processing %d batches...", length(obj.list)), log_list)
  
  # SCTransform each batch
  npcs <- as.numeric(get_param("rpca_npcs", config_list, 50))
  nfeatures <- as.numeric(get_param("rpca_nfeatures", config_list, 3000))
  method <- get_param("sct_method", config_list, "glmGamPoi")
  vst_flavor <- get_param("sct_vst_flavor", config_list, "v2")
  conserve_memory <- get_param("sct_conserve_memory", config_list, TRUE)
  
  obj.list <- lapply(obj.list, function(x) {
    DefaultAssay(x) <- "RNA"
    x <- SCTransform(x, 
                    assay = "RNA",
                    new.assay.name = "SCT",
                    method = method,
                    vst.flavor = vst_flavor,
                    variable.features.n = nfeatures,
                    conserve.memory = conserve_memory,
                    verbose = FALSE)
    x <- RunPCA(x, assay = "SCT", layer = "scale.data", npcs = npcs, verbose = FALSE)
    x
  })
  
  # Select integration features
  log_message("Selecting integration features...", log_list)
  features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nfeatures)
  
  # PrepSCTIntegration
  log_message("Preparing SCT integration...", log_list)
  obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features, verbose = FALSE)
  
  # Find integration anchors
  log_message("Finding integration anchors...", log_list)
  dims_str <- get_param("rpca_dims", config_list, "1:50")
  # Parse dims string (e.g., "1:50" -> c(1:50))
  if (grepl(":", dims_str)) {
    dims <- eval(parse(text = dims_str))
  } else {
    dims <- as.integer(strsplit(dims_str, ",")[[1]])
  }
  
  anchors <- FindIntegrationAnchors(
    object.list = obj.list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca",
    dims = dims,
    verbose = TRUE
  )
  
  # Integrate data
  log_message("Integrating data...", log_list)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  DefaultAssay(integrated) <- "integrated"
  
  # Downstream analysis
  log_message("Running downstream analysis...", log_list)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
  
  dims_str <- get_param("rpca_dims", config_list, "1:30")
  # Parse dims string (e.g., "1:30" -> c(1:30))
  if (grepl(":", dims_str)) {
    dims <- eval(parse(text = dims_str))
  } else {
    dims <- as.integer(strsplit(dims_str, ",")[[1]])
  }
  
  resolution <- as.numeric(get_param("rpca_resolution", config_list, 0.6))
  
  integrated <- FindNeighbors(integrated, dims = dims, verbose = FALSE)
  integrated <- FindClusters(integrated, resolution = resolution, 
                             cluster.name = "clusters_rpca", verbose = FALSE)
  integrated <- RunUMAP(integrated, dims = dims, verbose = FALSE)
  
} else if (opt$method == "scVI") {
  log_message("Performing scVI integration...", log_list)
  
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    stop("SeuratWrappers package is required for scVI integration")
  }
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package is required for scVI integration")
  }
  
  # Setup Python environment
  conda_env <- get_param("scvi_conda_env", config_list, "/home/user3/miniconda3/envs/scvi-env/")
  python_path <- get_param("scvi_python_path", config_list, "/home/user3/miniconda3/envs/scvi-env/bin/python")
  
  # Unset Python environment variables
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.unsetenv("PYTHONHOME")
  Sys.unsetenv("PYTHONPATH")
  
  # Use conda environment
  reticulate::use_condaenv("scvi-env", required = TRUE)
  
  # Verify Python setup
  py_config <- reticulate::py_config()
  log_message(sprintf("Python: %s", py_config$python), log_list)
  
  # Join layers if Seurat v5 (scVI doesn't need split layers)
  if (packageVersion("Seurat") >= "5.0.0") {
    merged <- JoinLayers(merged, assay = "RNA")
  }
  
  DefaultAssay(merged) <- "RNA"
  
  # Check for batch column
  batch_col <- get_param("scvi_batch_col", config_list, "GEM")
  if (!batch_col %in% colnames(merged@meta.data)) {
    stop(sprintf("Batch column '%s' not found in metadata", batch_col))
  }
  
  # Run PCA first (required for scVI integration)
  log_message("Running PCA (required for scVI)...", log_list)
  npcs <- as.numeric(get_param("scvi_npcs", config_list, 50))
  nfeatures <- as.numeric(get_param("scvi_nfeatures", config_list, 3000))
  
  # Normalize and find variable features if not already done
  if (!"data" %in% Layers(merged, assay = "RNA")) {
    merged <- NormalizeData(merged, verbose = FALSE)
  }
  if (length(VariableFeatures(merged)) == 0) {
    merged <- FindVariableFeatures(merged, nfeatures = nfeatures, verbose = FALSE)
  }
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, npcs = npcs, verbose = FALSE)
  
  # Check PCA results
  if ("pca" %in% names(merged@reductions)) {
    pca_embeddings <- Embeddings(merged, reduction = "pca")
    log_message(sprintf("PCA completed: %d cells, %d dimensions", 
                       nrow(pca_embeddings), ncol(pca_embeddings)), log_list)
  } else {
    stop("PCA failed: pca reduction not found")
  }
  
  # Run scVI integration
  log_message("Running scVIIntegration...", log_list)
  # Note: scVI doesn't need split layers, uses batch column directly
  # IntegrateLayers is in Seurat package (v5), scVIIntegration is in SeuratWrappers
  merged <- IntegrateLayers(
    object = merged,
    method = SeuratWrappers::scVIIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.scvi",
    batch = batch_col,
    layers = "counts",
    conda_env = conda_env,
    python_path = python_path,
    verbose = TRUE
  )
  # Check cell counts before and after integration
  n_cells_before <- ncol(merged)
  log_message(sprintf("Cell count before integration: %d", n_cells_before), log_list)
  
  # Run scVI integration
  log_message("Running scVIIntegration...", log_list)
  merged <- tryCatch({
    IntegrateLayers(
      object = merged,
      method = SeuratWrappers::scVIIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.scvi",
      batch = batch_col,
      layers = "counts",
      conda_env = conda_env,
      python_path = python_path,
      verbose = TRUE
    )
  }, error = function(e) {
    log_message(sprintf("!!! ERROR in scVIIntegration: %s !!!", e$message), log_list, level = "ERROR")
    stop("scVIIntegration failed: ", e$message)
  })
  
  # Check cell counts after integration
  n_cells_after <- ncol(merged)
  log_message(sprintf("Cell count after integration: %d", n_cells_after), log_list)
  
  # Check integrated.scvi reduction
  if (!"integrated.scvi" %in% names(merged@reductions)) {
    stop("integrated.scvi reduction not found after integration")
  }
  
  scvi_embeddings <- Embeddings(merged, reduction = "integrated.scvi")
  log_message(sprintf("integrated.scvi reduction: %d cells, %d dimensions", 
                     nrow(scvi_embeddings), ncol(scvi_embeddings)), log_list)
  
  # Check for loadings (scVI typically doesn't create loadings, only embeddings)
  if ("integrated.scvi" %in% names(merged@reductions)) {
    if (is.null(merged@reductions$integrated.scvi@feature.loadings)) {
      log_message("Note: integrated.scvi does not have feature loadings (this is normal for scVI)", log_list)
    } else {
      log_message(sprintf("integrated.scvi loadings: %d features, %d dimensions", 
                         nrow(merged@reductions$integrated.scvi@feature.loadings),
                         ncol(merged@reductions$integrated.scvi@feature.loadings)), log_list)
    }
  }
  
  # Intermediate save after integration (before downstream analysis)
  output_path_scvi <- get_output_path(opt$run_id, opt$output_step,
                                      get_param("output_step6_integration_scvi", config_list, "step6_integration_scvi.qs"),
                                      output_base_dir)
  log_message(sprintf("scVIIntegration completed. Saving intermediate result to: %s", output_path_scvi), log_list)
  qs::qsave(merged, output_path_scvi)
  
  # Downstream analysis
  log_message("Running downstream analysis...", log_list)
  integrated <- tryCatch({
    dims_str <- get_param("scvi_dims", config_list, "1:30")
    # Ensure dims_str is a single character string
    if (length(dims_str) > 1) {
      dims_str <- dims_str[1]
      warning("Multiple values for scvi_dims, using first: ", dims_str)
    }
    dims_str <- as.character(dims_str)
    # Parse dims string (e.g., "1:30" -> c(1:30))
    if (grepl(":", dims_str)) {
      dims <- eval(parse(text = dims_str))
    } else {
      dims <- as.integer(strsplit(dims_str, ",")[[1]])
    }
    
    # Check if dims are within available dimensions
    max_dims <- ncol(Embeddings(merged, reduction = "integrated.scvi"))
    if (max(dims) > max_dims) {
      log_message(sprintf("Warning: Requested dims %d exceeds available dimensions %d. Using 1:%d", 
                         max(dims), max_dims, max_dims), log_list, level = "WARNING")
      dims <- 1:max_dims
    }
    
    resolution <- as.numeric(get_param("scvi_resolution", config_list, 0.6))
    
    log_message(sprintf("Finding neighbors using dims: %s", paste(dims, collapse = ", ")), log_list)
    merged <- FindNeighbors(merged, reduction = "integrated.scvi", dims = dims, 
                            graph.name = c("scvi.nn", "scvi.snn"), verbose = FALSE)
    
    log_message("Finding clusters...", log_list)
    merged <- FindClusters(merged, resolution = resolution, 
                           graph.name = "scvi.snn", verbose = FALSE)
    
    log_message("Running UMAP...", log_list)
    merged <- RunUMAP(merged, reduction = "integrated.scvi", dims = dims, 
                     reduction.name = "umap.scvi", verbose = FALSE)
    
    # Join layers after integration
    if (packageVersion("Seurat") >= "5.0.0") {
      merged <- JoinLayers(merged)
    }
    
    merged
  }, error = function(e) {
    log_message(sprintf("!!! ERROR in downstream analysis: %s !!!", e$message), log_list, level = "ERROR")
    log_message("Returning object with integration only (without downstream analysis)", log_list, level = "WARNING")
    # Return merged object even if downstream analysis fails
    return(merged)
  })
  
  log_message("Downstream analysis completed", log_list)
}

# Save final results
output_filename <- ifelse(opt$method == "RPCA",
                          get_param("output_step6_integration_rpca", config_list, "step6_integration_rpca.qs"),
                          get_param("output_step6_integration_scvi", config_list, "step6_integration_scvi.qs"))
output_path <- get_output_path(opt$run_id, opt$output_step, output_filename, output_base_dir)

# Check if intermediate save path is different from final save path
if (opt$method == "scVI" && exists("output_path_scvi")) {
  if (output_path != output_path_scvi) {
    log_message(sprintf("Final save path differs from intermediate. Final: %s", output_path), log_list)
  }
}

log_message(sprintf("Saving final integrated object to: %s", output_path), log_list)
save_intermediate(integrated, output_path, log_list)

# Final summary
log_message(sprintf("Step %d completed: Integration finished with %d cells", 
                   opt$output_step, ncol(integrated)), log_list)
if ("integrated.scvi" %in% names(integrated@reductions) || "integrated.rpca" %in% names(integrated@reductions)) {
  reduction_name <- ifelse(opt$method == "RPCA", "integrated.rpca", "integrated.scvi")
  reduction_embeddings <- Embeddings(integrated, reduction = reduction_name)
  log_message(sprintf("Final reduction '%s': %d cells, %d dimensions", 
                     reduction_name, nrow(reduction_embeddings), ncol(reduction_embeddings)), log_list)
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Integrated object has %d cells.\n", 
           opt$output_step, ncol(integrated)))

