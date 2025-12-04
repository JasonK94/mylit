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
  make_option(c("--config", "-c"),
    type = "character", default = NULL,
    help = "Path to config.csv file", metavar = "character"
  ),
  make_option(c("--run_id", "-r"),
    type = "character", default = "run1",
    help = "Run ID", metavar = "character"
  ),
  make_option(c("--input_step", "-i"),
    type = "integer", default = 5,
    help = "Input step number", metavar = "integer"
  ),
  make_option(c("--output_step", "-o"),
    type = "integer", default = 6,
    help = "Output step number", metavar = "integer"
  ),
  make_option(c("--method", "-m"),
    type = "character", default = "RPCA",
    help = "Integration method: RPCA or scVI", metavar = "character"
  ),
  make_option(c("--input_dir", "-I"),
    type = "character", default = NULL,
    help = "Input directory (overrides run_id-based path)", metavar = "character"
  ),
  make_option(c("--output_dir", "-O"),
    type = "character", default = NULL,
    help = "Output directory (overrides run_id-based path)", metavar = "character"
  )
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
# Use consistent approach: pipe_dir is parent of scripts directory
script_dir_for_source <- dirname(normalizePath(opt$config))
pipe_dir <- dirname(script_dir_for_source)
source(file.path(pipe_dir, "myR/R/pipe_utils.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting integration (%s)", opt$output_step, opt$method), log_list)

# Load config
config_list <- load_config(opt$config)
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")

# Setup future for parallelization
log_message("Setting up parallel processing...", log_list)
# Check for environment variable first (for resource-limited execution)
# Then check config, then use default
n_workers <- as.numeric(Sys.getenv(
  "INTEGRATION_WORKERS",
  get_param("integration_workers", config_list, 4)
))
# Safety limit: don't use more than 8 workers per process
n_workers <- min(n_workers, 8)
log_message(sprintf("Using %d workers for parallel processing", n_workers), log_list)
future::plan(future::multisession, workers = n_workers)
options(future.globals.maxSize = 200 * 1024^3)

# Load input from previous step (Step 5: Doublet Detection)
# Use input_dir if provided, otherwise use run_id-based path
if (!is.null(opt$input_dir)) {
  input_filename <- get_param("output_step5_doublet", config_list, "step5_doublet_list.qs")
  input_path <- file.path(opt$input_dir, sprintf("step%d", opt$input_step), input_filename)
  log_message(sprintf("Using custom input directory: %s", opt$input_dir), log_list)
} else {
  input_path <- get_output_path(
    opt$run_id, opt$input_step,
    get_param("output_step5_doublet", config_list, "step5_doublet_list.qs"),
    output_base_dir
  )
}
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

# ------------------------------------------------------------------
# Metadata Joining & Enhancement (Robust)
# ------------------------------------------------------------------
log_message("Enhancing metadata...", log_list)

# Load required libraries
if (!requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
if (!requireNamespace("tibble", quietly = TRUE)) library(tibble)
if (!requireNamespace("readr", quietly = TRUE)) library(readr)

# 1. HTO/SNP Mapping Logic (Preserved for ID standardization)
hto_id_map <- c("1_IS_72" = "57",   "2_IS_24" = "58",
                "22_IS_24" = "87", "23_IS_24" = "88",
                "28_IS_24" = "92","4_IS_24" = "64",
                "4_IS_72" = "66", "8_IS_72" = "74",
                "14_IS_24" = "80", "18_IS_24" = "83",
                "2_IS_72" = "60", "10_ICH_24"="77",
                "13_SAH_24"="77", "15_SAH_24"="81",
                "17_SAH_24"="82", "19_ICH_24"="84",
                "24_ICH_24"="89", "3_ICH_72"="63",
                "5_SAH_24"="67", "5_SAH_72"="68",
                "6_SAH_24"="69",  "6_SAH_72"="70",
                "9_ICH_72"="76", "11_ICH_24"="78",
                "20_SAH_24"="85", "21_SAH_24"="86",
                "25_SAH_24"="90", "27_ICH_24"="91",
                "3_ICH_24"="61","3_ICH_48"="62", "7_SAH_72"="72")

hto_gem_map <- c("Samples1_4" = "GEM9", "Samples5_9" = "GEM10",
                 "GEM1_2" = "GEM11", "GEM2_2" = "GEM12")

meta <- merged@meta.data

# Apply mappings
if ("sample_id" %in% colnames(meta)) {
  matches <- meta$sample_id %in% names(hto_id_map)
  if (any(matches)) meta$sample_id[matches] <- hto_id_map[meta$sample_id[matches]]
  
  matches_gem <- meta$GEM %in% names(hto_gem_map)
  if (any(matches_gem)) meta$GEM[matches_gem] <- hto_gem_map[meta$GEM[matches_gem]]
  
  meta$join_key <- meta$sample_id
} else {
  if ("Best_Sample" %in% colnames(meta)) {
    meta$join_key <- meta$Best_Sample
  } else {
    meta$join_key <- meta$sample_name
  }
}
merged@meta.data <- meta

# 2. Join with Clinical Metadata (Robust Method with Conflict Handling)
clinical_meta_path <- get_param("dir_meta_data", config_list, "/home/user3/Clinical information/meta_data_prep_251113.csv")

if (file.exists(clinical_meta_path)) {
  log_message(sprintf("Joining with clinical metadata from: %s", clinical_meta_path), log_list)
  
  meta_clinical <- readr::read_csv(clinical_meta_path, show_col_types = FALSE)
  
  # Prepare Seurat meta
  current_meta <- merged@meta.data
  
  # Ensure join keys are character
  current_meta$join_key <- as.character(current_meta$join_key)
  meta_clinical$sample_no <- as.character(meta_clinical$sample_no)
  
  # Identify overlapping columns (excluding join key)
  dup_cols <- intersect(names(current_meta), names(meta_clinical))
  dup_cols <- setdiff(dup_cols, c("join_key", "sample_no")) # Exclude keys
  
  # Rename duplicates in clinical data to avoid automatic .x/.y suffixing by join
  # We want controlled suffixing: _sobj (original), _csv (new)
  # Actually, left_join does .x and .y. We can rename them afterwards.
  
  # Join
  # We use join_key (Seurat) = sample_no (Clinical)
  new_meta <- current_meta %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(meta_clinical, by = c("join_key" = "sample_no"), suffix = c("_sobj", "_csv")) %>%
    tibble::column_to_rownames("barcodes")
  
  # Robust Conflict Handling
  # Iterate over duplicated columns (now with suffixes)
  for (col in dup_cols) {
    col_sobj <- paste0(col, "_sobj")
    col_csv <- paste0(col, "_csv")
    
    if (col_sobj %in% names(new_meta) && col_csv %in% names(new_meta)) {
      val_sobj <- as.character(new_meta[[col_sobj]])
      val_csv <- as.character(new_meta[[col_csv]])
      
      # 1. If identical (ignoring NA), keep one (sobj preferred or csv? usually csv is ground truth)
      # 2. If one is NA and other is valid, keep valid.
      # 3. If both valid and different, keep both.
      
      # Simple logic:
      # If identical -> remove _csv, rename _sobj back to col
      # If different -> keep both
      
      # Check identity (handling NAs)
      is_identical <- identical(val_sobj, val_csv)
      # Or more robust:
      # is_identical <- all(val_sobj == val_csv | (is.na(val_sobj) & is.na(val_csv)), na.rm = TRUE) 
      # But identical() handles NAs well.
      
      if (is_identical) {
        # Keep sobj version, remove csv version
        new_meta[[col]] <- new_meta[[col_sobj]]
        new_meta[[col_sobj]] <- NULL
        new_meta[[col_csv]] <- NULL
      } else {
        # Conflict: Keep original as 'col' and new as 'col_csv'
        log_message(sprintf("Conflict or difference in column '%s'. Keeping original as '%s' and new as '%s_csv'.", col, col, col), log_list)
        new_meta[[col]] <- new_meta[[col_sobj]]
        new_meta[[col_sobj]] <- NULL
        # col_csv remains as is
      }
    }
  }
  
  merged@meta.data <- new_meta
  log_message("Clinical metadata joined successfully.", log_list)
  
} else {
  log_message("WARNING: Clinical metadata file not found.", log_list, level = "WARNING")
}

# 3. Add Metadata from Config (e.g., SET, SET_DETAIL)
log_message("Adding metadata from config...", log_list)
if (!is.null(config_list$config)) {
  config_samples <- config_list$config
  
  # Columns to add (exclude standard ones)
  std_cols <- c("no", "name", "patient_name", "contamination_risk", "sample_name", "gem_name", 
                "demultiplex_id", "sample_col_meta_data", "gem_col_meta_data", "dir_meta_data",
                "dir_input_filtered_barcode_matrix", "dir_input_raw_barcode_matrix", "dir_input_bam",
                "dir_demultiplex_output", "multiplex_method", "demultiplex_method", "demultiplex_annotation_col",
                "output_step1_demulti", "output_step2_nmz", "output_step3_soupx", "output_step4_sct", 
                "output_step5_doublet", "output_step6_integration_rpca", "output_step6_integration_scvi")
  
  cols_to_add <- setdiff(names(config_samples), std_cols)
  
  if (length(cols_to_add) > 0) {
    log_message(sprintf("Adding columns from config: %s", paste(cols_to_add, collapse = ", ")), log_list)
    
    # Ensure sample_name is character for joining
    config_samples$sample_name <- as.character(config_samples$sample_name)
    merged$sample_name <- as.character(merged$sample_name) # Ensure this exists
    
    # If sample_name not in merged, try orig.ident
    if (!"sample_name" %in% colnames(merged@meta.data)) {
       merged$sample_name <- merged$orig.ident
    }
    
    # Join
    current_meta <- merged@meta.data
    current_meta <- current_meta %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(config_samples %>% dplyr::select(sample_name, dplyr::all_of(cols_to_add)), by = "sample_name") %>%
      tibble::column_to_rownames("barcodes")
      
    merged@meta.data <- current_meta
  }
}

# 4. Create SET columns (Fallback if not in config)
if (!"SET" %in% colnames(merged@meta.data) && "GEM" %in% colnames(merged@meta.data)) {
  log_message("Creating SET column from GEM (Fallback)...", log_list)
  merged$SET <- dplyr::case_when(
    merged$GEM %in% c("GEM1", "GEM2", "GEM3", "GEM4") ~ "SET1",
    merged$GEM %in% c("GEM5", "GEM6", "GEM7", "GEM8") ~ "SET2",
    merged$GEM %in% c("GEM9", "GEM10", "GEM11", "GEM12") ~ "SET3",
    TRUE ~ "Other"
  )
}

log_message(sprintf("SET distribution: %s", paste(names(table(merged$SET)), table(merged$SET), sep = "=", collapse = ", ")), log_list)

# ------------------------------------------------------------------
# Optional: Remove Doublets
# ------------------------------------------------------------------
remove_doublets <- get_param("integration_remove_doublets", config_list, FALSE)
if (remove_doublets) {
  if ("scDblFinder.class" %in% colnames(merged@meta.data)) {
    log_message("Removing doublets...", log_list)
    n_before <- ncol(merged)
    merged <- subset(merged, subset = scDblFinder.class == "singlet")
    n_after <- ncol(merged)
    log_message(sprintf("Removed %d doublets. Cells: %d -> %d", n_before - n_after, n_before, n_after), log_list)
  } else {
    log_message("WARNING: scDblFinder.class not found. Skipping doublet removal.", log_list, level = "WARNING")
  }
}

# Set default assay
DefaultAssay(merged) <- "RNA"

# ------------------------------------------------------------------
# Integration Logic (Simplified & Robust)
# ------------------------------------------------------------------

if (opt$method == "RPCA") {
  log_message("Performing RPCA integration (Default: SCT + FindIntegrationAnchors)...", log_list)

  # Check for GEM column
  batch_col <- get_param("metadata_gem_col", config_list, "GEM")
  if (!batch_col %in% colnames(merged@meta.data)) {
    stop(sprintf("Batch column '%s' not found in metadata", batch_col))
  }

  # Split by batch
  log_message(sprintf("Splitting by batch column: %s", batch_col), log_list)
  obj.list <- SplitObject(merged, split.by = batch_col)

  # Filter out small batches
  min_cells_per_batch <- 50
  obj.list <- obj.list[vapply(obj.list, function(x) as.integer(ncol(x)), integer(1)) >= min_cells_per_batch]

  if (length(obj.list) < 2) {
    stop(sprintf("Need at least 2 batches for integration, but only %d batch(es) with >= %d cells", length(obj.list), min_cells_per_batch))
  }

  # SCTransform each batch
  npcs <- as.numeric(get_param("rpca_npcs", config_list, 50))
  nfeatures <- as.numeric(get_param("rpca_nfeatures", config_list, 3000))
  
  obj.list <- lapply(obj.list, function(x) {
    DefaultAssay(x) <- "RNA"
    x <- SCTransform(x, assay = "RNA", new.assay.name = "SCT", method = "glmGamPoi", vst.flavor = "v2", variable.features.n = nfeatures, conserve.memory = TRUE, verbose = FALSE)
    x <- RunPCA(x, assay = "SCT", layer = "scale.data", npcs = npcs, verbose = FALSE)
    x
  })

  # Select features & Prep
  features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nfeatures)
  obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features, verbose = FALSE)

  # Find Anchors & Integrate
  log_message("Finding integration anchors (RPCA)...", log_list)
  anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:npcs, verbose = TRUE)
  
  log_message("Integrating data...", log_list)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  DefaultAssay(integrated) <- "integrated"
  
  # Downstream
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
  integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
  integrated <- FindClusters(integrated, resolution = 0.6, cluster.name = "clusters_rpca", verbose = FALSE)
  integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)

} else if (opt$method == "scVI") {
  log_message("Performing scVI integration (Default: IntegrateLayers)...", log_list)

  if (!requireNamespace("SeuratWrappers", quietly = TRUE) || !requireNamespace("reticulate", quietly = TRUE)) {
    stop("SeuratWrappers and reticulate are required for scVI integration")
  }

  # Setup Python
  conda_env <- get_param("scvi_conda_env", config_list, "/home/user3/miniconda3/envs/scvi-env/")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.unsetenv("PYTHONHOME")
  Sys.unsetenv("PYTHONPATH")
  reticulate::use_condaenv("scvi-env", required = TRUE)
  
  DefaultAssay(merged) <- "RNA"
  
  # Pre-processing for scVI (Normalize -> Scale -> PCA)
  # scVI uses raw counts, but Seurat workflow often initializes with PCA
  log_message("Running pre-processing (Normalize, Scale, PCA)...", log_list)
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = 3000, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, verbose = FALSE)

  # Batch Column
  batch_col <- get_param("scvi_batch_col", config_list, "GEM")
  if (!batch_col %in% colnames(merged@meta.data)) {
     stop(sprintf("Batch column '%s' not found. Please check config.", batch_col))
  }
  
  # Ensure counts layer exists (scVI needs it)
  # If v5, JoinLayers might be needed if split
  if (packageVersion("Seurat") >= "5.0.0") {
      merged <- JoinLayers(merged, assay = "RNA")
  }

  # Run scVI
  log_message(sprintf("Running scVIIntegration with batch='%s'...", batch_col), log_list)
  
  # No fallback logic. If it fails, let it fail with a clear error.
  merged <- IntegrateLayers(
    object = merged,
    method = SeuratWrappers::scVIIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.scvi",
    batch = batch_col,
    layers = "counts",
    conda_env = conda_env,
    verbose = TRUE
  )

  # Downstream
  merged <- RunUMAP(merged, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
}
  # Check for loadings (scVI typically doesn't create loadings, only embeddings)
  if ("integrated.scvi" %in% names(merged@reductions)) {
    if (is.null(merged@reductions$integrated.scvi@feature.loadings)) {
      log_message("Note: integrated.scvi does not have feature loadings (this is normal for scVI)", log_list)
    } else {
      log_message(sprintf(
        "integrated.scvi loadings: %d features, %d dimensions",
        nrow(merged@reductions$integrated.scvi@feature.loadings),
        ncol(merged@reductions$integrated.scvi@feature.loadings)
      ), log_list)
    }
  }

  # Intermediate save after integration (before downstream analysis)
  output_path_scvi <- get_output_path(
    opt$run_id, opt$output_step,
    get_param("output_step6_integration_scvi", config_list, "step6_integration_scvi.qs"),
    output_base_dir
  )
  log_message(sprintf("scVIIntegration completed. Saving intermediate result to: %s", output_path_scvi), log_list)
  qs::qsave(merged, output_path_scvi)

  # Downstream analysis
  log_message("Running downstream analysis...", log_list)
  integrated <- tryCatch(
    {
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
        log_message(sprintf(
          "Warning: Requested dims %d exceeds available dimensions %d. Using 1:%d",
          max(dims), max_dims, max_dims
        ), log_list, level = "WARNING")
        dims <- 1:max_dims
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
    log_message(sprintf(
      "Warning: Requested dims %d exceeds available dimensions %d. Using 1:%d",
      max(dims), max_dims, max_dims
    ), log_list, level = "WARNING")
    dims <- 1:max_dims
  }

  resolution <- as.numeric(get_param("scvi_resolution", config_list, 0.6))

  log_message(sprintf("Finding neighbors using dims: %s", paste(dims, collapse = ", ")), log_list)
  # Downstream
  merged <- FindNeighbors(merged, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
  merged <- FindClusters(merged, resolution = 0.6, verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)

  # Join layers after integration
  if (packageVersion("Seurat") >= "5.0.0") {
    merged <- JoinLayers(merged)
  }

  integrated <- merged # Assign merged to integrated for consistency with RPCA path

  log_message("Downstream analysis completed", log_list)
}

# Save final results
output_filename <- ifelse(opt$method == "RPCA",
  get_param("output_step6_integration_rpca", config_list, "step6_integration_rpca.qs"),
  get_param("output_step6_integration_scvi", config_list, "step6_integration_scvi.qs")
)
output_path <- get_output_path(opt$run_id, opt$output_step, output_filename, output_base_dir, opt$output_dir)

# Check if intermediate save path is different from final save path
if (opt$method == "scVI" && exists("output_path_scvi")) {
  if (output_path != output_path_scvi) {
    log_message(sprintf("Final save path differs from intermediate. Final: %s", output_path), log_list)
  }
}

log_message(sprintf("Saving final integrated object to: %s", output_path), log_list)
save_intermediate(integrated, output_path, log_list)

# Final summary
log_message(sprintf(
  "Step %d completed: Integration finished with %d cells",
  opt$output_step, ncol(integrated)
), log_list)
if ("integrated.scvi" %in% names(integrated@reductions) || "integrated.rpca" %in% names(integrated@reductions)) {
  reduction_name <- ifelse(opt$method == "RPCA", "integrated.rpca", "integrated.scvi")
  reduction_embeddings <- Embeddings(integrated, reduction = reduction_name)
  log_message(sprintf(
    "Final reduction '%s': %d cells, %d dimensions",
    reduction_name, nrow(reduction_embeddings), ncol(reduction_embeddings)
  ), log_list)
}

close_logging(log_list)
cat(sprintf(
  "Step %d completed successfully. Integrated object has %d cells.\n",
  opt$output_step, ncol(integrated)
))
